"""Score EVERY possible peptide-end window in the secretome, streaming.

``10_7_pepend_windows.R`` builds a training set whose ``all`` split holds only
the MOTIF-anchored candidate windows -- dibasic, amidation, precursor termini,
plus a random sample, 41,870 of them for the C terminus. A known peptide end
that sits at no such motif (25 of the 86 C ends) can therefore never be ranked
by it, and neither can a novel end anywhere else. This module scores every
mature residue of every precursor instead: 2,980,852 C anchors over 5,616
precursors.

That set cannot travel through the R array contract -- as a list of 36x26
tibbles it is tens of GB -- so the windows are built here, a chunk of
precursors at a time, straight out of the per-precursor residue feature
matrices, scored, and thrown away. Only the scalar window score survives, as a
mean and an sd over the ensemble members.

Window geometry is ``lf_pepend_slice`` in R/pepend_windows.R, re-derived rather
than reused: window position ``p`` (1-based) of the window anchored at
precursor residue ``a`` is precursor residue ``a + p - anchor``, where
``anchor`` is 28 for a C window and 8 for an N window; residues outside the
mature precursor are all-zero with the ``padding`` channel set to 1.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import PEPEND_ANCHOR, Config
from .io import MEMBERS_DIR, member_paths


@dataclass
class Residues:
    """Per-precursor residue features, concatenated into one array.

    ``feat`` is ``(total_residues, n_channels)``; row ``offset[i] + r - 1`` is
    residue ``r`` (1-based, precursor coordinates) of precursor ``i``.
    ``n_prot``/``c_prot`` are the mature range, 1-based inclusive, as
    9.2_add_contact_data.R defines it -- everything before ``n_prot`` is the
    signal peptide and counts as outside.
    """

    feat: np.ndarray
    offset: np.ndarray
    n_prot: np.ndarray
    c_prot: np.ndarray
    accession: np.ndarray | None = None

    @classmethod
    def from_npz(cls, path) -> "Residues":
        with np.load(Path(path).expanduser(), allow_pickle=False) as z:
            acc = z["accession"] if "accession" in z.files else None
            return cls(
                feat=np.ascontiguousarray(z["feat"], dtype="float32"),
                offset=z["offset"].astype("int64"),
                n_prot=z["n_prot"].astype("int64"),
                c_prot=z["c_prot"].astype("int64"),
                accession=acc,
            )

    def __len__(self) -> int:
        return int(self.offset.size)

    @property
    def n_anchors(self) -> np.ndarray:
        """Anchors per precursor: every mature residue is a possible end."""
        return self.c_prot - self.n_prot + 1

    def validate(self, cfg: Config) -> "Residues":
        if self.feat.ndim != 2 or self.feat.shape[1] != cfg.n_channels:
            raise ValueError(
                f"feat is {self.feat.shape}, expected (total, {cfg.n_channels})")
        for name, a in (("offset", self.offset), ("n_prot", self.n_prot),
                        ("c_prot", self.c_prot)):
            if a.shape != (len(self),):
                raise ValueError(f"{name} is {a.shape}, expected {(len(self),)}")
        if (self.n_anchors <= 0).any():
            bad = int(np.flatnonzero(self.n_anchors <= 0)[0])
            raise ValueError(f"precursor {bad} has an empty mature range "
                             f"({self.n_prot[bad]}..{self.c_prot[bad]})")
        # every row a window can reach must exist
        if int((self.offset + self.c_prot).max()) > self.feat.shape[0]:
            raise ValueError("offset + c_prot runs past the end of feat")
        return self


def build_windows(res: Residues, i: int, anchors: np.ndarray, cfg: Config,
                  anchor_pos: int, pad_col: int) -> np.ndarray:
    """The ``(len(anchors), seq_len, n_channels)`` windows of precursor ``i``.

    Mirrors ``lf_pepend_slice``: a position outside ``[n_prot, c_prot]`` is
    every channel 0 with ``padding`` 1.
    """
    a = np.asarray(anchors, dtype="int64")
    p = np.arange(1, cfg.seq_len + 1, dtype="int64")
    coord = a[:, None] + p[None, :] - anchor_pos            # (A, seq_len)
    lo, hi = int(res.n_prot[i]), int(res.c_prot[i])
    ok = (coord >= lo) & (coord <= hi)
    x = np.zeros((a.size, cfg.seq_len, cfg.n_channels), dtype="float32")
    rows = int(res.offset[i]) + coord - 1
    x[ok] = res.feat[rows[ok]]
    x[~ok, pad_col] = 1.0
    return x


def load_member_models(out_dir, cfg: Config, term: str, n_seeds: int | None = None):
    """Rebuild every ensemble member from the weights the isolated trainer left.

    ``python -m lf_dcnn train --isolated`` writes
    ``<out_dir>/members/member_NNN_<term>.weights.h5`` per member (cli.py), so
    the whole ensemble can score new windows -- not just member 0, which is all
    R's ``lf_dcnn_run`` rebuilds.
    """
    from .model import build_model

    d = Path(out_dir).expanduser() / MEMBERS_DIR
    ws = sorted(d.glob(f"member_*_{term}.weights.h5"))
    if not ws:
        raise FileNotFoundError(f"no member weights for terminus {term} under {d}")
    if n_seeds is not None and len(ws) != n_seeds:
        raise ValueError(f"found {len(ws)} member weight files, expected {n_seeds}")
    models = []
    for w in ws:
        m = build_model(cfg)
        m.load_weights(w)
        models.append(m)
    return models


def scan(res: Residues, models, cfg: Config, term: str,
         chunk_windows: int = 200_000, batch_size: int = 4096,
         top_k: int = 5, thresholds=(0.2,),
         progress=None) -> dict[str, np.ndarray]:
    """Score every mature position of every precursor.

    Windows are built once per chunk and scored by every member in turn, so the
    feature slicing is paid once rather than ``n_seeds`` times.

    Per window it returns the ensemble mean, the across-member sd (``ddof=1``,
    matching :attr:`Result.pred_sd`), the mean of the ``top_k`` members that
    score THAT window highest, and how many members exceed each of
    ``thresholds`` -- the same summaries the candidate-window tables carry, so
    the exhaustive scan and those tables can be read in one currency.

    The per-member scores of a chunk are held together (``n_members`` x chunk)
    to compute the order statistics, which a running sum/sum-of-squares cannot
    give. That is ~16 MB at the default chunk, against 20 x 3M floats for the
    whole scan.
    """
    if term not in PEPEND_ANCHOR:
        raise ValueError(f"term must be one of {sorted(PEPEND_ANCHOR)}, got {term!r}")
    res.validate(cfg)
    anchor_pos = PEPEND_ANCHOR[term]
    pad_col = cfg.channel_names.index("padding")
    n_mem = len(models)
    kk = max(1, min(int(top_k), n_mem))
    ths = [float(t) for t in thresholds]

    counts = res.n_anchors
    total = int(counts.sum())
    starts = np.concatenate(([0], np.cumsum(counts)))[:-1]

    mean_all = np.zeros(total, dtype="float32")
    sd_all = np.zeros(total, dtype="float32")
    top_all = np.zeros(total, dtype="float32")
    n_gt = {t: np.zeros(total, dtype="int16") for t in ths}
    # The masked 3-way insertion head, when the members were built with one.
    # Accumulated as running sums per chunk rather than held per member: the
    # full (n_members, chunk, 3) block buys nothing here, since no order
    # statistic is wanted -- only the mean, the spread of `inserting`, and how
    # many members call each class.
    ins_names = [str(c) for c in (getattr(cfg, "ins_class_names", None) or [])]
    ins_mean = ins_sd_ins = ins_votes = None
    ins_idx = ins_names.index("inserting") if "inserting" in ins_names else 0
    prot_idx = np.zeros(total, dtype="int32")
    anchor_of = np.zeros(total, dtype="int32")
    for i in range(len(res)):
        s, n = int(starts[i]), int(counts[i])
        prot_idx[s:s + n] = i
        anchor_of[s:s + n] = np.arange(int(res.n_prot[i]), int(res.c_prot[i]) + 1)

    done = 0
    i = 0
    while i < len(res):
        # take precursors until the chunk is big enough; always at least one,
        # so a precursor longer than chunk_windows still goes through
        j, k = i, 0
        while j < len(res) and (k == 0 or k + counts[j] <= chunk_windows):
            k += int(counts[j])
            j += 1
        lo_row, hi_row = int(starts[i]), int(starts[i]) + k

        x = np.empty((k, cfg.seq_len, cfg.n_channels), dtype="float32")
        at = 0
        for p in range(i, j):
            a = anchor_of[int(starts[p]):int(starts[p]) + int(counts[p])]
            x[at:at + a.size] = build_windows(res, p, a, cfg, anchor_pos, pad_col)
            at += int(a.size)

        G = np.empty((n_mem, k), dtype="float64")
        c_sum = c_sq = c_votes = None
        for mi, m in enumerate(models):
            pred = m.predict(x, verbose=0, batch_size=batch_size)
            G[mi] = np.asarray(pred["global"], dtype="float64").reshape(-1)
            if isinstance(pred, dict) and "ins_class" in pred:
                P = np.asarray(pred["ins_class"], dtype="float64")
                if c_sum is None:
                    c_sum = np.zeros((k, P.shape[1]))
                    c_sq = np.zeros(k)
                    c_votes = np.zeros((k, P.shape[1]), dtype="int16")
                c_sum += P
                c_sq += P[:, ins_idx] ** 2
                c_votes[np.arange(k), P.argmax(1)] += 1
        if c_sum is not None:
            if ins_mean is None:
                ins_mean = np.zeros((total, c_sum.shape[1]), dtype="float32")
                ins_sd_ins = np.zeros(total, dtype="float32")
                ins_votes = np.zeros((total, c_sum.shape[1]), dtype="int16")
            mu = c_sum / n_mem
            ins_mean[lo_row:hi_row] = mu
            if n_mem > 1:
                # sample sd of p(inserting) from the running sums, same ddof=1
                # as sd_all; clipped because catastrophic cancellation can push
                # an all-but-identical set a hair below zero.
                var = (c_sq - n_mem * mu[:, ins_idx] ** 2) / (n_mem - 1)
                ins_sd_ins[lo_row:hi_row] = np.sqrt(np.clip(var, 0.0, None))
            ins_votes[lo_row:hi_row] = c_votes
        mean_all[lo_row:hi_row] = G.mean(0)
        if n_mem > 1:
            sd_all[lo_row:hi_row] = G.std(0, ddof=1)
        # the kk highest members of each window, not the kk best members overall
        top_all[lo_row:hi_row] = (G if kk == n_mem else
                                  np.partition(G, n_mem - kk, axis=0)[n_mem - kk:]).mean(0)
        for t in ths:
            n_gt[t][lo_row:hi_row] = (G > t).sum(0)

        done += k
        if progress is not None:
            progress(done, total)
        i = j

    out = {
        "score": mean_all,
        "sd": sd_all,
        f"score_top{kk}": top_all,
        "prot_idx": prot_idx,
        "anchor": anchor_of,
    }
    for t in ths:
        out[f"n_seeds_gt_{t:g}"] = n_gt[t]
    if ins_mean is not None:
        for c, name in enumerate(ins_names):
            out[f"ins_p_{name}"] = ins_mean[:, c]
            out[f"ins_votes_{name}"] = ins_votes[:, c]
        out["ins_sd_inserting"] = ins_sd_ins
    return out


def scan_to_npz(residues_npz, out_dir, term, out_npz, n_seeds=None,
                calibrator=None, chunk_windows=200_000, batch_size=4096,
                top_k=5, thresholds=(0.2,), verbose=True) -> dict[str, np.ndarray]:
    """Read residues, rebuild the ensemble, scan, write ``out_npz``.

    ``out_dir`` is the isolated trainer's output directory; its sibling
    ``../in/config.json`` holds the Config the members were trained with, so the
    scan cannot silently use a different geometry than the training did.
    """
    out_dir = Path(out_dir).expanduser()
    cfg_p = out_dir / "config.json"
    if not cfg_p.exists():
        cfg_p = out_dir.parent / "in" / "config.json"
    cfg = Config.from_json(cfg_p)

    res = Residues.from_npz(residues_npz)
    models = load_member_models(out_dir, cfg, term, n_seeds=n_seeds)
    if verbose:
        print(f"scan: {len(res)} precursors, {int(res.n_anchors.sum()):,} {term} "
              f"anchors, {len(models)} members, config {cfg_p}", flush=True)

    def prog(done, total):
        if verbose:
            print(f"  {done:,}/{total:,} ({100 * done / total:.1f}%)", flush=True)

    got = scan(res, models, cfg, term, chunk_windows=chunk_windows,
               batch_size=batch_size, top_k=top_k, thresholds=thresholds,
               progress=prog if verbose else None)
    if calibrator is not None:
        got["score_cal"] = np.asarray(
            calibrator.predict(got["score"].astype("float64")), dtype="float32")
    extra = {}
    if res.accession is not None:
        extra["accession"] = res.accession
    np.savez_compressed(Path(out_npz).expanduser(), term=np.array(term), **got, **extra)
    if verbose:
        print(f"scan: wrote {out_npz}", flush=True)
    return got
