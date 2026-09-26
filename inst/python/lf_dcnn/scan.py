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
         progress=None) -> dict[str, np.ndarray]:
    """Score every mature position of every precursor.

    Windows are built once per chunk and scored by every member in turn, so the
    feature slicing is paid once rather than ``n_seeds`` times. Returns the
    ensemble mean and the across-member sd (``ddof=1``, matching
    :attr:`Result.pred_sd`), plus the precursor index and anchor of each row.
    """
    if term not in PEPEND_ANCHOR:
        raise ValueError(f"term must be one of {sorted(PEPEND_ANCHOR)}, got {term!r}")
    res.validate(cfg)
    anchor_pos = PEPEND_ANCHOR[term]
    pad_col = cfg.channel_names.index("padding")
    n_mem = len(models)

    counts = res.n_anchors
    total = int(counts.sum())
    starts = np.concatenate(([0], np.cumsum(counts)))[:-1]

    score_sum = np.zeros(total, dtype="float64")
    score_sq = np.zeros(total, dtype="float64")
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

        for m in models:
            g = np.asarray(m.predict(x, verbose=0, batch_size=batch_size)["global"],
                           dtype="float64").reshape(-1)
            score_sum[lo_row:hi_row] += g
            score_sq[lo_row:hi_row] += g * g

        done += k
        if progress is not None:
            progress(done, total)
        i = j

    mean = score_sum / n_mem
    if n_mem > 1:
        var = (score_sq - n_mem * mean * mean) / (n_mem - 1)
        sd = np.sqrt(np.maximum(var, 0.0))
    else:
        sd = np.zeros_like(mean)
    return {
        "score": mean.astype("float32"),
        "sd": sd.astype("float32"),
        "prot_idx": prot_idx,
        "anchor": anchor_of,
    }


def scan_to_npz(residues_npz, out_dir, term, out_npz, n_seeds=None,
                calibrator=None, chunk_windows=200_000, batch_size=4096,
                verbose=True) -> dict[str, np.ndarray]:
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
               batch_size=batch_size, progress=prog if verbose else None)
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
