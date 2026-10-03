"""Per-residue features and labels for the window-free residue models.

The window model slices a 36-residue frame and hands it to a CNN. These models
take ONE residue at a time: the same 26 channels, optionally with the immediate
neighbours bolted on as extra columns, and nothing else. No frame, no
convolution, no pooling.

Row order is the row order of ``lf_pepend_residues.npz`` throughout -- row
``offset[i] + r - 1`` is residue r of precursor i -- which is also the order
``10_9p_export_residue_labels.R`` writes its labels in, so features and labels
line up by index and never need a join.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: the R port's "c" variants append the neighbouring residues' features; ±2 is
#: what ``10_score_AA.R`` uses for the model it actually trains (nn4c).
DEFAULT_CONTEXT = (-2, -1, 1, 2)


@dataclass
class Residues:
    x: np.ndarray            # (n_rows, n_feat) float32
    feat_names: list[str]
    y: dict[str, np.ndarray]  # label name -> (n_rows,) int8
    prot_idx: np.ndarray     # (n_rows,) int32, which precursor each row is in
    resno: np.ndarray        # (n_rows,) int32, 1-based residue number
    mature: np.ndarray       # (n_rows,) int8
    labelled: np.ndarray     # (n_prot,) int8, precursor has any docked peptide
    accession: np.ndarray
    gene: np.ndarray

    def __len__(self) -> int:
        return int(self.x.shape[0])


def _shift_within(a: np.ndarray, prot_idx: np.ndarray, k: int) -> np.ndarray:
    """Shift rows by k WITHIN each precursor, zero-filling across the boundary.

    A plain np.roll would carry the last residues of one precursor into the
    first of the next, inventing neighbours that do not exist -- which is the
    sort of thing that silently inflates a score.
    """
    out = np.zeros_like(a)
    if k == 0:
        return a.copy()
    if k > 0:
        out[k:] = a[:-k]
        bad = prot_idx[k:] != prot_idx[:-k]
        out[k:][bad] = 0.0
    else:
        out[:k] = a[-k:]
        bad = prot_idx[:k] != prot_idx[-k:]
        out[:k][bad] = 0.0
    return out


def load(residues_npz, labels_npz, context=DEFAULT_CONTEXT,
         drop_constant: bool = True) -> Residues:
    """Read features and labels, optionally widening each row with neighbours."""
    with np.load(Path(residues_npz).expanduser(), allow_pickle=False) as z:
        feat = np.ascontiguousarray(z["feat"], dtype="float32")
        offset = z["offset"].astype("int64")
    with np.load(Path(labels_npz).expanduser(), allow_pickle=True) as z:
        lab = {k: z[k] for k in z.files}

    prot_idx = lab["prot_idx"].astype("int32")
    if feat.shape[0] != prot_idx.shape[0]:
        raise ValueError(
            f"features have {feat.shape[0]} rows but labels have {prot_idx.shape[0]}; "
            "they must be the same residue table")

    names = [f"c{i}" for i in range(feat.shape[1])]
    chan = lab.get("channel_names")
    if chan is not None and len(chan) == feat.shape[1]:
        names = [str(c) for c in chan]

    cols = [feat]
    out_names = list(names)
    for k in context:
        cols.append(_shift_within(feat, prot_idx, k))
        tag = f"lag{-k}" if k < 0 else f"lead{k}"
        out_names += [f"{n}_{tag}" for n in names]
    x = np.concatenate(cols, axis=1) if len(cols) > 1 else cols[0]

    if drop_constant:
        keep = x.std(axis=0) > 0
        if not keep.all():
            dropped = [n for n, k in zip(out_names, keep) if not k]
            print(f"dropping {len(dropped)} constant feature(s): "
                  f"{', '.join(dropped[:6])}{' ...' if len(dropped) > 6 else ''}",
                  flush=True)
            x = x[:, keep]
            out_names = [n for n, k in zip(out_names, keep) if k]

    return Residues(
        x=x, feat_names=out_names,
        y={"pep_pocket": lab["y_pocket"].astype("int8"),
           "ct_context": lab["y_ct"].astype("int8")},
        prot_idx=prot_idx, resno=lab["resno"].astype("int32"),
        mature=lab["mature"].astype("int8"), labelled=lab["labelled"].astype("int8"),
        accession=np.asarray(lab["accession"]), gene=np.asarray(lab["gene"]),
    )


def training_rows(d: Residues, n_ctrl: int = 100, seed: int = 42):
    """Rows to train on: every labelled precursor plus `n_ctrl` random controls.

    All 2.98M mature residues would be 0.017% positive. ``10_score_AA.R`` does
    the same thing -- the known genes plus 100 sampled controls -- which lifts
    the rate to a few percent and keeps the controls honest negatives rather
    than residues of a protein whose peptides simply were not docked.
    """
    rng = np.random.default_rng(seed)
    lab_prot = np.flatnonzero(d.labelled == 1)
    ctrl_pool = np.flatnonzero(d.labelled == 0)
    ctrl = rng.choice(ctrl_pool, size=min(n_ctrl, ctrl_pool.size), replace=False)
    keep_prot = np.concatenate([lab_prot, ctrl])
    rows = np.flatnonzero(np.isin(d.prot_idx, keep_prot) & (d.mature == 1))
    return rows, lab_prot, ctrl


def split_by_precursor(d: Residues, rows: np.ndarray, val_frac: float = 0.4,
                       seed: int = 42):
    """Train/val split on WHOLE precursors.

    Splitting on rows would put residues of one protein on both sides, and
    neighbouring residues are near-duplicates of each other once context columns
    are attached -- the val score would then be measuring memorisation.
    """
    rng = np.random.default_rng(seed)
    prots = np.unique(d.prot_idx[rows])
    rng.shuffle(prots)
    n_val = max(1, int(round(val_frac * prots.size)))
    val_prot = set(prots[:n_val].tolist())
    is_val = np.array([p in val_prot for p in d.prot_idx[rows]])
    return rows[~is_val], rows[is_val]


def cache_matrix(d: Residues, path) -> dict:
    """Write the built feature matrix once, for memory-mapped reuse.

    An isolated trainer pays the feature build in every subprocess otherwise --
    reading 3M x 26 and widening it to 125 columns, ~30 s a time, 80 times. Saved
    as a plain .npy so each member can ``np.load(..., mmap_mode="r")`` and touch
    only the rows it needs.
    """
    from pathlib import Path
    import json
    p = Path(path).expanduser()
    p.parent.mkdir(parents=True, exist_ok=True)
    np.save(p, d.x)
    meta = {"feat_names": d.feat_names, "shape": list(d.x.shape)}
    p.with_suffix(".json").write_text(json.dumps(meta))
    return meta
