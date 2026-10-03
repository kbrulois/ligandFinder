"""Write a residues npz with the PLM channels appended, for scanning a PLM trunk.

`lf_pepend_residues.npz` carries the 26 hand-built channels. A trunk trained with
--plm expects 58, and `build_windows` slices straight out of that table, so the
exhaustive scan needs a table of the same width. This joins the per-residue PLM
parquet onto it by (accession, residue index), in the table's own row order.

Proteins with no embedding get zero channels, the same policy 10_8 applies at
training time, and the count is printed rather than assumed.

    python -m lf_winxgb.plm_residues --plm ~/AF2_analysis/lf_plm/esm_c_pca32.parquet
"""
from __future__ import annotations

import argparse
import os

import numpy as np


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lf_winxgb.plm_residues")
    ap.add_argument("--residues-npz", default="~/AF2_analysis/lf_pepend_residues.npz")
    ap.add_argument("--labels-npz", default="~/AF2_analysis/lf_resid_labels.npz",
                    help="for its accession and resno arrays, in residue row order")
    ap.add_argument("--plm", default="~/AF2_analysis/lf_plm/esm_c_pca32.parquet")
    ap.add_argument("--out", default="~/AF2_analysis/lf_pepend_residues_esmc.npz")
    a = ap.parse_args(argv)

    with np.load(os.path.expanduser(a.residues_npz), allow_pickle=False) as z:
        feat = np.ascontiguousarray(z["feat"], dtype="float32")
        offset, n_prot, c_prot = z["offset"], z["n_prot"], z["c_prot"]
    lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
    acc = np.asarray([str(v) for v in lz["accession"]])
    prot_idx, resno = lz["prot_idx"].astype("int64"), lz["resno"].astype("int64")
    print(f"residues {feat.shape[0]:,} x {feat.shape[1]} channels", flush=True)

    import pandas as pd
    df = pd.read_parquet(os.path.expanduser(a.plm))
    cols = [c for c in df.columns if c.startswith("plm_")]
    print(f"plm table {len(df):,} rows x {len(cols)} channels", flush=True)

    # one integer key per (protein, residue) on both sides -- a string join over
    # 3M rows is far slower and no safer here
    code = {s: i for i, s in enumerate(acc)}
    df["pi"] = df["accession"].map(code)
    df = df[df["pi"].notna()]
    key_plm = (df["pi"].to_numpy(dtype="int64") << 20) + df["index"].to_numpy(dtype="int64")
    key_res = (prot_idx << 20) + resno

    order = np.argsort(key_plm, kind="stable")
    kp, vp = key_plm[order], df[cols].to_numpy(dtype="float32")[order]
    pos = np.searchsorted(kp, key_res)
    pos_ok = np.clip(pos, 0, kp.size - 1)
    hit = kp[pos_ok] == key_res

    plm = np.zeros((feat.shape[0], len(cols)), dtype="float32")
    plm[hit] = vp[pos_ok[hit]]
    n_prot_cov = np.unique(prot_idx[hit]).size
    print(f"matched {hit.sum():,} of {hit.size:,} residues "
          f"({n_prot_cov} of {acc.size} precursors); the rest are zeros", flush=True)

    out = np.concatenate([feat, plm], axis=1)
    np.savez_compressed(os.path.expanduser(a.out), feat=out, offset=offset,
                        n_prot=n_prot, c_prot=c_prot)
    print(f"-> {a.out}  ({out.shape[0]:,} x {out.shape[1]})", flush=True)


if __name__ == "__main__":
    main()
