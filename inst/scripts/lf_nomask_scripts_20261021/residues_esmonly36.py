"""Residue table for scanning the ESM-C-ONLY trunks (Kevin, 2026-10-24).

WHY THIS EXISTS. `lf_winxgb.plm_residues` appends the PLM channels to ALL 26
hand-built ones, giving 62 columns -- right for a trunk trained with `--plm`
alone, wrong for this arm. These trunks were trained with `--drop-channels`
as well, so they expect exactly 37:

    ['padding', 'plm_01', ..., 'plm_36']

and `lf_dcnn.scan.build_windows` does `x[ok] = res.feat[rows[ok]]`, a PURELY
POSITIONAL copy into a `(k, seq_len, cfg.n_channels)` array. A table of any
other width raises, and -- worse -- a table of the right width in the wrong
order would scan silently against the wrong channels. So this script does not
trust its own construction: it reads the trained model's channel_names and
asserts the column order it built matches, name for name.

The join and the zero-fill policy are `plm_residues`' own, unchanged: proteins
with no embedding get zeros, as 10_8 does at training time, and the count is
printed rather than assumed. For an ESM-ONLY model those zero rows carry no
information at all, which is why the number matters here more than it ever has.

    ~/.virtualenvs/r-tensorflow/bin/python residues_esmonly36.py
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="residues_esmonly36")
    ap.add_argument("--residues-npz", default=f"{A2}/lf_pepend_residues.npz")
    ap.add_argument("--labels-npz", default=f"{A2}/lf_resid_labels.npz")
    ap.add_argument("--plm", default=f"{A2}/lf_plm/esm_c_pca36.parquet")
    ap.add_argument("--config", default=f"{A2}/lf_pepend_run_C_cv5_f1_nomask_flatten_esmonly36_isolated/in/config.json",
                    help="a trained trunk of THIS arm; its channel_names is the contract")
    ap.add_argument("--curated-config", default=f"{A2}/lf_pepend_run_C_cv5_f1_nomask_flatten_isolated/in/config.json",
                    help="the curated arm, to locate `padding` in the residue table")
    ap.add_argument("--out", default=f"{A2}/lf_pepend_residues_esmonly36.npz")
    a = ap.parse_args(argv)

    want = json.load(open(os.path.expanduser(a.config)))["channel_names"]
    base = json.load(open(os.path.expanduser(a.curated_config)))["channel_names"]
    pad_col = base.index("padding")
    print(f"model wants {len(want)} channels, first three {want[:3]}", flush=True)
    print(f"residue table is in the curated order; padding at column {pad_col}", flush=True)

    with np.load(os.path.expanduser(a.residues_npz), allow_pickle=False) as z:
        feat = np.ascontiguousarray(z["feat"], dtype="float32")
        offset, n_prot, c_prot = z["offset"], z["n_prot"], z["c_prot"]
    if feat.shape[1] != len(base):
        raise SystemExit(f"residue table has {feat.shape[1]} channels, curated config says {len(base)}")
    lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
    acc = np.asarray([str(v) for v in lz["accession"]])
    prot_idx, resno = lz["prot_idx"].astype("int64"), lz["resno"].astype("int64")
    print(f"residues {feat.shape[0]:,} x {feat.shape[1]} channels", flush=True)

    import pandas as pd
    df = pd.read_parquet(os.path.expanduser(a.plm))
    cols = [c for c in df.columns if c.startswith("plm_")]
    print(f"plm table {len(df):,} rows x {len(cols)} channels", flush=True)

    # one integer key per (protein, residue) on both sides -- plm_residues' own
    # approach; a string join over 3M rows is far slower and no safer
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
    cov = np.unique(prot_idx[hit]).size
    print(f"matched {hit.sum():,} of {hit.size:,} residues "
          f"({cov} of {acc.size} precursors); the rest are zeros", flush=True)
    print(f"ZERO-EMBEDDING PRECURSORS: {acc.size - cov} -- for an ESM-only model "
          f"these train and score on no information at all", flush=True)

    # padding first, then the PLM channels: the trained order, not a guess
    out = np.concatenate([feat[:, [pad_col]], plm], axis=1)
    built = ["padding"] + cols
    if built != list(want):
        raise SystemExit("column order does not match the model:\n"
                         f"  built {built[:4]} ... {built[-2:]}\n"
                         f"  want  {list(want)[:4]} ... {list(want)[-2:]}")
    print(f"column order matches the model's channel_names, all {len(built)}", flush=True)

    np.savez_compressed(os.path.expanduser(a.out), feat=out, offset=offset,
                        n_prot=n_prot, c_prot=c_prot)
    print(f"-> {a.out}  ({out.shape[0]:,} x {out.shape[1]})", flush=True)


if __name__ == "__main__":
    main()
