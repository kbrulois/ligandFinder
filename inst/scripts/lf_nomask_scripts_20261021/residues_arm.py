"""Residue table for ANY arm, assembled from the trained model's own contract.

Kevin, 2026-10-26.

WHY A THIRD BUILDER. The two that existed each hard-code one shape:
`lf_winxgb.plm_residues` appends the parquet to ALL 26 curated channels, and
`residues_esmonly36.py` builds `padding` + the 36 ESM channels. Neither can make
a table for an arm that DROPS some curated channels and adds PLM ones -- for
example `curated_AAswap`, which drops the four AA properties and adds the
21-letter one-hot, so its trunk expects 22 curated + 21 = 43.

This builder does not know any arm's shape. It reads the trained trunk's
`channel_names` and resolves each name, in that exact order, from whichever
source has it:

    a curated name  -> that column of lf_pepend_residues.npz
    a plm_* name    -> that column of the parquet, joined on (accession, index)

so it covers every arm built so far and every one that drops or adds channels,
with no new code. It refuses to write if a name resolves from neither source or
from both.

WHY THE ORDER IS THE WHOLE POINT. `lf_dcnn.scan.build_windows` does
`x[ok] = res.feat[rows[ok]]`, a purely POSITIONAL copy into a
`(k, seq_len, n_channels)` array. A wrong-width table raises -- loudly,
harmlessly. A right-width table in the WRONG ORDER scans silently against the
wrong channels and produces a full proteome of plausible nonsense. Resolving by
name against the model's own contract is what makes that impossible here;
check_resid_order.py then re-verifies the written file against the data.

Note the channel order 10_8 produces: it DROPS first and attaches the PLM
channels afterwards, so the kept curated channels keep their original relative
order and the plm_* ones follow. This script does not assume that -- it simply
follows channel_names -- but that is the order it will see.

    ARM=curated_AAswap ~/.virtualenvs/r-tensorflow/bin/python residues_arm.py
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="residues_arm")
    ap.add_argument("--arm", default=os.environ.get("ARM", "curated_AAswap"))
    ap.add_argument("--plm", default=os.environ.get("PLM", f"{A2}/lf_plm/aa_onehot.parquet"))
    ap.add_argument("--residues-npz", default=f"{A2}/lf_pepend_residues.npz",
                    help="the curated 26-channel table")
    ap.add_argument("--labels-npz", default=f"{A2}/lf_resid_labels.npz",
                    help="its channel_names, and (prot_idx, resno) in row order")
    ap.add_argument("--config", default=None,
                    help="a trained trunk of THIS arm; its channel_names is the contract")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    cfgp = a.config or (f"{A2}/lf_pepend_run_C_cv5_f1_nomask_flatten_"
                        f"{a.arm}_isolated/in/config.json")
    out = a.out or f"{A2}/lf_pepend_residues_{a.arm}.npz"

    want = [str(c) for c in json.load(open(os.path.expanduser(cfgp)))["channel_names"]]
    lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
    base = [str(c) for c in lz["channel_names"]]
    acc = np.asarray([str(v) for v in lz["accession"]])
    prot_idx = lz["prot_idx"].astype("int64")
    resno = lz["resno"].astype("int64")

    with np.load(os.path.expanduser(a.residues_npz), allow_pickle=False) as z:
        feat = np.ascontiguousarray(z["feat"], dtype="float32")
        offset, n_prot, c_prot = z["offset"], z["n_prot"], z["c_prot"]
    n_res = feat.shape[0]
    print(f"model wants {len(want)} channels: {want[:3]} ... {want[-2:]}", flush=True)
    print(f"curated table {n_res:,} x {feat.shape[1]}", flush=True)

    import pandas as pd
    df = pd.read_parquet(os.path.expanduser(a.plm))
    pcols = [c for c in df.columns if c.startswith("plm_")]
    print(f"plm table {len(df):,} rows x {len(pcols)} channels", flush=True)

    # ---- resolve every wanted name, and refuse anything ambiguous ----------
    from_cur = [c for c in want if c in base]
    from_plm = [c for c in want if c in pcols]
    both = sorted(set(from_cur) & set(from_plm))
    missing = [c for c in want if c not in base and c not in pcols]
    if both:
        raise SystemExit(f"channel(s) in BOTH sources, cannot resolve: {both}")
    if missing:
        raise SystemExit(f"channel(s) in neither the curated table nor the parquet: {missing}")
    print(f"resolving {len(from_cur)} curated + {len(from_plm)} plm channel(s)", flush=True)
    dropped = [c for c in base if c not in want]
    if dropped:
        print(f"curated channels this arm DROPS: {dropped}", flush=True)

    # ---- the plm join: one integer key per (protein, residue), as
    # plm_residues does -- a string join over 3M rows is far slower, no safer
    plm = None
    if from_plm:
        code = {s: i for i, s in enumerate(acc)}
        df = df.assign(pi=df["accession"].map(code))
        df = df[df["pi"].notna()]
        key_plm = (df["pi"].to_numpy(dtype="int64") << 20) + df["index"].to_numpy(dtype="int64")
        key_res = (prot_idx << 20) + resno
        order = np.argsort(key_plm, kind="stable")
        kp = key_plm[order]
        vp = df[from_plm].to_numpy(dtype="float32")[order]
        pos = np.clip(np.searchsorted(kp, key_res), 0, kp.size - 1)
        hit = kp[pos] == key_res
        plm = np.zeros((n_res, len(from_plm)), dtype="float32")
        plm[hit] = vp[pos[hit]]
        cov = np.unique(prot_idx[hit]).size
        print(f"matched {hit.sum():,} of {hit.size:,} residues "
              f"({cov} of {acc.size} precursors); the rest are zeros", flush=True)

    # ---- assemble IN THE MODEL'S ORDER ------------------------------------
    cols = []
    for c in want:
        if c in base:
            cols.append(feat[:, base.index(c)])
        else:
            cols.append(plm[:, from_plm.index(c)])
    arr = np.ascontiguousarray(np.column_stack(cols), dtype="float32")
    if arr.shape != (n_res, len(want)):
        raise SystemExit(f"assembled {arr.shape}, expected {(n_res, len(want))}")

    np.savez_compressed(os.path.expanduser(out), feat=arr, offset=offset,
                        n_prot=n_prot, c_prot=c_prot)
    print(f"-> {out}  ({arr.shape[0]:,} x {arr.shape[1]})", flush=True)
    print("column order follows the trunk's channel_names by construction; "
          "run check_resid_order.py to verify it against the data", flush=True)


if __name__ == "__main__":
    main()
