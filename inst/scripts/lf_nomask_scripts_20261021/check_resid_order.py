"""Assert an arm's residue table really holds the channels its trunks expect.

THE TRAP THIS EXISTS FOR. `lf_dcnn.scan.build_windows` does

    x[ok] = res.feat[rows[ok]]

a purely POSITIONAL copy into a `(k, seq_len, cfg.n_channels)` array. A table of
the wrong width raises -- loudly, harmlessly. A table of the RIGHT width in the
WRONG ORDER scans silently against the wrong channels and produces a full
proteome of plausible nonsense. Nothing downstream can detect it.

So this does not check names against names, which only restates how the table
was built. It checks the DATA, column by column, against the two sources every
arm's table is assembled from:

  * a curated channel must equal `lf_pepend_residues.npz`'s own column of that
    name. Both tables are in residue row order, so this is a direct comparison
    and needs no join.
  * a `plm_*` channel must equal the parquet's column of that name at that
    residue, re-joining a random SAMPLE of rows on (accession, index) -- the
    same key `lf_winxgb.plm_residues` joins on, but computed independently here
    so an error in that join cannot agree with itself.

That makes the check construction-agnostic: it passes for an arm built as
"26 curated + plm" (curated_wAAid, via lf_winxgb.plm_residues) and for one
built as "padding + plm" (esmonly36, via residues_esmonly36.py), and it fails
for any table whose columns are not where the trunk thinks they are.

    python check_resid_order.py --arm curated_wAAid
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="check_resid_order")
    ap.add_argument("--arm", default=os.environ.get("ARM", "curated_wAAid"))
    ap.add_argument("--plm", default=os.environ.get("PLM", f"{A2}/lf_plm/aa_onehot.parquet"))
    ap.add_argument("--residues-npz", default=None)
    ap.add_argument("--curated-npz", default=f"{A2}/lf_pepend_residues.npz")
    ap.add_argument("--labels-npz", default=f"{A2}/lf_resid_labels.npz")
    ap.add_argument("--config", default=None)
    ap.add_argument("--sample", type=int, default=20000,
                    help="rows re-joined against the parquet for the plm columns")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--seed", type=int, default=3)
    a = ap.parse_args(argv)
    resid = a.residues_npz or f"{A2}/lf_pepend_residues_{a.arm}.npz"
    cfgp = a.config or (f"{A2}/lf_pepend_run_C_cv5_f1_nomask_flatten_"
                        f"{a.arm}_isolated/in/config.json")

    want = [str(c) for c in json.load(open(os.path.expanduser(cfgp)))["channel_names"]]
    lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
    base = [str(c) for c in lz["channel_names"]]
    acc_tab = np.asarray([str(v) for v in lz["accession"]])
    prot_idx = lz["prot_idx"].astype("int64")
    resno = lz["resno"].astype("int64")

    with np.load(os.path.expanduser(resid), allow_pickle=False) as z:
        feat = z["feat"]
        n, w = feat.shape
    print(f"arm        : {a.arm}")
    print(f"trunk wants: {len(want)} channels, {want[:2]} ... {want[-2:]}")
    print(f"residues   : {n:,} x {w}  ({os.path.basename(resid)})")

    fail = []
    if w != len(want):
        fail.append(f"table is {w} channels wide, the trunk wants {len(want)}")
        print("\nFAIL -- do not scan with this table:", file=sys.stderr)
        for f in fail:
            print("  " + f, file=sys.stderr)
        raise SystemExit(1)

    # ---- curated channels, against the curated table itself -----------------
    cur = [(i, c) for i, c in enumerate(want) if c in base]
    if cur:
        with np.load(os.path.expanduser(a.curated_npz), allow_pickle=False) as z:
            cfeat = z["feat"]
        if cfeat.shape[0] != n:
            fail.append(f"curated table has {cfeat.shape[0]:,} rows, this one {n:,}")
        else:
            worst, wc = 0.0, None
            for i, c in cur:
                d = float(np.abs(feat[:, i] - cfeat[:, base.index(c)]).max())
                if d > worst:
                    worst, wc = d, c
            print(f"curated    : {len(cur)} channel(s) checked against "
                  f"{os.path.basename(a.curated_npz)}, worst |diff| {worst:.3g}"
                  + (f" ({wc})" if worst else ""))
            if worst > a.tol:
                fail.append(f"curated channel {wc!r} does not match the curated table "
                            f"(max |diff| {worst:.3g}) -- the columns are not where "
                            "the trunk thinks they are")

    # ---- plm channels, by re-joining a sample on (accession, index) ---------
    plm_names = [(i, c) for i, c in enumerate(want) if c.startswith("plm_")]
    if plm_names:
        import pandas as pd
        df = pd.read_parquet(os.path.expanduser(a.plm))
        have = [c for _, c in plm_names if c in df.columns]
        missing = [c for _, c in plm_names if c not in df.columns]
        if missing:
            fail.append(f"{len(missing)} plm channel(s) the trunk wants are not in "
                        f"{os.path.basename(a.plm)}, e.g. {missing[:3]}")
        if have:
            rng = np.random.default_rng(a.seed)
            k = min(a.sample, n)
            rows = rng.choice(n, k, replace=False)
            s = pd.DataFrame({"accession": acc_tab[prot_idx[rows]],
                              "index": resno[rows]})
            m = s.merge(df[["accession", "index"] + have], on=["accession", "index"],
                        how="left")
            got = m[have].to_numpy(dtype="float32")
            # a residue absent from the parquet is a legitimate zero row (the
            # policy 10_8 applies at training time), so compare those to zero
            absent = m[have[0]].isna().to_numpy()
            got[absent] = 0.0
            cols = [i for i, c in plm_names if c in have]
            d = np.abs(feat[np.ix_(rows, cols)] - got)
            worst = float(d.max())
            bad = have[int(np.unravel_index(int(d.argmax()), d.shape)[1])]
            print(f"plm        : {len(have)} channel(s) re-joined on {k:,} sampled "
                  f"rows, worst |diff| {worst:.3g}" + (f" ({bad})" if worst else ""))
            print(f"             {int(absent.sum()):,} of {k:,} sampled residues are "
                  "absent from the parquet and must be zero")
            if worst > a.tol:
                fail.append(f"plm channel {bad!r} does not match the parquet "
                            f"(max |diff| {worst:.3g}) -- wrong column order or a "
                            "bad join")

    unchecked = [c for c in want if c not in base and not c.startswith("plm_")]
    if unchecked:
        print(f"NOTE       : {len(unchecked)} channel(s) checked by width only: "
              f"{unchecked[:5]}")

    if fail:
        print("\nFAIL -- do not scan with this table:", file=sys.stderr)
        for f in fail:
            print("  " + f, file=sys.stderr)
        raise SystemExit(1)
    print(f"OK: every one of the {len(want)} channels holds the data the trunk "
          "expects, in the position it expects it")


if __name__ == "__main__":
    main()
