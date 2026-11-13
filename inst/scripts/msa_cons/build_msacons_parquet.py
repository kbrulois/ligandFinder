"""The `curated_AAswap_msacons` arm's per-residue table (Kevin, 2026-11-10).

The arm is curated_AAswap with ONE change: the Aminode conservation channel
`cons_rs_n` is dropped and the MSA conservation score from msa_conservation.py
takes its place. So this table is aa_onehot.parquet (the 21 one-hot channels
AAswap already uses) plus one column, `plm_cons_msa`.

HOW `plm_cons_msa` IS MADE TO MEAN WHAT `cons_rs_n` MEANT:

  * direction and range: Capra-Singh JSD, higher = more conserved, min-max
    scaled to [0, 1] over the residues the curated table covers -- the same
    global (never per-protein) scaling 9.2 applies to cons_rs_n.
  * the residues EXCLUDED BY DESIGN (intracellular / TM / lumenal, and the
    long-protein tail past AF fragment 1) are 0, exactly as cons_rs_n is 0
    there. Excluded = cons_rs_n, min_afm, mean_afm and relASA all exactly 0 in
    lf_pepend_residues.npz. Filling those residues was measured to hurt
    (the abandoned all-residue rebuild); a conservation value there would be a
    partial version of that same change, smuggled in through one channel.
  * a precursor with no MSA would get 0 everywhere, as a residue with no
    Aminode value did -- but that is refused unless --allow-missing, so a
    partly-downloaded MSA set cannot be baked into a trained arm by accident.

THE CHECKS, which refuse to write:
  1. aa_onehot.parquet's rows are the residue table's rows, in order
     (accession, index equal row for row).
  2. the MSA table's residue letter equals the one-hot's at every joined row
     (the MSA query was already checked against the cache sequence; this
     re-checks after the join).

    ~/.virtualenvs/lf-plm/bin/python build_msacons_parquet.py
"""
from __future__ import annotations

import argparse
import os

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")
EXCL_CHANNELS = ["cons_rs_n", "min_afm", "mean_afm", "relASA"]


def main(argv=None):
    import pyarrow as pa
    import pyarrow.parquet as pq

    ap = argparse.ArgumentParser(prog="build_msacons_parquet")
    ap.add_argument("--onehot", default=f"{A2}/lf_plm/aa_onehot.parquet")
    ap.add_argument("--msa-cons", default=f"{A2}/lf_msa/msa_cons_raw.parquet")
    ap.add_argument("--residues-npz", default=f"{A2}/lf_pepend_residues.npz")
    ap.add_argument("--labels-npz", default=f"{A2}/lf_resid_labels.npz")
    ap.add_argument("--out", default=f"{A2}/lf_plm/aa_onehot_msacons.parquet")
    ap.add_argument("--scale-from", default=None,
                    help="reuse the raw-JSD min/max of an earlier build (its .scale.json). "
                         "REQUIRED when rebuilding the channel for arms already trained on it: "
                         "recomputing min-max over a larger MSA set would put the scan's channel "
                         "on a different scale from the one the trunks learned.")
    ap.add_argument("--allow-missing", action="store_true",
                    help="write even if non-excluded residues of some precursors have no MSA")
    a = ap.parse_args(argv)

    oh = pq.read_table(a.onehot)
    n = oh.num_rows
    oh_acc = oh.column("accession").to_numpy(zero_copy_only=False).astype(str)
    oh_idx = oh.column("index").to_numpy().astype("int64")
    oh_aa = oh.column("AA").to_numpy(zero_copy_only=False).astype(str)

    # ---- check 1: the one-hot is in residue-table order ---------------------
    lz = np.load(a.labels_npz, allow_pickle=True)
    acc_tab = np.asarray([str(v) for v in lz["accession"]])
    t_acc = acc_tab[lz["prot_idx"].astype("int64")]
    t_idx = lz["resno"].astype("int64")
    if t_acc.size != n or not (np.array_equal(t_acc, oh_acc) and np.array_equal(t_idx, oh_idx)):
        raise SystemExit("aa_onehot.parquet is not in lf_resid_labels.npz row order; refusing")
    chan = [str(c) for c in lz["channel_names"]]
    with np.load(a.residues_npz, allow_pickle=False) as z:
        feat = z["feat"][:, [chan.index(c) for c in EXCL_CHANNELS]]
    excluded = np.all(feat == 0, axis=1)
    print(f"{n:,} residues; excluded by design: {excluded.sum():,} ({excluded.mean():.1%})", flush=True)

    # ---- join the MSA scores onto the residue rows --------------------------
    mc = pq.read_table(a.msa_cons, columns=["accession", "index", "AA", "jsd"])
    m_acc = mc.column("accession").to_numpy(zero_copy_only=False).astype(str)
    m_idx = mc.column("index").to_numpy().astype("int64")
    m_aa = mc.column("AA").to_numpy(zero_copy_only=False).astype(str)
    m_jsd = mc.column("jsd").to_numpy().astype("float64")
    # key = accession + index; both tables are small enough for a dict of row numbers
    pos = {(ac, i): r for r, (ac, i) in enumerate(zip(m_acc.tolist(), m_idx.tolist()))}
    rows = np.fromiter((pos.get(k, -1) for k in zip(oh_acc.tolist(), oh_idx.tolist())),
                       dtype="int64", count=n)
    have = rows >= 0

    # ---- check 2: letters agree after the join ------------------------------
    bad = have & (m_aa[np.where(have, rows, 0)] != oh_aa)
    if bad.any():
        r = int(np.argmax(bad))
        raise SystemExit(f"{bad.sum():,} joined rows disagree on the residue letter, e.g. "
                         f"{oh_acc[r]} {oh_idx[r]}: one-hot '{oh_aa[r]}' vs MSA "
                         f"'{m_aa[rows[r]]}'; refusing")

    raw = np.where(have, m_jsd[np.where(have, rows, 0)], np.nan)
    covered = have & ~excluded
    if a.scale_from:
        import json
        with open(a.scale_from) as f:
            sc = json.load(f)
        lo, hi = sc["lo"], sc["hi"]
        print(f"scaling PINNED to {a.scale_from}: raw JSD [{lo:.4f}, {hi:.4f}]", flush=True)
        n_out = int(((raw[covered] < lo) | (raw[covered] > hi)).sum())
        if n_out:
            print(f"  {n_out:,} covered residues fall outside that range and are clipped to [0, 1]", flush=True)
    else:
        lo, hi = float(np.nanmin(raw[covered])), float(np.nanmax(raw[covered]))
    val = np.zeros(n, dtype="float32")
    val[covered] = np.clip((raw[covered] - lo) / (hi - lo), 0.0, 1.0).astype("float32")
    no_msa_acc = sorted(set(oh_acc[~have].tolist()))
    print(f"MSA score on {have.sum():,} rows; non-excluded rows with no MSA: "
          f"{(~have & ~excluded).sum():,} in {len(no_msa_acc)} precursors", flush=True)
    if no_msa_acc and (~have & ~excluded).any() and not a.allow_missing:
        raise SystemExit(f"{len(no_msa_acc)} precursors have no MSA score, e.g. {no_msa_acc[:10]}; "
                         "their residues would be written as 0. Finish the MSAs or pass --allow-missing")
    print(f"scaling over covered residues: raw JSD [{lo:.4f}, {hi:.4f}] -> [0, 1]", flush=True)
    print(f"plm_cons_msa on covered: mean {val[covered].mean():.3f} "
          f"median {np.median(val[covered]):.3f}; excluded rows set to 0", flush=True)

    out = oh.append_column("plm_cons_msa", pa.array(val))
    pq.write_table(out, a.out, compression="zstd")
    ## the scaling travels with the parquet, so a later rebuild can pin to it
    import json
    with open(a.out + ".scale.json", "w") as f:
        json.dump({"lo": float(lo), "hi": float(hi), "source": os.path.basename(a.msa_cons),
                   "n_covered": int(covered.sum()), "pinned_from": a.scale_from}, f, indent=2)
    print(f"-> {a.out} ({out.num_rows:,} x {out.num_columns - 3} channels)", flush=True)
    if no_msa_acc:
        print("   no MSA:", ", ".join(no_msa_acc[:30]), "..." if len(no_msa_acc) > 30 else "")


if __name__ == "__main__":
    main()
