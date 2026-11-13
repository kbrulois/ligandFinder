"""Three motif indicator channels as a per-residue "PLM" table (2026-10-26).

THE ARM. `motifind` is the usual 26 hand-built channels PLUS three binary
indicators of the residue's own identity at every one of the 36 window
positions: is it glycine, is it basic (K or R), is it proline. 26 + 3 = 29.

WHY. `curated_wAAid` (26 + a 21-letter one-hot, 47 channels) roughly doubled
average precision at the C terminus -- 4 of 4 non-degenerate cells, mean 2.07x
-- and halved the anchors needed for 50% recall, but made anchors-for-80%-recall
WORSE in 5 of 6 C cells. Attribution said the model keys on G at offset +1 and
R at +2/+3: the amidation glycine PAM consumes, then the dibasic convertase
site. The N terminus gained nothing, because its signature is "a basic residue
here" and the curated `AA_charge` channel already expresses that exactly.

So the C-terminal signal the one-hot found is a THREE-LETTER fact, carried by 21
channels. The hypothesis this arm tests is that 21 is far more width than the
signal needs and that the extra width is what costs the tail. G and basic are
the motif itself; P is included because proline is the one other letter whose
behaviour is structural rather than physicochemical (it blocks cleavage and
breaks helices) and the four AA properties place it unremarkably.

ALL-ZERO IS NOT AMBIGUOUS. Unlike the 21-letter one-hot, where an all-zero
block could only mean padding, here the great majority of real residues are
all-zero: they are simply none of G, K, R, P. That is fine -- the curated 26
carry an explicit `padding` channel, which is what distinguishes "outside the
sequence" from "an ordinary residue". The indicators do not have to.

WHY A PARQUET AND NOT A NEW FEATURE BUILDER. 10_8's `--plm` is the project's
existing, tested route for "extra per-residue channels": `lf_plm_attach` slots
them in beside the 26, and `lf_winxgb.plm_residues` appends the same columns to
the residue table the exhaustive scan slices. Writing the indicators in that
shape means neither the training script, the scan, nor the window builder needs
a single edit -- and `10_7_pepend_windows.R`, which builds the SHARED residue
cache, is never touched.

The columns must be named `plm_*` or `lf_plm_read`'s `grep("^plm_")` will not
see them, so they are `plm_is_G`, `plm_is_basic`, `plm_is_P`.

THE CHECK, unchanged from `aa_onehot_build.py`. A per-residue table is only
useful if its rows line up with `lf_pepend_residues.npz`, and a misalignment is
silent -- it scans against the wrong residues and produces plausible nonsense.
So this script does not trust its own join: it rebuilds the four scaled AA
property channels from the letters it assigned and compares them, all 3,056,424
rows, against the ones 9.2/10_7 already put in the residue table. Aligned
letters give max|diff| ~ 0. Anything else and it refuses to write.

    ~/.virtualenvs/lf-plm/bin/python motif_ind_build.py
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

A2 = os.path.expanduser("~/AF2_analysis")

# 10_7_pepend_windows.R's aa_props_raw, which is 9.2_add_contact_data.R's, in
# that script's own row order -- the min/max below are taken over these 21 rows,
# so the order does not matter but the membership does.
AA_PROPS_RAW = {
    "A": (1.8, 0.0, 71.08, 6.00),   "R": (-4.5, 1.0, 156.19, 10.76),
    "N": (-3.5, 0.0, 114.10, 5.41), "D": (-3.5, -1.0, 115.09, 2.77),
    "C": (2.5, 0.0, 103.14, 5.07),  "Q": (-3.5, 0.0, 128.13, 5.65),
    "E": (-3.5, -1.0, 129.12, 3.22), "G": (-0.4, 0.0, 57.05, 5.97),
    "H": (-3.2, 0.1, 137.14, 7.59), "I": (4.5, 0.0, 113.16, 6.02),
    "L": (3.8, 0.0, 113.16, 5.98),  "K": (-3.9, 1.0, 128.17, 9.74),
    "M": (1.9, 0.0, 131.19, 5.74),  "F": (2.8, 0.0, 147.18, 5.48),
    "P": (-1.6, 0.0, 97.12, 6.30),  "S": (-0.8, 0.0, 87.08, 5.68),
    "T": (-0.7, 0.0, 101.10, 5.60), "W": (-0.9, 0.0, 186.21, 5.89),
    "Y": (-1.3, 0.0, 163.18, 5.66), "V": (4.2, 0.0, 99.13, 5.96),
    "U": (2.5, -1.0, 150.04, 5.47),
}
PROP_COLS = ["AA_hydro", "AA_charge", "AA_mw", "AA_pI"]

# the indicators, in the order they become channels. Each is a set of letters;
# the name carries the `plm_` prefix lf_plm_read greps for.
INDICATORS = [
    ("plm_is_G", {"G"}),
    ("plm_is_basic", {"K", "R"}),
    ("plm_is_P", {"P"}),
]


def scaled_props() -> dict[str, np.ndarray]:
    """min-max over the 21 rows, exactly as 10_7 does with dplyr::across."""
    raw = np.array([AA_PROPS_RAW[k] for k in AA_PROPS_RAW], dtype="float64")
    lo, hi = raw.min(axis=0), raw.max(axis=0)
    sc = (raw - lo) / (hi - lo)
    return {k: sc[i] for i, k in enumerate(AA_PROPS_RAW)}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="motif_ind_build")
    ap.add_argument("--seqs", default=f"{A2}/lf_plm/precursor_seqs.tsv",
                    help="accession/gene/seq for all 5,616 precursors, exported "
                         "from lf_pepend_residue_cache.rds$prec")
    ap.add_argument("--residues-npz", default=f"{A2}/lf_pepend_residues.npz",
                    help="the curated 26-channel table, for the alignment check")
    ap.add_argument("--labels-npz", default=f"{A2}/lf_resid_labels.npz",
                    help="its (prot_idx, resno) in residue row order")
    ap.add_argument("--out", default=f"{A2}/lf_plm/motif_ind.parquet")
    ap.add_argument("--tol", type=float, default=1e-6)
    a = ap.parse_args(argv)

    seqs = pd.read_csv(os.path.expanduser(a.seqs), sep="\t", dtype=str)
    print(f"sequences: {len(seqs):,} precursors", flush=True)

    lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
    acc_tab = np.asarray([str(v) for v in lz["accession"]])
    prot_idx = lz["prot_idx"].astype("int64")
    resno = lz["resno"].astype("int64")
    n_res = prot_idx.size
    print(f"residue table: {n_res:,} rows over {acc_tab.size:,} precursors", flush=True)

    # ---- the letter of every residue-table row, in its own row order --------
    # residues are contiguous per protein and 1-based, so one pass over the
    # sequences fills the whole column with no join at all.
    by_acc = dict(zip(seqs["accession"], seqs["seq"]))
    miss = [str(x) for x in acc_tab if x not in by_acc]
    if miss:
        raise SystemExit(f"{len(miss)} precursors have no sequence, e.g. {miss[:5]}")
    letters = np.empty(n_res, dtype="U1")
    for pi, acc in enumerate(acc_tab):
        sel = prot_idx == pi
        r = resno[sel]
        s = np.frombuffer(by_acc[acc].encode("ascii"), dtype="S1").astype("U1")
        if r.max() > s.size:
            raise SystemExit(f"{acc}: residue {r.max()} past the {s.size}-residue sequence")
        letters[sel] = s[r - 1]
    alpha = sorted(set(letters.tolist()))
    print(f"alphabet: {len(alpha)} letters {''.join(alpha)}", flush=True)
    unknown = sorted(set(alpha) - set(AA_PROPS_RAW))
    if unknown:
        raise SystemExit(f"letters with no property row: {unknown}")

    # ---- THE ALIGNMENT CHECK ------------------------------------------------
    # rebuild the four scaled property channels from these letters and compare
    # with what is already in the table. This is what catches a bad join.
    with np.load(os.path.expanduser(a.residues_npz), allow_pickle=False) as z:
        feat = z["feat"]
    chan = [str(c) for c in lz["channel_names"]]
    cols = [chan.index(c) for c in PROP_COLS]
    sc = scaled_props()
    want = np.array([sc[c] for c in letters], dtype="float32")
    dev = np.abs(want - feat[:, cols]).max(axis=0)
    print("alignment check, max |diff| per property channel:", flush=True)
    for c, d in zip(PROP_COLS, dev):
        print(f"    {c:10s} {d:.3g}", flush=True)
    if dev.max() > a.tol:
        bad = int(np.argmax(np.abs(want - feat[:, cols]).max(axis=1)))
        raise SystemExit(
            f"letters do NOT line up with the residue table (max |diff| {dev.max():.3g}).\n"
            f"  first worst row {bad}: {acc_tab[prot_idx[bad]]} residue {resno[bad]} "
            f"-> '{letters[bad]}'\n"
            "  the indicators would scan against the wrong residues; refusing to write.")
    print(f"letters agree with the table's own AA properties on all {n_res:,} rows", flush=True)

    # ---- write ---------------------------------------------------------------
    ind_cols = [n for n, _ in INDICATORS]
    ind = np.zeros((n_res, len(INDICATORS)), dtype="float32")
    for j, (_, members) in enumerate(INDICATORS):
        ind[:, j] = np.isin(letters, sorted(members)).astype("float32")

    out = pd.DataFrame({"accession": acc_tab[prot_idx], "index": resno.astype("int32"),
                        "AA": letters})
    out = pd.concat([out, pd.DataFrame(ind, columns=ind_cols)], axis=1)
    os.makedirs(os.path.dirname(os.path.expanduser(a.out)), exist_ok=True)
    out.to_parquet(os.path.expanduser(a.out), index=False, compression="zstd")
    sz = os.path.getsize(os.path.expanduser(a.out)) / 1e6
    print(f"-> {a.out}  ({len(out):,} x {len(ind_cols)} indicator channels, {sz:.1f} MB)",
          flush=True)
    ## report the marginal rate of each indicator and how often none fires: a
    ## typo in a letter set shows up here as a rate that is obviously wrong,
    ## and the all-zero share is the thing a reader will ask about first.
    for j, (nm, members) in enumerate(INDICATORS):
        print(f"   {nm:14s} {''.join(sorted(members)):5s} "
              f"{int(ind[:, j].sum()):>9,} residues  "
              f"{100 * ind[:, j].mean():5.2f}%", flush=True)
    none_fired = int((ind.sum(1) == 0).sum())
    print(f"   none of the three: {none_fired:,} residues "
          f"({100 * none_fired / n_res:.2f}%) -- expected, and distinguished from "
          "padding by the curated `padding` channel", flush=True)
    print(f"   row sums: min {ind.sum(1).min():.0f} max {ind.sum(1).max():.0f} "
          "(max 1 means the three sets are disjoint)", flush=True)


if __name__ == "__main__":
    main()
