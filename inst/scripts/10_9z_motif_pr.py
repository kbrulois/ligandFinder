"""Precision-recall for the window model as prior knowledge narrows the search.

The scan ranks every anchor in the proteome. Prior knowledge about where a
convertase can cut removes most of them before the model is consulted at all.
This asks what that is worth: six search spaces, from no prior to the amidation
motif, each scored by the same 20-seed ensemble.

The one thing that makes the comparison honest is the recall denominator: it is
ALWAYS every known end in the scan, never the ends that survive the restriction.
A prior that throws away two-thirds of the known ends is charged for them, so
the curves share an axis and the trade reads directly off the plot.

AP and lift are reported but are NOT comparable between the spaces -- the
prevalence floor moves with every restriction. The comparable number is the
precision at a matched recall, which is what the companion plot shows.

  Rscript inst/scripts/10_9z_motif_codes.R --term C --tables <tables>
  python inst/scripts/10_9z_motif_pr.py --tables <tables> --term C

<tables> also supplies knowns_near.csv (accession, gene, term, anchor, split,
pep_name), the same table 10_9z_locus_pr.py reads.
"""
from __future__ import annotations
import argparse, csv, json, os
import numpy as np

RECS = (0.05, 0.10, 0.20, 0.30)


def universes(motif, L, n):
    di, gdi, mono, term = L["dibasic"], L["G + dibasic"], L["monobasic"], L["terminus"]
    return [
        ("all anchors",              "no prior",                    np.ones(n, bool)),
        ("basic at +1 (any)",        "mono-, di- or G-dibasic",     np.isin(motif, [di, gdi, mono])),
        ("dibasic or G-dibasic",     "canonical PC site",           np.isin(motif, [di, gdi])),
        ("G-dibasic or protein end", "amidation motif or true end", np.isin(motif, [gdi, term])),
        ("dibasic",                  "K/R at +1,+2",                motif == di),
        ("G-dibasic (amidation)",    "G-K/R-K/R at +1,+2,+3",       motif == gdi),
    ]


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--tables", required=True, help="dir with prec_order.csv, knowns_near.csv")
    ap_.add_argument("--term", default="C", choices=["N", "C"])
    ap_.add_argument("--a2", default=os.path.expanduser("~/AF2_analysis"))
    ap_.add_argument("--score-key", default="score", help="score | score_top5 | score_max")
    a = ap_.parse_args()
    OUT = os.path.join(a.a2, "pepend_motif"); os.makedirs(OUT, exist_ok=True)

    z = np.load(os.path.join(a.a2, f"lf_pepend_scanmax_{a.term}_near.npz"))
    score = z[a.score_key]
    prot, anc = z["prot_idx"].astype(np.int64), z["anchor"].astype(np.int64)
    mz = np.load(os.path.join(a.a2, f"lf_pepend_motif_{a.term}.npz"))
    motif, lv = mz["motif"].astype(np.int64), [str(x) for x in mz["levels"]]
    L = {s: i + 1 for i, s in enumerate(lv)}
    if len(motif) != len(score):
        raise SystemExit(f"motif has {len(motif)} rows, scan has {len(score)} -- rebuild the codes")

    ## locate the known ends in the scan's row order
    pr = {r["accession"]: i for i, r in
          enumerate(csv.DictReader(open(os.path.join(a.tables, "prec_order.csv"))))}
    key = prot * 100000 + anc
    order = np.argsort(key, kind="stable"); ks = key[order]
    pos_rows, pos_split, missing = [], [], []
    for r in csv.DictReader(open(os.path.join(a.tables, "knowns_near.csv"))):
        if r["term"] != a.term:
            continue
        i = pr.get(r["accession"])
        if i is None:
            missing.append(r["gene"]); continue
        j = np.searchsorted(ks, i * 100000 + int(r["anchor"]))
        if j >= len(ks) or ks[j] != i * 100000 + int(r["anchor"]):
            missing.append(r["gene"]); continue
        pos_rows.append(order[j]); pos_split.append(r["split"])
    pos_rows, pos_split = np.array(pos_rows), np.array(pos_split)
    print(f"known {a.term} ends in the scan: {len(pos_rows)} "
          f"(held out {int(np.sum(pos_split == 'val'))}), outside the scan range: {len(missing)}")

    res, curves = {}, []
    for tag, split in (("heldout", "val"), ("all", None)):
        pos = pos_rows if split is None else pos_rows[pos_split == split]
        npos = len(pos)
        lab = np.zeros(len(score), bool); lab[pos] = True
        tab = []
        for name, desc, sel in universes(motif, L, len(score)):
            s, y = score[sel], lab[sel]
            y = y[np.argsort(-s, kind="stable")]
            tp = np.cumsum(y); rank = np.arange(1, len(y) + 1)
            prec, rec = tp / rank, tp / npos
            row = dict(universe=name, desc=desc, split=tag, n_anchors=len(y),
                       n_pos=int(y.sum()), n_pos_all=npos, recall_ceiling=int(y.sum()) / npos,
                       prevalence=int(y.sum()) / len(y),
                       ap=float(prec[y].sum() / npos))
            row["lift"] = row["ap"] / row["prevalence"]
            for r in RECS:
                i = int(np.searchsorted(rec, r - 1e-12))
                hit = i < len(rec) and rec[i] >= r - 1e-12
                row[f"n_at_r{int(r*100)}"] = int(rank[i]) if hit else None
                row[f"p_at_r{int(r*100)}"] = float(prec[i]) if hit else None
            tab.append(row)
            ## thin for plotting: every rank where recall advances, plus a log grid
            idx = np.unique(np.concatenate([
                np.flatnonzero(y),
                np.clip(np.unique(np.logspace(0, np.log10(len(y)), 1500).astype(np.int64)) - 1,
                        0, len(y) - 1)]))
            curves += [dict(split=tag, universe=name, rank=int(i + 1),
                            precision=float(prec[i]), recall=float(rec[i])) for i in idx]
        res[tag] = tab

    with open(os.path.join(OUT, "motif_pr.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "motif_curves.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, ["split", "universe", "rank", "precision", "recall"])
        w.writeheader(); w.writerows(curves)
    with open(os.path.join(OUT, "motif_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(res["heldout"][0].keys())); w.writeheader()
        for t in res.values(): w.writerows(t)

    for tag in res:
        print(f"\n=== {tag} {a.term} ends (denominator {res[tag][0]['n_pos_all']}) ===")
        print(f"{'universe':26s}{'anchors':>10s}{'ends':>6s}{'ceil':>6s}{'prev':>10s}"
              f"{'AP':>7s}{'lift':>7s}  {'to r=.10':>9s}{'prec':>7s}  {'to r=.20':>9s}{'prec':>7s}")
        for r in res[tag]:
            def f(k, fmt):
                v = r[k]
                return "--" if v is None else format(v, fmt)
            print(f"{r['universe']:26s}{r['n_anchors']:>10,d}{r['n_pos']:>6d}"
                  f"{r['recall_ceiling']:>6.2f}{r['prevalence']:>10.1e}{r['ap']:>7.3f}"
                  f"{r['lift']:>7.0f}  {f('n_at_r10', ',d'):>9s}{f('p_at_r10', '.3f'):>7s}"
                  f"  {f('n_at_r20', ',d'):>9s}{f('p_at_r20', '.3f'):>7s}")


if __name__ == "__main__":
    main()
