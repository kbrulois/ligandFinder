"""Positive predictive value of each cleavage motif, with and without the model.

PPV is precision under a decision rule rather than a ranking: of the anchors a
rule flags, what fraction are real peptide C ends. For a motif alone the rule is
"flag every anchor in this class", so its PPV is the class's prevalence. Adding
the model turns the rule into "flag the top k of this class", and the PPV is the
precision at rank k.

Two ground truths, because the labelled set is small and very incomplete:

  knowns     the 86 curated C ends the model was built around
  annotated  those plus every ligand_list and UniProt peptide/propeptide C end

PPV against `knowns` is a severe underestimate -- most anchors scored as false
positives have simply never been annotated either way. `annotated` is closer to
what you would find if you looked, and is still a lower bound.

The motif-alone PPV involves no model, so nothing is contaminated. The top-k
columns are a different matter: scored against `annotated` the model is being
credited for ends it trained on. `heldout` is the honest arm -- positives are
the 29 held-out ends, and train and unlabelled ends are dropped from the
ranking rather than counted against it, as in 10_9z_locus_pr.py.

  python inst/scripts/10_9z_motif_ppv.py --tables <tables> --term C
"""
from __future__ import annotations
import argparse, csv, json, math, os
import numpy as np

TOPK = (10, 50, 100, 500)


def wilson(k, n, z=1.96):
    """Wilson score interval -- the normal approximation is useless at k ~ 10."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--tables", required=True)
    ap_.add_argument("--term", default="C", choices=["N", "C"])
    ap_.add_argument("--a2", default=os.path.expanduser("~/AF2_analysis"))
    ap_.add_argument("--score-key", default="score")
    a = ap_.parse_args()
    OUT = os.path.join(a.a2, "pepend_motif"); os.makedirs(OUT, exist_ok=True)

    z = np.load(os.path.join(a.a2, f"lf_pepend_scanmax_{a.term}_near.npz"))
    score = z[a.score_key]
    prot, anc = z["prot_idx"].astype(np.int64), z["anchor"].astype(np.int64)
    mz = np.load(os.path.join(a.a2, f"lf_pepend_motif_{a.term}.npz"))
    motif, lv = mz["motif"].astype(np.int64), [str(x) for x in mz["levels"]]

    pr = {r["accession"]: i for i, r in
          enumerate(csv.DictReader(open(os.path.join(a.tables, "prec_order.csv"))))}
    key = prot * 100000 + anc
    order = np.argsort(key, kind="stable"); ks = key[order]

    def mark(arr, i, c):
        k = i * 100000 + c
        j = np.searchsorted(ks, k)
        if j < len(ks) and ks[j] == k:
            arr[order[j]] = True
            return True
        return False

    known = np.zeros(len(score), bool)
    heldout = np.zeros(len(score), bool)
    for r in csv.DictReader(open(os.path.join(a.tables, "knowns_near.csv"))):
        if r["term"] != a.term:
            continue
        i = pr.get(r["accession"])
        if i is None:
            continue
        mark(known, i, int(r["anchor"]))
        if r["split"] == "val":
            mark(heldout, i, int(r["anchor"]))
    annotated = known.copy()
    for fn, keep in (("ligand_list.csv", None),
                     ("uniprot_peptides.csv", ("peptide", "propeptide"))):
        for r in csv.DictReader(open(os.path.join(a.tables, fn))):
            if keep is not None and r["type"] not in keep:
                continue
            i = pr.get(r["accession"])
            if i is not None:
                mark(annotated, i, int(r["start"] if a.term == "N" else r["end"]))

    base = {"knowns": known, "annotated": annotated, "heldout": heldout}
    for nm, v in base.items():
        print(f"{nm:10s} positives in the scan: {int(v.sum()):,d}  "
              f"overall PPV {v.mean():.2e}")

    rows = []
    for cls in lv + ["ALL (whole scan)"]:
        sel = np.ones(len(score), bool) if cls.startswith("ALL") else motif == lv.index(cls) + 1
        n = int(sel.sum())
        row = dict(motif=cls, n_anchors=n)
        for nm, v in base.items():
            k = int(v[sel].sum())
            lo, hi = wilson(k, n)
            row[f"n_{nm}"] = k
            row[f"ppv_{nm}"] = k / n
            row[f"lo_{nm}"], row[f"hi_{nm}"] = lo, hi
            row[f"enrich_{nm}"] = (k / n) / v.mean()
        ## and with the model: PPV of "top k of this class"
        o = np.argsort(-score[sel], kind="stable")
        for nm in ("knowns", "annotated"):
            cum = np.cumsum(base[nm][sel][o])
            for kk in TOPK:
                row[f"ppv_top{kk}_{nm}"] = float(cum[kk - 1] / kk) if n >= kk else None
        ## honest arm: held-out positives, train and unlabelled ends removed
        hsel = sel & ~(annotated & ~heldout)
        ho = np.argsort(-score[hsel], kind="stable")
        hcum = np.cumsum(heldout[hsel][ho])
        row["n_anchors_heldout_arm"] = int(hsel.sum())
        for kk in TOPK:
            row[f"ppv_top{kk}_heldout"] = (float(hcum[kk - 1] / kk)
                                           if int(hsel.sum()) >= kk else None)
        rows.append(row)

    with open(os.path.join(OUT, "motif_ppv.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    print("\n  top-k columns: left block = vs every annotated end (model saw many of them),"
          " right block = held-out ends only, train/unlabelled dropped\n")
    print(f"{'motif class':14s}{'anchors':>11s}"
          f"{'knowns':>8s}{'PPV':>10s}{'annot':>8s}{'PPV':>9s}{'95% CI':>18s}{'enrich':>8s}"
          f"   {'top10':>7s}{'top50':>7s}{'t100':>7s}{'t500':>7s}"
          f"  |{'top10':>7s}{'top50':>7s}{'t100':>7s}{'t500':>7s}")
    for r in rows:
        ci = f"{r['lo_annotated']:.4f}-{r['hi_annotated']:.4f}"
        tk = "".join(f"{r[f'ppv_top{k}_annotated']:>7.2f}"
                     if r[f"ppv_top{k}_annotated"] is not None else f"{'--':>7s}"
                     for k in TOPK)
        th = "".join(f"{r[f'ppv_top{k}_heldout']:>7.2f}"
                     if r[f"ppv_top{k}_heldout"] is not None else f"{'--':>7s}"
                     for k in TOPK)
        print(f"{r['motif']:14s}{r['n_anchors']:>11,d}"
              f"{r['n_knowns']:>8d}{r['ppv_knowns']:>10.2e}"
              f"{r['n_annotated']:>8d}{r['ppv_annotated']:>9.4f}{ci:>18s}"
              f"{r['enrich_annotated']:>8.1f}   {tk}  |{th}")


if __name__ == "__main__":
    main()
