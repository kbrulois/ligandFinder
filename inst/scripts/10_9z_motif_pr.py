"""Precision-recall for the window model within each cleavage-motif group.

The scan ranks every anchor in the proteome. Prior knowledge about where a
convertase can cut splits those anchors into groups before the model is
consulted at all. This asks how cleanly the model ranks inside each group:

  all anchors           no prior, the whole scan
  monobasic             K/R at +1 only (or G at +1 with K/R at +2)
  dibasic               K/R at +1 and +2
  G-dibasic (amidation) G at +1, then K/R at +2 and +3

Recall is WITHIN the group -- denominated by the known ends that group contains,
so every curve runs to 100% and "precision at 50% recall" means the same kind of
thing everywhere. That normalisation deliberately hides what each prior costs in
coverage, so the coverage is reported separately (motif_coverage.csv): the
groups are narrow, and between them they do not hold every known end. recall_all
is kept in the curve file for the global-denominator view.

Two arms, as in 10_9z_locus_pr.py. An anchor that is a real peptide end from
OUTSIDE the positive set -- a train end, a ligand_list end, a UniProt
peptide/propeptide end -- is either counted as a false positive or dropped from
the ranking entirely. The second is the fair one: the trunks memorised their
train ends, so those anchors sit near the top of every ranking, and counting
them as mistakes punishes the model hardest exactly where the group is smallest.
Only exact anchors are dropped; a neighbour of a real end stays a false
positive, which is what the near-end negatives trained it to be.

AP and lift are still not comparable between groups -- the prevalence floor
moves with the group -- so the companion plot leads with precision at a matched
within-group recall.

  Rscript inst/scripts/10_9z_motif_codes.R --term C --tables <tables>
  python  inst/scripts/10_9z_motif_pr.py   --term C --tables <tables>
  Rscript inst/scripts/10_9z_motif_pr_plot.R

<tables> also supplies knowns_near.csv (accession, gene, term, anchor, split,
pep_name), ligand_list.csv and uniprot_peptides.csv -- the same tables
10_9z_locus_pr.py reads.
"""
from __future__ import annotations
import argparse, csv, json, os
import numpy as np

RECS = (0.25, 0.50, 0.75, 1.00)
ARMS = ("counted", "excluded")


def groups(motif, L, n):
    return [
        ("all anchors",           "no prior",              np.ones(n, bool)),
        ("monobasic",             "K/R at +1 only",        motif == L["monobasic"]),
        ("dibasic",               "K/R at +1 and +2",      motif == L["dibasic"]),
        ("G-dibasic (amidation)", "G then K/R at +2,+3",   motif == L["G + dibasic"]),
    ]


def load_other_ends(tables, term, pr):
    """Real ends that are not in the positive set: ligand_list, UniProt."""
    out = set()
    for fn, keep in (("ligand_list.csv", None),
                     ("uniprot_peptides.csv", ("peptide", "propeptide"))):
        for r in csv.DictReader(open(os.path.join(tables, fn))):
            if keep is not None and r["type"] not in keep:
                continue
            i = pr.get(r["accession"])
            if i is not None:
                out.add((i, int(r["start"] if term == "N" else r["end"])))
    return out


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

    def find(i, c):
        k = i * 100000 + c
        j = np.searchsorted(ks, k)
        return int(order[j]) if j < len(ks) and ks[j] == k else None

    pos_rows, pos_split, missing = [], [], []
    for r in csv.DictReader(open(os.path.join(a.tables, "knowns_near.csv"))):
        if r["term"] != a.term:
            continue
        i = pr.get(r["accession"])
        j = None if i is None else find(i, int(r["anchor"]))
        if j is None:
            missing.append(r["gene"]); continue
        pos_rows.append(j); pos_split.append(r["split"])
    pos_rows, pos_split = np.array(pos_rows), np.array(pos_split)
    print(f"known {a.term} ends in the scan: {len(pos_rows)} "
          f"(held out {int(np.sum(pos_split == 'val'))}), outside the scan range: {len(missing)}")

    other = np.zeros(len(score), bool)
    for i, c in load_other_ends(a.tables, a.term, pr):
        j = find(i, c)
        if j is not None:
            other[j] = True
    train_rows = pos_rows[pos_split != "val"]
    print(f"unlabelled real {a.term} ends located in the scan: {int(other.sum()):,d}")

    ## ---- coverage: what the groups hold, and what they miss ---------------------
    cov = []
    for tag, split in (("heldout", "val"), ("all", None)):
        sel = pos_rows if split is None else pos_rows[pos_split == split]
        for name in lv:
            cov.append(dict(split=tag, motif=name,
                            n_anchors=int((motif == L[name]).sum()),
                            n_ends=int((motif[sel] == L[name]).sum()),
                            n_ends_all=len(sel)))
    with open(os.path.join(OUT, "motif_coverage.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(cov[0].keys())); w.writeheader(); w.writerows(cov)

    res, curves = {}, []
    for tag, split in (("heldout", "val"), ("all", None)):
        pos = pos_rows if split is None else pos_rows[pos_split == split]
        npos = len(pos)
        lab = np.zeros(len(score), bool); lab[pos] = True
        ## anchors that are real ends but not positives here; never drop a positive
        drop = other.copy()
        if split == "val":
            drop[train_rows] = True
        drop &= ~lab
        tab = []
        for (name, desc, base), arm in ((g, m) for g in groups(motif, L, len(score))
                                        for m in ARMS):
            sel = base & ~drop if arm == "excluded" else base
            s, y = score[sel], lab[sel]
            y = y[np.argsort(-s, kind="stable")]
            ng = int(y.sum())
            if ng == 0:
                continue
            tp = np.cumsum(y); rank = np.arange(1, len(y) + 1)
            prec = tp / rank
            rec = tp / ng            # within the group: always reaches 1
            rec_all = tp / npos      # the global denominator, kept for reference
            row = dict(group=name, desc=desc, split=tag, arm=arm, n_anchors=len(y),
                       n_dropped=int((base & drop).sum()), n_ends=ng, n_ends_all=npos,
                       share_of_ends=ng / npos, prevalence=ng / len(y),
                       ap=float(prec[y].sum() / ng))
            row["lift"] = row["ap"] / row["prevalence"]
            for r in RECS:
                i = int(np.searchsorted(rec, r - 1e-12))
                row[f"n_at_r{int(r*100)}"] = int(rank[i])
                row[f"p_at_r{int(r*100)}"] = float(prec[i])
            tab.append(row)
            idx = np.unique(np.concatenate([
                np.flatnonzero(y),
                np.clip(np.unique(np.logspace(0, np.log10(len(y)), 1500).astype(np.int64)) - 1,
                        0, len(y) - 1)]))
            curves += [dict(split=tag, arm=arm, group=name, rank=int(i + 1),
                            precision=float(prec[i]), recall=float(rec[i]),
                            recall_all=float(rec_all[i])) for i in idx]
        res[tag] = tab

    with open(os.path.join(OUT, "motif_pr.json"), "w") as f:
        json.dump(res, f, indent=1)
    with open(os.path.join(OUT, "motif_curves.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, ["split", "arm", "group", "rank", "precision",
                               "recall", "recall_all"])
        w.writeheader(); w.writerows(curves)
    with open(os.path.join(OUT, "motif_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(res["heldout"][0].keys())); w.writeheader()
        for t in res.values(): w.writerows(t)

    for tag in res:
        for arm in ARMS:
            print(f"\n=== {tag} {a.term} ends, recall within group, "
                  f"other real ends {arm} ===")
            print(f"{'group':24s}{'anchors':>10s}{'dropped':>8s}{'ends':>6s}{'share':>7s}"
                  f"{'prev':>10s}{'AP':>7s}{'lift':>7s}"
                  f"  {'to 50%':>8s}{'prec':>7s}  {'to 100%':>9s}{'prec':>7s}")
            for r in (x for x in res[tag] if x["arm"] == arm):
                print(f"{r['group']:24s}{r['n_anchors']:>10,d}{r['n_dropped']:>8,d}"
                      f"{r['n_ends']:>6d}{r['share_of_ends']:>7.0%}{r['prevalence']:>10.1e}"
                      f"{r['ap']:>7.3f}{r['lift']:>7.0f}"
                      f"  {r['n_at_r50']:>8,d}{r['p_at_r50']:>7.3f}"
                      f"  {r['n_at_r100']:>9,d}{r['p_at_r100']:>7.3f}")


if __name__ == "__main__":
    main()
