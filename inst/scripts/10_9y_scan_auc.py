"""ROC / PR for each trunk over the FULL scan: every mature anchor is a negative.

Held-out positives only; train positives are removed rather than counted as
negatives. With ~3M negatives, sklearn on every bootstrap replicate is far too
slow, so AUCs come from the positives' ranks among the sorted negatives -- the
same quantities, O(P log N) a replicate. `check` validates both against sklearn
on the full data before any of it is used.
"""
import numpy as np, sys, csv, json
from sklearn.metrics import roc_auc_score, average_precision_score

def auc_from_ranks(pos, neg_sorted):
    """(roc, ap) for `pos` scored against the fixed, ascending `neg_sorted`."""
    n_neg = neg_sorted.size
    below = np.searchsorted(neg_sorted, pos, side="left")      # negatives strictly below
    equal = np.searchsorted(neg_sorted, pos, side="right") - below
    roc = float(np.mean((below + 0.5 * equal) / n_neg))
    order = np.argsort(-pos, kind="stable")
    above = n_neg - (below + equal)                            # negatives strictly above
    above = above[order] + 0.5 * equal[order]                  # ties: half above
    k = np.arange(1, pos.size + 1)
    return roc, float(np.mean(k / (k + above)))

def check(y, p):
    pos = p[y == 1]; neg = np.sort(p[y == 0])
    r, a = auc_from_ranks(pos, neg)
    return r, a, roc_auc_score(y, p), average_precision_score(y, p)

def boot_pos(pos, neg_sorted, n=2000, seed=0):
    """Resample the positives only -- with 29-33 of them against millions of
    negatives, they carry essentially all the variance."""
    rng = np.random.default_rng(seed); R, A = [], []
    for _ in range(n):
        s = rng.choice(pos, pos.size, replace=True)
        r, a = auc_from_ranks(s, neg_sorted); R.append(r); A.append(a)
    q = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    return q(R), q(A)

if __name__ == "__main__":
    S2 = sys.argv[1]
    d = json.load(open(f"{S2}/scan_inputs.json"))
    rows, curves = [], []
    for e in d:
        z = np.load(e["npz"])
        pos = np.asarray(z["pos"], dtype="float64")
        neg = np.sort(np.asarray(z["neg"], dtype="float64"))
        roc, ap = auc_from_ranks(pos, neg)
        (rlo, rhi), (alo, ahi) = boot_pos(pos, neg)
        prev = pos.size / (pos.size + neg.size)
        rows.append(dict(model=e["label"], head="CNN attention head", n=pos.size + neg.size,
                         positives=pos.size, prevalence=round(prev, 8),
                         roc_auc=round(roc, 4), roc_lo=round(rlo, 4), roc_hi=round(rhi, 4),
                         pr_auc=round(ap, 4), pr_lo=round(alo, 4), pr_hi=round(ahi, 4),
                         pr_lift=round(ap / prev)))
        ## curve points: thin the negatives, keep every positive
        y = np.r_[np.ones(pos.size), np.zeros(neg.size)]
        p = np.r_[pos, neg]
        from sklearn.metrics import roc_curve, precision_recall_curve
        fpr, tpr, _ = roc_curve(y, p); rec, pre = precision_recall_curve(y, p)[1], precision_recall_curve(y, p)[0]
        k = max(1, fpr.size // 600)
        for a_, b_ in zip(fpr[::k], tpr[::k]):
            curves.append(dict(model=e["label"], curve="ROC", x=round(float(a_), 6), y=round(float(b_), 6)))
        k = max(1, rec.size // 600)
        for a_, b_ in zip(rec[::k], pre[::k]):
            curves.append(dict(model=e["label"], curve="PR", x=round(float(a_), 6), y=round(float(b_), 6)))
        print(f"  {e['label']:34s} pos {pos.size:3d}  neg {neg.size:9,d}  ROC {roc:.4f}  PR {ap:.4f}  lift {ap/prev:,.0f}x", flush=True)
    with open(f"{S2}/scan_auc_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(f"{S2}/scan_auc_curves.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(curves[0])); w.writeheader(); w.writerows(curves)
    print(f"{len(rows)} trunks -> scan_auc_table.csv")
