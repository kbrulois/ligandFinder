"""ROC over the full scan for the page series, both heads."""
import numpy as np, json, sys, csv
sys.path.insert(0, sys.argv[1])
from scan_auc import auc_from_ranks, boot_pos
from sklearn.metrics import roc_curve
S2 = sys.argv[1]
rows, curves = [], []
for e in json.load(open("/Users/kbrulois/AF2_analysis/pepend_auc_inputs/scan_inputs.json")):
    z = np.load(e["npz"]); pos = z["pos"].astype("float64"); neg = np.sort(z["neg"].astype("float64"))
    roc, _ = auc_from_ranks(pos, neg)
    (lo, hi), _ = boot_pos(pos, neg)
    rows.append(dict(series=e["label"], head=e["head"], n=int(pos.size + neg.size),
                     positives=int(pos.size), prevalence=pos.size / (pos.size + neg.size),
                     roc_auc=round(roc, 4), roc_lo=round(lo, 4), roc_hi=round(hi, 4)))
    y = np.r_[np.ones(pos.size), np.zeros(neg.size)]; p = np.r_[pos, neg]
    fpr, tpr, _ = roc_curve(y, p); k = max(1, fpr.size // 600)
    for a_, b_ in zip(fpr[::k], tpr[::k]):
        curves.append(dict(series=e["label"], head=e["head"], x=round(float(a_), 6), y=round(float(b_), 6)))
    print(f"  {e['label']:24s} {e['head']:20s} ROC {roc:.4f} [{lo:.4f}-{hi:.4f}]", flush=True)
for nm, d in (("page_roc_scan.csv", rows), ("page_roc_scan_curves.csv", curves)):
    with open(f"{S2}/{nm}", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(d[0])); w.writeheader(); w.writerows(d)
