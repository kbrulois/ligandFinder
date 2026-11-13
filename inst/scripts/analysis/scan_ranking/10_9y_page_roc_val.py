"""ROC on the VALIDATION negatives, for the four trunks the pages draw."""
import numpy as np, sys, csv
from sklearn.metrics import roc_auc_score, roc_curve
S2 = sys.argv[1]; rng = np.random.default_rng(0)
PAGE = [("N_near", "N terminus"), ("C_near", "C terminus"),
        ("N_ntpad_slack8", "N terminus, NT-padded"), ("C_ctpad_slack8", "C terminus, CT-padded")]
def boot(y, p, n=2000):
    ip, ineg = np.flatnonzero(y == 1), np.flatnonzero(y == 0); v = []
    for _ in range(n):
        i = np.concatenate([rng.choice(ip, ip.size, True), rng.choice(ineg, ineg.size, True)])
        v.append(roc_auc_score(y[i], p[i]))
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))
rows, curves = [], []
for tag, label in PAGE:
    z = np.load(f"{S2}/auc/{tag}.npz", allow_pickle=True)
    y = z["_val__y"].astype(int)
    for head, hname in (("cnn", "CNN attention head"), ("xgb_pm", "xgb window head")):
        p = z[f"_val__{head}"].astype("float64")
        lo, hi = boot(y, p)
        rows.append(dict(series=label, head=hname, n=int(y.size), positives=int(y.sum()),
                         prevalence=round(float(y.mean()), 5), roc_auc=round(float(roc_auc_score(y, p)), 4),
                         roc_lo=round(lo, 4), roc_hi=round(hi, 4)))
        fpr, tpr, _ = roc_curve(y, p)
        k = max(1, fpr.size // 400)
        for a_, b_ in zip(fpr[::k], tpr[::k]):
            curves.append(dict(series=label, head=hname, x=round(float(a_), 5), y=round(float(b_), 5)))
for nm, d in (("page_roc_val.csv", rows), ("page_roc_val_curves.csv", curves)):
    with open(f"{S2}/{nm}", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(d[0])); w.writeheader(); w.writerows(d)
for r in rows: print(f"  {r['series']:24s} {r['head']:20s} ROC {r['roc_auc']:.3f} [{r['roc_lo']:.3f}-{r['roc_hi']:.3f}]")
