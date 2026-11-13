"""ROC / PR AUC for every window and window-free model, on its own validation split."""
import numpy as np, json, sys, csv
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

S2 = sys.argv[1]; A = "/Users/kbrulois/AF2_analysis/"
rng = np.random.default_rng(0)

def boot(y, p, n=2000):
    """Stratified bootstrap CI: resample positives and negatives separately, so
    every replicate keeps the prevalence that sets the PR AUC floor."""
    ip, ineg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    r, pr = [], []
    for _ in range(n):
        idx = np.concatenate([rng.choice(ip, ip.size, True), rng.choice(ineg, ineg.size, True)])
        r.append(roc_auc_score(y[idx], p[idx])); pr.append(average_precision_score(y[idx], p[idx]))
    q = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    return q(r), q(pr)

rows, curves = [], []
def add(kind, model, split, y, p):
    y = np.asarray(y).astype(int); p = np.asarray(p, dtype="float64")
    roc, prc = roc_auc_score(y, p), average_precision_score(y, p)
    (rlo, rhi), (plo, phi) = boot(y, p)
    prev = float(y.mean())
    rows.append(dict(kind=kind, model=model, split=split, n=int(y.size), positives=int(y.sum()),
                     prevalence=round(prev, 5), roc_auc=round(float(roc), 4),
                     roc_lo=round(rlo, 4), roc_hi=round(rhi, 4),
                     pr_auc=round(float(prc), 4), pr_lo=round(plo, 4), pr_hi=round(phi, 4),
                     pr_lift=round(float(prc) / prev, 1)))
    fpr, tpr, _ = roc_curve(y, p)
    pre, rec, _ = precision_recall_curve(y, p)
    k = max(1, fpr.size // 400)
    for a_, b_ in zip(fpr[::k], tpr[::k]): curves.append(dict(kind=kind, model=model, curve="ROC", x=round(float(a_),5), y=round(float(b_),5)))
    k = max(1, rec.size // 400)
    for a_, b_ in zip(rec[::k], pre[::k]): curves.append(dict(kind=kind, model=model, curve="PR", x=round(float(a_),5), y=round(float(b_),5)))

## ---- window models: the global head's call on held-out candidate windows -----
WIN = [("C_near", "C: current"), ("C_near_slack8", "C: mask fixed"),
       ("C_ctpad", "C: CT-padded, old mask"), ("C_ctpad_slack8", "C: CT-padded, mask fixed"),
       ("N_near", "N: current"), ("N_ntpad_slack8", "N: NT-padded, mask fixed")]
for tag, label in WIN:
    z = np.load(f"{S2}/auc/{tag}.npz", allow_pickle=True)
    y = z["_val__y"]
    for head, hname in (("cnn", "CNN attention head"), ("xgb_pm", "xgb window head")):
        add("window", f"{label} — {hname}", "val windows", y, z[f"_val__{head}"])

## ---- window-free residue models: one residue's own features -----------------
z = np.load(A + "lf_resid_preds.npz", allow_pickle=True)
v = z["in_val"].astype(bool)
for t, tname in (("pep_pocket", "pep_pocket"), ("pep_other", "pep_other"),
                 ("ct_context", "CT context"), ("nt_context", "NT context")):
    y = z["y_" + t][v]
    for m, mname in (("mlp", "MLP"), ("xgb", "XGB")):
        add("no window", f"{tname} — {mname}", "val residues", y, z[f"{t}__{m}"][v])

with open(f"{S2}/auc_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(f"{S2}/auc_curves.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(curves[0])); w.writeheader(); w.writerows(curves)
print(f"{len(rows)} models -> auc_table.csv; {len(curves)} curve points")
for r in rows:
    print(f"  {r['kind']:9s} {r['model']:46s} ROC {r['roc_auc']:.3f} [{r['roc_lo']:.3f}-{r['roc_hi']:.3f}]  "
          f"PR {r['pr_auc']:.3f} [{r['pr_lo']:.3f}-{r['pr_hi']:.3f}]  prev {r['prevalence']:.4f}  lift {r['pr_lift']}x")
