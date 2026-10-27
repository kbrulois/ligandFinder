#!/usr/bin/env python3
"""Why do ends that look alike in the UMAP get different calls? (Kevin, 2026-10-22)

1. Shift test: how far is an end's window from its own neighbours on the same
   precursor (anchor +/-1..3) -- the windows the model must tell it apart from --
   compared with the distance to the nearest OTHER end? In the z-scored 936-d
   space, in PCA-50, and after UMAP (.transform of the shifted windows).
2. Separability: can TP be told from FN from the end window alone? 5-fold
   GroupKFold by gene, ROC AUC of (a) logistic regression on all 936 z-scored
   features, (b) on PCA-50, (c) kNN (k = 15) in the UMAP plane. If (a) >> (c),
   the UMAP is throwing information away; if all are ~0.5, the call is not a
   property of the end window by itself.
Same preprocessing and UMAP settings as end_umap.py.
"""
import csv, os, warnings
import numpy as np
import umap
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")

A2 = os.path.expanduser("~/AF2_analysis"); S = A2 + "/lf_nomask_scripts_20261021"
z = np.load(S + "/end_windows.npz"); SH = list(z["shifts"])
meta = list(csv.DictReader(open(S + "/end_windows_meta.csv")))
rows_out = []
for T in "CN":
    X = z[f"X_{T}"].astype("float64"); Xs = z[f"Xshift_{T}"].astype("float64")
    m = [r for r in meta if r["term"] == T]; n = len(m)
    sd = X.std(0); keep = sd > 1e-9; mu = X[:, keep].mean(0)
    Z = (X[:, keep] - mu) / sd[keep]; Zs = (Xs[:, :, keep] - mu) / sd[keep]
    pca = PCA(n_components=50, random_state=42).fit(Z)
    P = pca.transform(Z); Ps = pca.transform(Zs.reshape(-1, Z.shape[1])).reshape(n, len(SH), -1)
    um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit(P)
    U = um.embedding_; Us = um.transform(Ps.reshape(-1, 50)).reshape(n, len(SH), 2)

    # ---- 1. shift test -------------------------------------------------------
    print(f"\n== {T} terminus: {n} ends ==")
    for nm, A, As in (("936-d z-scored", Z, Zs), ("PCA-50", P, Ps), ("UMAP 2-d", U, Us)):
        D = np.sqrt(((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
        nn_other = D.min(1)
        line = f"  {nm:<15} nearest other end: median {np.median(nn_other):6.2f}   own neighbour at"
        for j, d in enumerate(SH):
            if abs(d) in (1, 3) and d > 0 or d == -1:
                dd = np.sqrt(((A - As[:, j]) ** 2).sum(-1))
                line += f"  {d:+d}: {np.median(dd):6.2f} ({100 * np.mean(dd < nn_other):.0f}% closer)"
        print(line)
        rows_out.append(dict(term=T, space=nm, nn_other=float(np.median(nn_other)),
                             **{f"shift{d:+d}": float(np.median(np.sqrt(((A - As[:, j]) ** 2).sum(-1))))
                                for j, d in enumerate(SH)}))

    # ---- 2. separability -----------------------------------------------------
    y = np.array([r["call"] == "TP" for r in m]).astype(int)
    g = np.array([r["gene"] for r in m])
    internal = Xs.reshape(n, len(SH), 36, -1)  # placeholder to keep shapes obvious
    padded = X.reshape(n, 36, -1)[:, :, -1].max(1) > 0
    for subset, mask in (("all ends", np.ones(n, bool)), ("internal ends only", ~padded)):
        ii = np.flatnonzero(mask); yy, gg = y[ii], g[ii]
        res = {}
        for nm, feats, mk in (("logistic, 936-d", Z, lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=5000))),
                              ("logistic, PCA-50", P, lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.05, max_iter=5000))),
                              ("kNN-15, UMAP 2-d", U, lambda: KNeighborsClassifier(15))):
            pr = np.zeros(len(ii))
            for tr, te in GroupKFold(5).split(feats[ii], yy, gg):
                pr[te] = mk().fit(feats[ii][tr], yy[tr]).predict_proba(feats[ii][te])[:, 1]
            res[nm] = roc_auc_score(yy, pr)
        print(f"  TP vs FN, {subset:<19} (n={len(ii)}, {yy.mean() * 100:.0f}% TP): " +
              "   ".join(f"{k} AUC {v:.2f}" for k, v in res.items()))
        rows_out.append(dict(term=T, space=f"TPvsFN {subset}", n=len(ii), **res))
with open(A2 + "/lf_headcmp/end_window_diag.csv", "w", newline="") as f:
    keys = sorted({k for r in rows_out for k in r}, key=lambda k: (k not in ("term", "space"), k))
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows_out)
