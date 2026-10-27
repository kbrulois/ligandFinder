#!/usr/bin/env python3
"""UMAP of the flattened input windows of every known peptide end
(Kevin, 2026-10-22), faceted terminus x seen / never seen.

Per terminus, all ends together (seen = the 5-fold CV knowns; never seen =
peptides and large cleavage products): the 36 x 26 window flattened to 936
features, constant columns dropped, each column z-scored, PCA to 50, UMAP
(15 neighbours, min_dist 0.3, seed 42). One embedding per terminus, so the
three panels of a row share coordinates. Points are coloured by the linear
head's call at T = 0.2 (end_calls.R).

Also reports, in PCA space, each end's distance to the nearest SEEN end (seen
ends: leave-one-out) -- is a never-seen window further from the training
windows than a seen one is from its neighbours?

    python3 end_umap.py            (system python: umap-learn, sklearn)
"""
import csv, os
import numpy as np
import umap
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
import matplotlib.patheffects as pe
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = os.path.expanduser("~/AF2_analysis"); S = A2 + "/lf_nomask_scripts_20261021"
INK, BODY, MUTED, RULE, BG = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9", "#DCE1E6"
TPC, FNC = "#4a3aa7", "#eda100"          # violet / amber: apart in lightness as well as hue
SETS = ["seen", "never seen: peptides", "never seen: large products"]
TITLE = {"seen": "seen (cross-validated knowns)", "never seen: peptides": "never seen: peptides",
         "never seen: large products": "never seen: large products"}

z = np.load(S + "/end_windows.npz")
meta = list(csv.DictReader(open(S + "/end_windows_meta.csv")))
ch = list(z["channels"]); pad = ch.index("padding"); nch = len(ch)
coords, stats = {}, []
for T in "CN":
    X = z[f"X_{T}"].astype("float64"); m = [r for r in meta if r["term"] == T]
    assert len(m) == X.shape[0]
    sd = X.std(0); keep = sd > 1e-9
    Z = (X[:, keep] - X[:, keep].mean(0)) / sd[keep]
    P = PCA(n_components=50, random_state=42).fit_transform(Z)
    U = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=42).fit_transform(P)
    ## name the islands by what sets them apart: where the window runs off the
    ## precursor (padding), or channels that are all zero (intracellular / TM)
    W3 = X.reshape(len(m), 36, nch); anc = 28 if T == "C" else 8
    lab = DBSCAN(eps=1.2, min_samples=5).fit_predict(U); names = []
    for c in sorted(set(lab) - {-1}):
        ii = np.flatnonzero(lab == c)
        if len(ii) < 20: continue
        pd_ = W3[ii, :, pad]; real = 1 - pd_
        ## the DSSP one-hots sum to 1 on every residue that has structure channels
        ss = [k for k, c in enumerate(ch) if c.startswith("SS_")]
        zero = (W3[ii][:, :, ss].sum(-1) * real).sum() / real.sum() < 0.5
        before = (pd_[:, :anc - 1].max(1) > 0).mean(); after = (pd_[:, anc:].max(1) > 0).mean()
        txt = ("all-zero channels\n(intracellular / TM)" if zero else
               "end at the precursor's\nC-terminus" if after > 0.8 else
               "end right after\nthe signal peptide" if before > 0.8 else "internal end")
        names.append((U[ii, 0].mean(), U[ii, 1].max(), txt))
    coords[T] = (U, m, names)
    ## distance to the nearest seen end, in PCA space
    seen = np.array([r["set"] == "seen" for r in m])
    D = np.sqrt(((P[:, None, :] - P[None, seen, :]) ** 2).sum(-1))
    si = np.flatnonzero(seen)
    for j, s_idx in enumerate(si): D[s_idx, j] = np.inf          # leave-one-out
    nn = D.min(1)
    padded = X.reshape(len(m), 36, nch)[:, :, pad].max(1) > 0
    for i, r in enumerate(m):
        r.update(umap1=f"{U[i, 0]:.4f}", umap2=f"{U[i, 1]:.4f}", nn_seen=f"{nn[i]:.3f}", padded=int(padded[i]))
    for st in SETS:
        ii = [i for i, r in enumerate(m) if r["set"] == st]
        for c in ("TP", "FN", "all"):
            jj = [i for i in ii if c == "all" or m[i]["call"] == c]
            if jj:
                stats.append(dict(term=T, set=st, call=c, n=len(jj), nn_seen_median=float(np.median(nn[jj])),
                                  padded_pct=100 * float(np.mean(padded[jj]))))
    print(f"{T}: {X.shape[0]} windows, {keep.sum()} of {X.shape[1]} features vary; "
          f"PCA 50 explains {100 * PCA(50).fit(Z).explained_variance_ratio_.sum():.0f}%")

with open(S + "/end_umap_coords.csv", "w", newline="") as f:
    rows = coords["C"][1] + coords["N"][1]
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(A2 + "/lf_headcmp/end_umap_nn_seen.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(stats[0])); w.writeheader(); w.writerows(stats)
for s in stats:
    print(f"  {s['term']} {s['set']:<27} {s['call']:<3} n={s['n']:<4} median dist to nearest seen {s['nn_seen_median']:.1f}"
          f"   padded {s['padded_pct']:.0f}%")
st_by = {(s["term"], s["set"], s["call"]): s for s in stats}

fig, axs = plt.subplots(2, 3, figsize=(13.6, 9.4), dpi=180)
fig.patch.set_facecolor("white")
for i, T in enumerate("CN"):
    U, m, names = coords[T]
    for j, st in enumerate(SETS):
        a = axs[i][j]
        a.scatter(U[:, 0], U[:, 1], s=9, c=BG, lw=0, zorder=1)
        for c, col, mk, sz in (("FN", FNC, "D", 20), ("TP", TPC, "o", 22)):
            ii = [k for k, r in enumerate(m) if r["set"] == st and r["call"] == c]
            a.scatter(U[ii, 0], U[ii, 1], s=sz, c=col, marker=mk, edgecolors="white", linewidths=0.5, zorder=3)
        n_tp = st_by.get((T, st, "TP"), {}).get("n", 0); n_fn = st_by.get((T, st, "FN"), {}).get("n", 0)
        a.set_title(f"{T} terminus  |  {TITLE[st]}", fontsize=10.5, color=INK, loc="left", pad=16)
        dist = st_by[(T, st, "all")]["nn_seen_median"]
        a.text(0, 1.015, f"{n_tp} TP, {n_fn} FN at T = 0.2   ·   median distance to nearest seen end {dist:.1f}",
               transform=a.transAxes, fontsize=8, color=MUTED, va="bottom")
        for x0, y0, txt in names:
            a.annotate(txt, (x0, y0), xytext=(0, 6), textcoords="offset points", ha="center", va="bottom",
                       fontsize=7.6, color=BODY, path_effects=[pe.withStroke(linewidth=2.5, foreground="white")], zorder=5)
        a.margins(0.08, 0.14)
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values(): s.set_color(RULE)
        if j == 0: a.set_ylabel("UMAP 2", color=BODY, fontsize=9.5)
        if i == 1: a.set_xlabel("UMAP 1", color=BODY, fontsize=9.5)
fig.suptitle("UMAP of the flattened input windows of known peptide ends", x=0.04, y=0.995, ha="left",
             fontsize=15, color=INK, fontweight="bold")
fig.text(0.04, 0.952, "36 residues x 26 input channels per end window, flattened, z-scored, PCA 50 -> UMAP; "
         "one embedding per terminus, so a row's panels share coordinates.  Grey: every end of that terminus.",
         fontsize=9.2, color=MUTED, ha="left")
fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=TPC, mec="white", label="TP: peak on the end at T = 0.2"),
                    Line2D([], [], marker="D", ls="", ms=5.5, mfc=FNC, mec="white", label="FN: no peak centred on the end"),
                    Line2D([], [], marker="o", ls="", ms=5, mfc=BG, mec=BG, label="other ends, same terminus")],
           loc="upper left", bbox_to_anchor=(0.035, 0.94), ncol=3, frameon=False, fontsize=9.3)
fig.tight_layout(rect=(0, 0.035, 1, 0.905))
fig.text(0.04, 0.006, "Calls: CNN linear head; seen ends on the profile of the fold that held each one out, never-seen ends on the "
         "all-50 profile.  Distance: Euclidean in the 50-PC space, to the nearest seen end (seen ends: leave-one-out).",
         fontsize=8.6, color=MUTED, ha="left", va="bottom")
p = A2 + "/end_window_umap.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); print("wrote", p)
