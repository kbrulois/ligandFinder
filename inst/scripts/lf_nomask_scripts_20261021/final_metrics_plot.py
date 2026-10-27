#!/usr/bin/env python3
"""The final window model's two heads, side by side: the LINEAR head (the CNN
window score off the flattened per-residue class surface) and the xgb head
(per-member xgboost on that same surface). Same trunk, same members -- only the
head differs. Attention is not shown: it is retired.

Reads metrics_final.R's two CSVs (flatten-trunk rows only) and writes
final_model_metrics.png (table) and final_recall_rank.png (small multiples).

Colour = head (categorical slots 1/2 of the reference palette, validated on
white: CVD dE 24.7, normal 33.6, contrast 4.4 / 3.2); the xgb head is also
dashed, so identity never rests on colour alone. Text stays in ink.
The trunk-vs-trunk version is final_metrics_plot_trunks.py.
"""
import os
import csv, math, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = pathlib.Path.home() / "AF2_analysis"
D = A2 / "lf_headcmp"
OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else A2
INK, BODY, MUTED, RULE, GRID = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9", "#E6E5E1"
TRUNK = os.environ.get("TRUNK", "flatten")   # "" default = the curated flatten arm
SUF = "" if TRUNK == "flatten" else f"_{TRUNK}"
HEADS = ["cnn", "xgb_pm"]
COL = {"cnn": "#2a78d6", "xgb_pm": "#eb6834"}
DASH = {"cnn": "-", "xgb_pm": (0, (4, 2.5))}
NAME = {"cnn": "linear head", "xgb_pm": "xgb head"}
SETS = ["cross-validated", "never seen: peptides", "never seen: large products"]

rows = [r for r in csv.DictReader(open(D / "final_metrics.csv")) if r["trunk"] == TRUNK]
curv = [c for c in csv.DictReader(open(D / "final_recall_rank.csv")) if c["trunk"] == TRUNK]
num = lambda r, k: float(r[k]) if r[k] not in ("", "NA") else float("nan")
by = {(r["set"], r["term"], r["head"]): r for r in rows}

# ---------------------------------------------------------------- the table
COLS = [("Ends", "n", "{:.0f}", None), ("ROC AUC", "roc_auc", "{:.4f}", +1),
        ("Norm. PR-AUC", "norm_pr_auc", "{:.3f}", +1), ("AP", "ap", "{:.4f}", +1),
        ("Anchors for\n50% recall", "n50", "{:,.0f}", -1), ("Anchors for\n80% recall", "n80", "{:,.0f}", -1)]
W = 12.6; L, R = 0.55, W - 0.55
ROW, SEP = 0.36, 0.30
nrows = len(SETS) * 2 * len(HEADS)
H = 1.6 + 0.6 + nrows * ROW + len(SETS) * 2 * SEP + 0.95
fig = plt.figure(figsize=(W, H), dpi=180); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10.5, c=BODY, w="normal", ha="left": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
T(L, H - 0.45, "Window model: linear head vs xgb head", 17, INK, "bold")
T(L, H - 0.85, "Same trunk and members (5-fold CV x 10 seeds per terminus), scored over every mature "
  "anchor of the proteome (~3M per terminus).  Other annotated ends are left out of the pool.", 10.5, MUTED)
X = [L, L + 2.4, L + 2.85]                               # set, term, head
XC = [L + 5.2, L + 6.5, L + 7.85, L + 8.95, L + 10.4, L + 11.5]
y = H - 1.5
T(X[0], y, "Set", 9.5, MUTED, "bold"); T(X[1], y, "Term", 9.5, MUTED, "bold")
T(X[2], y, "Head", 9.5, MUTED, "bold")
for (lab, _, _, _), x in zip(COLS, XC):
    T(x, y, lab, 9.5, MUTED, "bold", "right")
y -= 0.28; ax.plot([L, R], [y, y], color=INK, lw=0.8); y -= 0.32
for st in SETS:
    for term in ("C", "N"):
        for hd in HEADS:
            r = by.get((st, term, hd))
            if r is None: continue
            o = by.get((st, term, HEADS[1 - HEADS.index(hd)]))
            if hd == HEADS[0]:
                T(X[0], y, st, 10, INK, "bold"); T(X[1], y, term, 10, INK, "bold")
            # a swatch carries the head's identity; the text stays ink
            ax.plot([X[2], X[2] + 0.30], [y, y], color=COL[hd], lw=2.2, ls=DASH[hd], solid_capstyle="butt")
            T(X[2] + 0.42, y, NAME[hd], 10, BODY)
            for (lab, key, fmt, better), x in zip(COLS, XC):
                v = num(r, key); win = False
                if better is not None and o is not None and not math.isnan(v):
                    ov = num(o, key)
                    win = (v - ov) * better > 0 and fmt.format(v) != fmt.format(ov)
                T(x, y, "-" if math.isnan(v) else fmt.format(v), 10,
                  INK if win else BODY, "bold" if win else "normal", "right")
            y -= ROW
        ax.plot([L, R], [y + ROW / 2 - 0.02, y + ROW / 2 - 0.02], color=RULE, lw=0.6)
        y -= SEP - 0.04
T(L, y + 0.05, "Bold marks the better head on that row (higher AUC / PR; fewer anchors to reach the recall).  "
  "Normalised PR-AUC is AP re-derived at a common 1% prevalence", 9.5, MUTED)
T(L, y - 0.22, "by reweighting the negatives along the same ranking, so the sets compare.  Cross-validated: "
  "each end scored only by the fold that held it out; never seen: all 50 models.", 9.5, MUTED)
p = OUT / f"final_model_metrics{SUF}.png"; fig.savefig(p, facecolor="white"); plt.close(fig)
print("wrote", p)

# ------------------------------------------------- recall vs rank, 3 x 2 panels
fig, axs = plt.subplots(len(SETS), 2, figsize=(11.6, 10.4), dpi=180, sharex=True, sharey=True)
fig.patch.set_facecolor("white")
for i, st in enumerate(SETS):
    for j, term in enumerate(("C", "N")):
        a = axs[i][j]
        for hd in HEADS:
            pts = [(float(c["N"]), float(c["recall"])) for c in curv
                   if c["set"] == st and c["term"] == term and c["head"] == hd]
            if not pts: continue
            xs, ys = zip(*pts)
            a.plot(xs, ys, color=COL[hd], lw=1.9, ls=DASH[hd], label=NAME[hd])
        a.set_xscale("log"); a.set_xlim(1, 3e6); a.set_ylim(0, 1.0)
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(RULE)
        a.tick_params(colors=MUTED, labelsize=9)
        r0 = by.get((st, term, "cnn"))
        n = int(float(r0["n"])) if r0 else 0
        a.set_title(f"{term} terminus  |  {st}  (n = {n})", fontsize=10.5, color=INK, loc="left")
        if j == 0: a.set_ylabel("recall", color=BODY, fontsize=10)
        if i == len(SETS) - 1: a.set_xlabel("anchors examined, best-scored first (log)", color=BODY, fontsize=10)
h, l = axs[0][0].get_legend_handles_labels()
fig.suptitle("Recall vs rank: linear head vs xgb head", x=0.055, y=0.995, ha="left",
             fontsize=15, color=INK, fontweight="bold")
fig.legend(h, l, loc="upper left", ncol=2, frameon=False, fontsize=10,
           bbox_to_anchor=(0.045, 0.968), handlelength=2.6)
fig.tight_layout(rect=(0, 0, 1, 0.935))
p = OUT / f"final_recall_rank{SUF}.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); plt.close(fig)
print("wrote", p)
