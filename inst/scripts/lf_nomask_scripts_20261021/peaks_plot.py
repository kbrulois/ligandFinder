#!/usr/bin/env python3
"""Peaks, not anchors: table + budget curves from peaks_final.R.

Colour = head (validated categorical slots 1/2), line style = threshold
(solid 0.3, dashed 0.5), so identity never rests on colour alone. The budget
curve (ends recovered against peaks examined) is the comparison that holds
across thresholds; ROC / PR over peaks do not compare across T, because each T
emits a different set of peaks.
"""
import csv, math, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = pathlib.Path.home() / "AF2_analysis"; D = A2 / "lf_headcmp"
OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else A2
INK, BODY, MUTED, RULE, GRID = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9", "#E6E5E1"
HEADS = ["cnn", "xgb_pm"]; THR = ["0.3", "0.5"]
COL = {"cnn": "#2a78d6", "xgb_pm": "#eb6834"}
LS = {"0.3": "-", "0.5": (0, (4, 2.5))}
NAME = {"cnn": "linear head", "xgb_pm": "xgb head"}
SETS = ["cross-validated", "never seen: peptides", "never seen: large products"]

rows = list(csv.DictReader(open(D / "peak_metrics.csv")))
ranks = list(csv.DictReader(open(D / "peak_ranks.csv")))
tnorm = lambda t: f"{float(t):.1f}"
num = lambda r, k: float(r[k]) if r[k] not in ("", "NA") else float("nan")
by = {(r["set"], r["term"], r["head"], tnorm(r["thr"])): r for r in rows}

# ---------------------------------------------------------------- the table
COLS = [("Ends", "n_ends", "{:.0f}", None), ("Peaks\n(search space)", "peaks", "{:,.0f}", None),
        ("Ends with a\npeak (ceiling)", "ends_with_peak", "{:.0f}", +1),
        ("Peaks for\n50% of ends", "p50", "{:,.0f}", -1),
        ("ROC AUC\n(over peaks)", "roc_auc", "{:.3f}", None), ("Norm. PR-AUC\n(over peaks)", "norm_pr_auc", "{:.3f}", None)]
W = 12.8; L, R = 0.55, W - 0.55; ROW, SEP = 0.34, 0.28
H = 1.75 + 0.7 + len(SETS) * 2 * 4 * ROW + len(SETS) * 2 * SEP + 1.1
fig = plt.figure(figsize=(W, H), dpi=180); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10.5, c=BODY, w="normal", ha="left": ax.text(x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
T(L, H - 0.45, "Peaks, not anchors: the final window model scored over peaks of its profile", 16, INK, "bold")
T(L, H - 0.85, "Peak = a run of consecutive anchors with score >= T; its centre is the run's maximum.  Positive if the centre is "
  "within 2 residues of a true end; peaks near other real ends are held out.", 10, MUTED)
X = [L, L + 2.4, L + 2.85, L + 4.45]
XC = [L + 5.55, L + 6.85, L + 8.15, L + 9.35, L + 10.6, L + 11.7]
y = H - 1.55
for lab, x in zip(["Set", "Term", "Head", "T"], X): T(x, y, lab, 9.5, MUTED, "bold")
for (lab, _, _, _), x in zip(COLS, XC): T(x, y, lab, 9.5, MUTED, "bold", "right")
y -= 0.33; ax.plot([L, R], [y, y], color=INK, lw=0.8); y -= 0.3
for st in SETS:
    for term in ("C", "N"):
        first = True
        for hd in HEADS:
            for t in THR:
                r = by.get((st, term, hd, t))
                if r is None: continue
                o = by.get((st, term, hd, THR[1 - THR.index(t)]))
                if first: T(X[0], y, st, 10, INK, "bold"); T(X[1], y, term, 10, INK, "bold"); first = False
                ax.plot([X[2], X[2] + 0.30], [y, y], color=COL[hd], lw=2.2, ls=LS[t], solid_capstyle="butt")
                T(X[2] + 0.42, y, NAME[hd], 10, BODY); T(X[3], y, t, 10, BODY)
                for (lab, key, fmt, better), x in zip(COLS, XC):
                    v = num(r, key); win = False
                    if better is not None and o is not None and not math.isnan(v):
                        ov = num(o, key)
                        win = (not math.isnan(ov)) and (v - ov) * better > 0 and fmt.format(v) != fmt.format(ov)
                    T(x, y, "unreached" if math.isnan(v) else fmt.format(v), 10,
                      INK if win else (MUTED if math.isnan(v) else BODY), "bold" if win else "normal", "right")
                y -= ROW
        ax.plot([L, R], [y + ROW / 2 - 0.02, y + ROW / 2 - 0.02], color=RULE, lw=0.6); y -= SEP - 0.04
T(L, y + 0.05, "Bold marks the better threshold for that head (more ends reachable; fewer peaks to reach half of them).  "
  "'Peaks for 50%': each end's rank is 1 + the negative peaks above", 9.2, MUTED)
T(L, y - 0.21, "its best matching peak; 'unreached' when fewer than half the ends have any peak.  ROC and PR over peaks are "
  "shown for reference but do not compare across T (each T emits different peaks).", 9.2, MUTED)
T(L, y - 0.47, "Cross-validated: per fold, each end scored only by the fold that held it out (peak counts are per fold); never seen: "
  "all 50 models.", 9.2, MUTED)
p = OUT / "peak_metrics.png"; fig.savefig(p, facecolor="white"); plt.close(fig); print("wrote", p)

# ------------------------------------------------ budget curves, 3 x 2 panels
NG = sorted(set([int(round(10 ** (i / 30))) for i in range(0, 121)]))
fig, axs = plt.subplots(len(SETS), 2, figsize=(11.6, 10.2), dpi=180, sharex=True, sharey=True)
fig.patch.set_facecolor("white")
for i, st in enumerate(SETS):
    for j, term in enumerate(("C", "N")):
        a = axs[i][j]
        for hd in HEADS:
            for t in THR:
                rk = [float(c["rank"]) for c in ranks if c["set"] == st and c["term"] == term
                      and c["head"] == hd and tnorm(c["thr"]) == t]
                if not rk: continue
                ys = [sum(1 for v in rk if v <= n) / len(rk) for n in NG]
                a.plot(NG, ys, color=COL[hd], ls=LS[t], lw=1.9, label=f"{NAME[hd]}, T = {t}")
        a.set_xscale("log"); a.set_xlim(1, 1e4); a.set_ylim(0, 1.0)
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(RULE)
        a.tick_params(colors=MUTED, labelsize=9)
        r0 = by.get((st, term, "cnn", "0.3")); n = int(float(r0["n_ends"])) if r0 else 0
        a.set_title(f"{term} terminus  |  {st}  (n = {n})", fontsize=10.5, color=INK, loc="left")
        if j == 0: a.set_ylabel("fraction of ends recovered", color=BODY, fontsize=10)
        if i == len(SETS) - 1: a.set_xlabel("negative peaks examined, highest first (log)", color=BODY, fontsize=10)
h, l = axs[0][0].get_legend_handles_labels()
fig.suptitle("Ends recovered against peaks examined", x=0.055, y=0.995, ha="left", fontsize=15, color=INK, fontweight="bold")
fig.legend(h, l, loc="upper left", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.045, 0.968), handlelength=2.6)
fig.tight_layout(rect=(0, 0, 1, 0.935))
p = OUT / "peak_recall_budget.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); plt.close(fig); print("wrote", p)
