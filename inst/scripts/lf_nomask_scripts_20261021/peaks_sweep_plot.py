#!/usr/bin/env python3
"""Threshold sweep for peak calling (peaks_final.R with THRESH=0.1..0.5):
an expanded table (one row per set x terminus x T, the two heads side by side)
and a chart of the two numbers that pick T -- how many ends have a peak, and
how many negative peaks it takes to reach half of them.

ROC / PR over peaks are left out on purpose: each T emits a different set of
peaks, so they do not compare across T (they are in peak_metrics.csv).
Colour = head (validated slots 1/2); line style = terminus (C solid, N dashed).
"""
import csv, math, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = pathlib.Path.home() / "AF2_analysis"; D = A2 / "lf_headcmp"
OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else A2
INK, BODY, MUTED, RULE, GRID = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9", "#E6E5E1"
HEADS = ["cnn", "xgb_pm"]
COL = {"cnn": "#2a78d6", "xgb_pm": "#eb6834"}
NAME = {"cnn": "linear head", "xgb_pm": "xgb head"}
SETS = ["cross-validated", "never seen: peptides", "never seen: large products"]
TLS = {"C": "-", "N": (0, (4, 2.5))}

rows = list(csv.DictReader(open(D / "peak_metrics.csv")))
num = lambda r, k: float(r[k]) if r[k] not in ("", "NA") else float("nan")
THR = sorted({round(float(r["thr"]), 2) for r in rows})
by = {(r["set"], r["term"], r["head"], round(float(r["thr"]), 2)): r for r in rows}

# ---------------------------------------------------------------- the table
MET = [("Peaks", "peaks", lambda v, r: f"{v:,.0f}", None),
       ("Ends with\na peak", "ends_with_peak", lambda v, r: f"{v:.0f} ({100 * v / num(r, 'n_ends'):.0f}%)", +1),
       ("Peaks for\n50% of ends", "p50", lambda v, r: f"{v:,.0f}", -1),
       ("Peaks for\n80% of ends", "p80", lambda v, r: f"{v:,.0f}", -1),
       ## over each threshold's OWN peaks: they compare the heads at one T, not
       ## one T against another, so they are never bolded
       ("ROC AUC\n(over peaks)", "roc_auc", lambda v, r: f"{v:.3f}", None),
       ("Norm. PR-AUC\n(over peaks)", "norm_pr_auc", lambda v, r: f"{v:.3f}", None)]
W = 18.4; L, R = 0.55, W - 0.55; ROW, SEP = 0.27, 0.30
nrow = len(SETS) * 2 * len(THR)
H = 2.2 + nrow * ROW + len(SETS) * 2 * SEP + 1.15
fig = plt.figure(figsize=(W, H), dpi=170); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10, c=BODY, w="normal", ha="left": ax.text(x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
T(L, H - 0.42, "Peak calling: threshold sweep", 17, INK, "bold")
T(L, H - 0.78, "Peak = a run of consecutive anchors with score >= T; centre = the run's maximum; positive if within 2 residues "
  "of a true end, peaks near other real ends held out.", 9.8, MUTED)
X0 = [L, L + 2.45, L + 2.95]                         # set, term, T
BLK = {"cnn": L + 3.6, "xgb_pm": L + 10.6}           # left edge of each head's block
OFF = [0.85, 2.15, 3.3, 4.45, 5.5, 6.7]                         # right edges inside a block
yh = H - 1.32
for hd in HEADS:
    x0 = BLK[hd]
    ax.plot([x0, x0 + 0.3], [yh + 0.05, yh + 0.05], color=COL[hd], lw=2.4, solid_capstyle="butt")
    T(x0 + 0.42, yh + 0.05, NAME[hd], 10.5, INK, "bold")
    ax.plot([x0, x0 + OFF[-1]], [yh - 0.17, yh - 0.17], color=RULE, lw=0.8)
y = yh - 0.48
T(X0[0], y, "Set", 9.3, MUTED, "bold"); T(X0[1], y, "Term", 9.3, MUTED, "bold"); T(X0[2], y, "T", 9.3, MUTED, "bold")
for hd in HEADS:
    for (lab, _, _, _), o in zip(MET, OFF): T(BLK[hd] + o, y, lab, 9.3, MUTED, "bold", "right")
y -= 0.33; ax.plot([L, R], [y, y], color=INK, lw=0.8); y -= 0.27
for st in SETS:
    for term in ("C", "N"):
        best = {}
        for hd in HEADS:
            for _, key, _, better in MET:
                if better is None: continue
                vals = [num(by[(st, term, hd, t)], key) for t in THR if (st, term, hd, t) in by]
                vals = [v for v in vals if not math.isnan(v)]
                best[(hd, key)] = (max(vals) if better > 0 else min(vals)) if vals else None
        for k, t in enumerate(THR):
            if k == 0: T(X0[0], y, st, 9.8, INK, "bold"); T(X0[1], y, term, 9.8, INK, "bold")
            T(X0[2], y, f"{t:.2f}", 9.8, BODY)
            for hd in HEADS:
                r = by.get((st, term, hd, t))
                if r is None: continue
                for (lab, key, fmt, better), o in zip(MET, OFF):
                    v = num(r, key)
                    win = better is not None and not math.isnan(v) and best.get((hd, key)) == v
                    T(BLK[hd] + o, y, "unreached" if math.isnan(v) else fmt(v, r), 9.8,
                      INK if win else (MUTED if math.isnan(v) else BODY), "bold" if win else "normal", "right")
            y -= ROW
        ax.plot([L, R], [y + ROW / 2 - 0.02, y + ROW / 2 - 0.02], color=RULE, lw=0.6); y -= SEP - 0.03
T(L, y + 0.04, "Bold marks each head's best threshold within a set and terminus (most ends with a peak; fewest peaks to "
  "reach 50% / 80% of ends).  'Peaks for x%': each end's rank is 1 + the negative", 9, MUTED)
T(L, y - 0.21, "peaks above its best matching peak; 'unreached' when fewer than x% of ends have any peak.  Cross-validated: "
  "per fold, each end scored only by the fold that held it out (peak counts per fold); never seen: all 50 models.", 9, MUTED)
T(L, y - 0.46, "ROC AUC and normalised PR-AUC (AP at a common 1% prevalence) are computed over each threshold's own peaks, ranked by "
  "height: they compare the two heads at the same T, but not one T against another, since each T emits a different set of peaks.", 9, MUTED)
p = OUT / "peak_metrics_sweep.png"; fig.savefig(p, facecolor="white"); plt.close(fig); print("wrote", p)

# --------------------------------------- the two numbers that pick T, by set
fig, axs = plt.subplots(4, len(SETS), figsize=(13.2, 13.2), dpi=180, sharex=True)
fig.patch.set_facecolor("white")
for j, st in enumerate(SETS):
    for i, (key, ylab) in enumerate([("ends_with_peak", "ends with a peak (%)"),
                                     ("p50", "negative peaks to reach\n50% of ends (log)"),
                                     ("roc_auc", "ROC AUC (over peaks)"),
                                     ("norm_pr_auc", "norm. PR-AUC (over peaks)")]):
        a = axs[i][j]
        for hd in HEADS:
            for term in ("C", "N"):
                xs, ys = [], []
                for t in THR:
                    r = by.get((st, term, hd, t))
                    if r is None: continue
                    v = num(r, key)
                    if key == "ends_with_peak": v = 100 * v / num(r, "n_ends")
                    if not math.isnan(v): xs.append(t); ys.append(v)
                if xs:
                    a.plot(xs, ys, color=COL[hd], ls=TLS[term], lw=1.9, marker="o", ms=4.5,
                           mfc="white", mew=1.2, label=f"{NAME[hd]}, {term} terminus")
        if key == "p50":
            ## one shared log range across the three sets, with plain labels, so
            ## the cost reads the same in every panel
            a.set_yscale("log"); a.set_ylim(200, 20000)
            a.yaxis.set_major_locator(FixedLocator([300, 1000, 3000, 10000]))
            a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
            a.yaxis.set_minor_formatter(NullFormatter())
        elif key == "roc_auc": a.set_ylim(0.5, 1.0)
        elif key == "norm_pr_auc": a.set_ylim(0, 0.5)
        else: a.set_ylim(0, 100)
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(RULE)
        a.tick_params(colors=MUTED, labelsize=9)
        a.set_xticks(THR); a.set_xticklabels([f"{t:.2f}".rstrip("0").rstrip(".") if t not in (0.1, 0.2, 0.3, 0.4, 0.5) else f"{t:.1f}" for t in THR], fontsize=8)
        if i == 0: a.set_title(st, fontsize=11, color=INK, loc="left", fontweight="bold")
        if j == 0: a.set_ylabel(ylab, color=BODY, fontsize=9.5)
        if i == 3: a.set_xlabel("peak threshold T", color=BODY, fontsize=10)
    if all(math.isnan(num(by[(st, tm, hd, t)], "p50")) for tm in ("C", "N") for hd in HEADS for t in THR if (st, tm, hd, t) in by):
        axs[1][j].text(0.5, 0.5, "fewer than half the ends\nhave a peak at any T", transform=axs[1][j].transAxes,
                       ha="center", va="center", fontsize=9.5, color=MUTED)
h, l = axs[0][0].get_legend_handles_labels()
fig.suptitle("Peak calling against the threshold: reach, cost, and ranking over peaks", x=0.055, y=0.995, ha="left",
             fontsize=15, color=INK, fontweight="bold")
fig.legend(h, l, loc="upper left", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.045, 0.975), handlelength=2.6)
fig.tight_layout(rect=(0, 0.025, 1, 0.945))
fig.text(0.055, 0.005, "A missing point in the second row: fewer than half of the ends have any peak at that threshold, "
         "so 50% recall is never reached.  ROC and PR over peaks compare the heads at one T, not one T against another.\n"
         "Cross-validated: each end scored only by the fold that held it out; never seen: all 50 models.",
         fontsize=9, color=MUTED, ha="left", va="bottom")
p = OUT / "peak_threshold_sweep.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); plt.close(fig); print("wrote", p)
