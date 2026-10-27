#!/usr/bin/env python3
"""Segment-level confusion table and ROC curves, from segments_roc.R (v2).

Threshold T as a classifier: runs of anchors >= T are positive segments, runs
< T negative segments. Ends are counted per end -- TP if a positive segment's
maximum is within TOL residues, FN otherwise (TP + FN = ends); positive segments
whose maximum is near no end are FP, negative segments holding no real end TN.
Per T the single-point ROC AUC (1 + TPR - FPR) / 2 (balanced accuracy).
Colour = head (validated categorical slots 1/2); text stays in ink.
"""
import csv, math, os, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
HALO = [pe.withStroke(linewidth=2.6, foreground="white")]
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = pathlib.Path.home() / "AF2_analysis"; D = A2 / "lf_headcmp"
OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else A2
SUF = os.environ.get("SUF", "_tol4")     # matches segments_roc.R's own default (TOL 4)
TOL = int(os.environ.get("TOL", "4"))   # only labels the caption; the CSV already encodes it
INK, BODY, MUTED, RULE, GRID = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9", "#E6E5E1"
HEADS = ["cnn", "xgb_pm"]
COL = {"cnn": "#2a78d6", "xgb_pm": "#eb6834"}
NAME = {"cnn": "linear head", "xgb_pm": "xgb head"}
SETS = ["cross-validated", "never seen: peptides", "never seen: large products"]
TT = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
DEF = ("Positive segment = a run of anchors with score >= T; negative segment = a run < T.  TP = an end with a positive "
       f"segment's maximum within {TOL} residues; FN = every other end (TP + FN = the ends); FP = a positive segment whose maximum "
       "is near no end; TN = a negative segment holding no real end.")

rows = list(csv.DictReader(open(D / f"segment_confusion{SUF}.csv")))
aucs = {(r["set"], r["term"], r["head"]): r for r in csv.DictReader(open(D / f"segment_roc_auc{SUF}.csv"))}
num = lambda r, k: float(r[k]) if r[k] not in ("", "NA", "NaN") else float("nan")
by = {(r["set"], r["term"], r["head"], round(float(r["thr"]), 3)): r for r in rows}
curve = lambda st, term, hd: sorted([r for r in rows if r["set"] == st and r["term"] == term and r["head"] == hd],
                                    key=lambda r: -float(r["thr"]))

# ---------------------------------------------------------------- the table
MET = [("TP", "TP", "{:,.0f}"), ("FN", "FN", "{:,.0f}"), ("FP", "FP", "{:,.0f}"), ("TN", "TN", "{:,.0f}"),
       ("TPR", "TPR", "{:.3f}"), ("FPR", "FPR", "{:.3f}"), ("Precision", "precision", "{:.3f}"),
       ("ROC AUC\nat T", "auc_T", "{:.3f}")]
W = 21.6; L, R = 0.55, W - 0.55; ROW, SEP = 0.27, 0.42
nrow = len(SETS) * 2 * len(TT)
H = 2.45 + nrow * ROW + len(SETS) * 2 * SEP + 1.95
fig = plt.figure(figsize=(W, H), dpi=150); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10, c=BODY, w="normal", ha="left": ax.text(x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
T(L, H - 0.42, "Threshold as a classifier over segments", 17, INK, "bold")
T(L, H - 0.78, DEF, 9.6, MUTED)
X0 = [L, L + 2.45, L + 2.95]
BLK = {"cnn": L + 3.55, "xgb_pm": L + 12.35}
OFF = [0.6, 1.45, 2.55, 3.6, 4.55, 5.5, 6.7, 7.85]
yh = H - 1.38
for hd in HEADS:
    x0 = BLK[hd]
    ax.plot([x0, x0 + 0.3], [yh + 0.05, yh + 0.05], color=COL[hd], lw=2.4, solid_capstyle="butt")
    T(x0 + 0.42, yh + 0.05, NAME[hd], 10.5, INK, "bold")
    ax.plot([x0, x0 + OFF[-1]], [yh - 0.17, yh - 0.17], color=RULE, lw=0.8)
    ## unit groups: TP / FN count ends, FP / TN count segments
    T(x0 + (OFF[0] + OFF[1]) / 2 - 0.2, yh - 0.36, "ends", 8.6, MUTED, "normal", "center")
    T(x0 + (OFF[2] + OFF[3]) / 2 - 0.25, yh - 0.36, "segments", 8.6, MUTED, "normal", "center")
y = yh - 0.7
T(X0[0], y, "Set", 9.3, MUTED, "bold"); T(X0[1], y, "Term", 9.3, MUTED, "bold"); T(X0[2], y, "T", 9.3, MUTED, "bold")
for hd in HEADS:
    for (lab, _, _), o in zip(MET, OFF): T(BLK[hd] + o, y, lab, 9.3, MUTED, "bold", "right")
y -= 0.33; ax.plot([L, R], [y, y], color=INK, lw=0.8); y -= 0.28
for st in SETS:
    for term in ("C", "N"):
        T(X0[0], y, st, 9.8, INK, "bold"); T(X0[1], y, term, 9.8, INK, "bold")
        for hd in HEADS:
            a = aucs.get((st, term, hd))
            if a: T(BLK[hd] + OFF[-1], y, f"best T over the full grid  {float(a['thr_auc_T_max']):g}:  "
                    f"ROC AUC at T {float(a['auc_T_max']):.3f}", 9.3, INK, "bold", "right")
        y -= ROW
        best = {hd: max(num(by[(st, term, hd, round(t, 3))], "auc_T") for t in TT if (st, term, hd, round(t, 3)) in by)
                for hd in HEADS}
        for t in TT:
            T(X0[2], y, f"{t:.2f}", 9.8, BODY)
            for hd in HEADS:
                r = by.get((st, term, hd, round(t, 3)))
                if r is None: continue
                for (lab, key, fmt), o in zip(MET, OFF):
                    v = num(r, key)
                    win = key == "auc_T" and v == best[hd]
                    T(BLK[hd] + o, y, "-" if math.isnan(v) else fmt.format(v), 9.8,
                      INK if win else BODY, "bold" if win else "normal", "right")
            y -= ROW
        ax.plot([L, R], [y + ROW / 2 - 0.02, y + ROW / 2 - 0.02], color=RULE, lw=0.6); y -= SEP - 0.16
T(L, y + 0.02, "TPR = TP / (TP + FN), over ends; FPR = FP / (FP + TN), over segments; precision = TP segments / positive segments.  "
  "ROC AUC at T = (1 + TPR - FPR) / 2, the area under the curve through (0,0), that point and (1,1); bold = each head's best T in this table.", 9, MUTED)
T(L, y - 0.23, "No area over T is given: FPR cannot pass ~0.5 at any T > 0 (positive and negative segments alternate, so FP ~ TN once peaks "
  "are everywhere) and TPR is not monotone in T, so that area depends mostly on how the gap to (1,1) is closed.", 9, MUTED)
T(L, y - 0.46, "An end inside a positive segment whose maximum lands more than 2 residues away is an FN, and that segment an FP.  "
  "A negative segment holding an evaluation end is neither TN nor FN.  Segments touching an annotated end outside the", 9, MUTED)
T(L, y - 0.69, "evaluation set (other folds' train ends, ligand_list, UniProt peptide / propeptide) are held out of every count.  "
  "Cross-validated: per fold, each end scored only by the fold that held it out; counts summed over the 5 folds.  "
  "Never seen: all 50 models.", 9, MUTED)
p = OUT / f"segment_confusion_table{SUF}.png"; fig.savefig(p, facecolor="white"); plt.close(fig); print("wrote", p)

# ------------------------------------------------------------ ROC curves
fig, axs = plt.subplots(2, len(SETS), figsize=(13.2, 8.6), dpi=180, sharex=True, sharey=True)
fig.patch.set_facecolor("white")
for j, st in enumerate(SETS):
    for i, term in enumerate(("C", "N")):
        a = axs[i][j]
        a.plot([0, 1], [0, 1], color=GRID, lw=1, ls=(0, (3, 3)), zorder=0)
        for hd in HEADS:
            pts = [(float(r["thr"]), num(r, "FPR"), num(r, "TPR")) for r in curve(st, term, hd)]
            pts = [p for p in pts if not (math.isnan(p[1]) or math.isnan(p[2]))]
            p0 = [p for p in pts if p[0] == 0]; pts = [p for p in pts if p[0] > 0]
            xs = [0] + [p[1] for p in pts]; ys = [0] + [p[2] for p in pts]
            a.plot(xs, ys, color=COL[hd], lw=1.9, label=NAME[hd])
            ## T = 0 is degenerate (one segment per precursor): dotted, not part of the curve
            if p0:
                a.plot([xs[-1], p0[0][1]], [ys[-1], p0[0][2]], color=COL[hd], lw=1.1, ls=(0, (1.5, 2)))
                pts = pts + p0
            mk = [(p[1], p[2], p[0]) for p in pts if round(p[0], 3) in TT or p[0] == 0]
            a.scatter([m[0] for m in mk], [m[1] for m in mk], s=16, facecolor="white", edgecolor=COL[hd], lw=1.1, zorder=3)
            if hd == "cnn":
                for fx, ty, t in mk:
                    if t in (0.1, 0.5):
                        a.annotate(f"T={t:g}", (fx, ty), xytext=(5, -13), textcoords="offset points", fontsize=7.8,
                                   color=INK, path_effects=HALO, zorder=5)
                    elif t == 0:
                        a.annotate("T=0", (fx, ty), xytext=(-6, -12), textcoords="offset points", fontsize=7.8, color=INK,
                                   ha="right", path_effects=HALO, zorder=5)
        a.set_xlim(0, 1); a.set_ylim(0, 1.0)
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(RULE)
        a.tick_params(colors=MUTED, labelsize=9)
        a.set_title(f"{term} terminus  |  {st}", fontsize=10.5, color=INK, loc="left")
        a.legend(loc="upper left", frameon=False, fontsize=8.5)
        if j == 0: a.set_ylabel("TPR (ends)", color=BODY, fontsize=10)
        if i == 1: a.set_xlabel("FPR (segments)", color=BODY, fontsize=10)
fig.suptitle("ROC over thresholds: ends caught against negative segments called", x=0.055, y=0.995, ha="left",
             fontsize=15, color=INK, fontweight="bold")
fig.text(0.055, 0.945, "Each point is one threshold T, from T = 1 at (0,0) down to 0.025 (marked: 0.10-0.50 in steps of 0.05; T=0.1 and 0.5 labelled on the linear head).\nFPR stalls near 0.5 "
         "for every T > 0: positive and negative segments alternate, so FP ~ TN.  Only T = 0 (dotted; one segment per precursor) "
         "reaches FPR = 1,\nso no area under these points is reported.", fontsize=9.2, color=MUTED, ha="left", va="top")
fig.tight_layout(rect=(0, 0, 1, 0.9))
p = OUT / f"segment_roc{SUF}.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); plt.close(fig); print("wrote", p)

# ------------------------------------- the threshold scan, as 2D plots (fine grid)
TLS = {"C": "-", "N": (0, (4, 2.5))}
PANELS = [("TPR", "TPR (ends)"), ("FPR", "FPR (segments)"), ("precision", "precision (segments, log)"),
          ("auc_T", "ROC AUC at T\n(1 + TPR - FPR) / 2")]
fig, axs = plt.subplots(len(PANELS), len(SETS), figsize=(13.2, 13.8), dpi=180, sharex=True)
fig.patch.set_facecolor("white")
for j, st in enumerate(SETS):
    for i, (key, ylab) in enumerate(PANELS):
        a = axs[i][j]
        if key == "auc_T": a.axhline(0.5, color=MUTED, lw=0.9, ls=(0, (3, 3)), zorder=1)
        for hd in HEADS:
            for term in ("C", "N"):
                pts = sorted([(float(r["thr"]), num(r, key)) for r in curve(st, term, hd) if float(r["thr"]) > 0])
                pts = [p for p in pts if not math.isnan(p[1]) and (key != "precision" or p[1] > 0)]
                a.plot([p[0] for p in pts], [p[1] for p in pts], color=COL[hd], ls=TLS[term], lw=1.8,
                       label=f"{NAME[hd]}, {term} terminus")
                if key == "auc_T" and pts:
                    tb, vb = max(pts, key=lambda p: p[1])
                    a.scatter([tb], [vb], s=22, facecolor="white", edgecolor=COL[hd], lw=1.3, zorder=4)
        if key == "precision":
            a.set_yscale("log"); a.set_ylim(2e-4, 0.3)
            a.yaxis.set_major_locator(FixedLocator([0.001, 0.003, 0.01, 0.03, 0.1]))
            a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            a.yaxis.set_minor_formatter(NullFormatter())
        elif key == "auc_T":
            a.set_ylim(0.4, 0.8)
        else:
            a.set_ylim(0, 1.0)
        a.set_xlim(0, 1)
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(RULE)
        a.tick_params(colors=MUTED, labelsize=9)
        if i == 0: a.set_title(st, fontsize=11, color=INK, loc="left", fontweight="bold")
        if j == 0: a.set_ylabel(ylab, color=BODY, fontsize=10)
        if i == len(PANELS) - 1: a.set_xlabel("threshold T", color=BODY, fontsize=10)
h, l = axs[0][0].get_legend_handles_labels()
fig.suptitle("Threshold as a classifier over segments: TPR, FPR, precision and ROC AUC against T", x=0.055, y=0.995,
             ha="left", fontsize=15, color=INK, fontweight="bold")
fig.legend(h, l, loc="upper left", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.045, 0.975), handlelength=2.6)
fig.tight_layout(rect=(0, 0.04, 1, 0.958))
fig.text(0.055, 0.004, DEF.replace("; FP = ", ";\nFP = ") + "  T from 0.025 to 1 (T = 0 is one segment per precursor).\nBottom row: open marker = each curve's best T; "
         "dashed line = 0.5, a classifier whose TPR equals its FPR.",
         fontsize=8.6, color=MUTED, ha="left", va="bottom")
p = OUT / f"segment_threshold_scan{SUF}.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); plt.close(fig); print("wrote", p)
