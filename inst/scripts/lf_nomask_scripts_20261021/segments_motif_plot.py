#!/usr/bin/env python3
"""The per-motif segment table at one threshold, from segments_motif.R.

Rows: set x terminus x cleavage motif; the two heads side by side.
TP / FN count the evaluation ENDS whose anchor carries that motif; FP / TN the
SEGMENTS overlapping at least one anchor of it. ROC AUC at T = (1+TPR-FPR)/2.
Colour = head (validated categorical slots 1/2); text stays in ink.
"""
import csv, math, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = os.path.expanduser("~/AF2_analysis"); D = A2 + "/lf_headcmp"
THR = sys.argv[1] if len(sys.argv) > 1 else "0.3"
TOL = int(sys.argv[2]) if len(sys.argv) > 2 else 2        # +/- residues for "on the end"
ARM = os.environ.get("ARM", "")          # "" = curated; "_esmonly36" = the ESM-C arm
SFX = ARM + ("" if TOL == 2 else f"_tol{TOL}")
OUT = A2
INK, BODY, MUTED, RULE = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9"
HEADS = ["cnn", "xgb_pm"]
COL = {"cnn": "#2a78d6", "xgb_pm": "#eb6834"}
NAME = {"cnn": "linear head", "xgb_pm": "xgb head"}
SETS = ["cross-validated", "never seen: peptides", "never seen: large products"]
CLASSES = ["terminus", "dibasic", "G + dibasic", "monobasic", "other"]

rows = list(csv.DictReader(open(f"{D}/segment_motif_T{THR}{SFX}.csv")))
num = lambda r, k: float(r[k]) if r[k] not in ("", "NA", "NaN") else float("nan")
by = {(r["set"], r["term"], r["head"], r["motif"]): r for r in rows}

MET = [("Ends", "ends", "{:,.0f}"), ("TP", "TP", "{:,.0f}"), ("FN", "FN", "{:,.0f}"),
       ("FP", "FP", "{:,.0f}"), ("TN", "TN", "{:,.0f}"),
       ("TPR", "TPR", "{:.3f}"), ("FPR", "FPR", "{:.3f}"), ("ROC AUC\nat T", "auc_T", "{:.3f}")]
W = 21.0; L, R = 0.55, W - 0.55; ROW, SEP = 0.28, 0.42
## only the motif classes that actually occur (G + dibasic is C-only)
nrow = sum(1 for st in SETS for term in ("C", "N") for m in CLASSES
           if any((st, term, hd, m) in by for hd in HEADS))
H = 2.5 + (nrow + len(SETS) * 2) * ROW + len(SETS) * 2 * SEP + 1.0
fig = plt.figure(figsize=(W, H), dpi=150); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10, c=BODY, w="normal", ha="left": ax.text(x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
T(L, H - 0.42, f"Segments by cleavage motif, at T = {THR}   (on the end = within {TOL} residues)", 17, INK, "bold")
T(L, H - 0.78, "Motif of an anchor (lf_pepend_motif_vec, the rule the knowns are labelled with): terminus beats dibasic beats "
  "G + dibasic beats monobasic, else other.  C looks forward from the anchor (the cut is after it), N looks back; "
  "G + dibasic is a C-only class.", 9.4, MUTED)
X0 = [L, L + 2.45, L + 2.95]
BLK = {"cnn": L + 4.55, "xgb_pm": L + 12.5}
OFF = [0.75, 1.5, 2.25, 3.3, 4.4, 5.3, 6.2, 7.4]
yh = H - 1.42
for hd in HEADS:
    x0 = BLK[hd]
    ax.plot([x0, x0 + 0.3], [yh + 0.05, yh + 0.05], color=COL[hd], lw=2.4, solid_capstyle="butt")
    T(x0 + 0.42, yh + 0.05, NAME[hd], 10.5, INK, "bold")
    ax.plot([x0, x0 + OFF[-1]], [yh - 0.17, yh - 0.17], color=RULE, lw=0.8)
    T(x0 + (OFF[0] + OFF[2]) / 2 - 0.25, yh - 0.36, "ends", 8.6, MUTED, "normal", "center")
    T(x0 + (OFF[3] + OFF[4]) / 2 - 0.3, yh - 0.36, "segments", 8.6, MUTED, "normal", "center")
y = yh - 0.7
T(X0[0], y, "Set", 9.3, MUTED, "bold"); T(X0[1], y, "Term", 9.3, MUTED, "bold"); T(X0[2], y, "Motif", 9.3, MUTED, "bold")
for hd in HEADS:
    for (lab, _, _), o in zip(MET, OFF): T(BLK[hd] + o, y, lab, 9.3, MUTED, "bold", "right")
y -= 0.33; ax.plot([L, R], [y, y], color=INK, lw=0.8); y -= 0.28
for st in SETS:
    for term in ("C", "N"):
        T(X0[0], y, st, 9.8, INK, "bold"); T(X0[1], y, term, 9.8, INK, "bold")
        ## the whole set, as the un-split table reports it
        tot = {}
        for hd in HEADS:
            got = [by[(st, term, hd, m)] for m in CLASSES if (st, term, hd, m) in by]
            s = {k: sum(num(r, k) for r in got) for k in ("ends", "TP", "FN", "FP", "TN")}
            tpr = s["TP"] / (s["TP"] + s["FN"]) if s["TP"] + s["FN"] else float("nan")
            fpr = s["FP"] / (s["FP"] + s["TN"]) if s["FP"] + s["TN"] else float("nan")
            tot[hd] = dict(s, TPR=tpr, FPR=fpr, auc_T=(1 + tpr - fpr) / 2)
            T(BLK[hd] + OFF[-1], y, f"pooled over classes  {tot[hd]['auc_T']:.3f}", 9.3, INK, "bold", "right")
        y -= ROW
        best = {hd: max((num(by[(st, term, hd, m)], "auc_T") for m in CLASSES
                         if (st, term, hd, m) in by and num(by[(st, term, hd, m)], "ends") >= 5),
                        default=float("nan")) for hd in HEADS}
        for m in CLASSES:
            if not any((st, term, hd, m) in by for hd in HEADS): continue
            T(X0[2], y, m, 9.8, BODY)
            for hd in HEADS:
                r = by.get((st, term, hd, m))
                if r is None: continue
                thin = num(r, "ends") < 5          # too few ends for TPR to mean much
                for (lab, key, fmt), o in zip(MET, OFF):
                    v = num(r, key)
                    win = key == "auc_T" and not thin and v == best[hd]
                    c = MUTED if (thin and key in ("TPR", "auc_T")) else (INK if win else BODY)
                    T(BLK[hd] + o, y, "-" if math.isnan(v) else fmt.format(v), 9.8, c,
                      "bold" if win else "normal", "right")
            y -= ROW
        ax.plot([L, R], [y + ROW / 2 - 0.02, y + ROW / 2 - 0.02], color=RULE, lw=0.6); y -= SEP - 0.16
T(L, y + 0.02, "TP / FN count the evaluation ENDS whose own anchor carries that motif, so TP + FN = that row's ends; a TP is an end with a "
  f"positive segment's maximum within {TOL} residues.  FP / TN count the SEGMENTS", 9, MUTED)
T(L, y - 0.23, "overlapping at least one anchor of that motif and holding no real end: FP when the run is >= T and peaks near no end, TN when "
  "the run is < T.  A segment spanning two motif classes is counted in both rows", 9, MUTED)
T(L, y - 0.46, "(a long negative run usually overlaps several).  'Pooled over classes' adds those per-class counts up, so it shares the true TPR "
  "(each end has one motif) but its FPR counts a segment once per class it", 9, MUTED)
T(L, y - 0.69, "overlaps -- it is NOT the un-split figure (that is segment_confusion.csv).  ROC AUC at T = (1 + TPR - FPR) / 2, the area under the "
  "curve through (0,0), this threshold's one point and (1,1).  Bold = each", 9, MUTED)
T(L, y - 0.92, "head's best motif within a set and terminus, among rows with at least 5 ends; a row with fewer is greyed, not ranked.  Segments touching "
  "an annotated end outside the evaluation set (other folds' train ends,", 9, MUTED)
T(L, y - 1.15, "ligand_list, UniProt peptide / propeptide) are held out of every count.  Cross-validated: per fold, each end scored only by the fold that "
  "held it out, counts summed over the 5 folds.  Never seen: all 50 models,", 9, MUTED)
T(L, y - 1.38, "and both never-seen sets hold out the same union of annotated ends, so they share one negative pool: their FP / TN columns are "
  "identical and only TPR separates them.", 9, MUTED)
p = f"{OUT}/segment_motif_table_T{THR}{SFX}.png"; fig.savefig(p, facecolor="white"); plt.close(fig); print("wrote", p)
