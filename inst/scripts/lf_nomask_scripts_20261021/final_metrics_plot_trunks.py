#!/usr/bin/env python3
"""One arm against another (ARMS, default curated vs ESM-C): a metrics table and
recall-vs-rank small multiples, from metrics_final.R's two CSVs.

Colour = trunk (categorical slots 1/2 of the reference palette, validated on
white: CVD dE 24.7, normal 33.6, contrast 4.4 / 3.2); line style = head, so
identity never rests on colour alone. Text stays in ink, never series colour.
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
## ARMS picks which two arms the figure compares; the default is the original
## pair named by ARMS; the title and the per-cell delta follow it.
ARMS = os.environ.get("ARMS", "flatten,esmonly36").split(",")
SUF = "_" + "_vs_".join(ARMS)
## Arm-keyed labels. These are LOOKUPS WITH A FALLBACK, not plain dicts: a new
## arm used to raise KeyError here and kill the figure after the other five had
## already been written, which is a silly way to lose a run. An unknown arm now
## gets a readable name derived from its own token and the next validated
## colour slot.
_COL = {"flatten": "#2a78d6", "attention": "#eb6834", "esmonly36": "#1b9e77",
        "curated_wAAid": "#eb6834"}
_SLOTS = ["#eb6834", "#1b9e77", "#9467bd", "#7f5539"]
DASH = {"cnn": "-", "xgb_pm": (0, (4, 2.5))}
HEAD = {"cnn": "CNN head", "xgb_pm": "xgb head"}
## short enough for the trunk column, which sits 1.55 units before the head
## column; the title carries the full names.
_TRUNK = {"flatten": "curated", "attention": "attention", "esmonly36": "ESM-C only",
          "curated_wAAid": "curated+AAid"}
_TRUNK_LONG = {"flatten": "curated (26 hand-built channels)", "attention": "attention trunk",
               "esmonly36": "ESM-C only (36 PCA channels)",
               "curated_wAAid": "curated + AA one-hot (47 channels)"}


def _name(arm):
    return _TRUNK.get(arm, arm.replace("_", " "))


def _long(arm):
    return _TRUNK_LONG.get(arm, _name(arm))


class _Fallback(dict):
    """dict that never raises on an unknown arm -- see the note above."""
    def __init__(self, base, default):
        super().__init__(base); self._d = default
    def __missing__(self, k):
        return self._d(k)


TRUNK = _Fallback(_TRUNK, _name)
TRUNK_LONG = _Fallback(_TRUNK_LONG, _long)
COL = _Fallback(_COL, lambda k: _SLOTS[(abs(hash(k)) % len(_SLOTS))])
SETS = ["cross-validated", "never seen: peptides", "never seen: large products"]
SETNOTE = {"cross-validated": "each end scored only by the fold that held it out",
           "never seen: peptides": "products < 60 aa; no fold trained on the gene (all 50)",
           "never seen: large products": "products of 60 aa and over (all 50)"}

rows = list(csv.DictReader(open(D / "final_metrics.csv")))
curv = list(csv.DictReader(open(D / "final_recall_rank.csv")))
num = lambda r, k: float(r[k]) if r[k] not in ("", "NA") else float("nan")
by = {(r["set"], r["term"], r["trunk"], r["head"]): r for r in rows}

# ---------------------------------------------------------------- the table
COLS = [("Ends", "n", "{:.0f}", None), ("ROC AUC", "roc_auc", "{:.4f}", +1),
        ("Norm. PR-AUC", "norm_pr_auc", "{:.3f}", +1), ("AP", "ap", "{:.4f}", +1),
        ("Anchors for\n50% recall", "n50", "{:,.0f}", -1), ("Anchors for\n80% recall", "n80", "{:,.0f}", -1)]
W = 13.2; L, R = 0.55, W - 0.55
ROW, SEP = 0.34, 0.30
nrows = len(SETS) * 2 * 4
H = 1.55 + 0.55 + nrows * ROW + len(SETS) * 2 * SEP + 0.95
fig = plt.figure(figsize=(W, H), dpi=180); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10.5, c=BODY, w="normal", ha="left": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
T(L, H - 0.45, f"Window models: {TRUNK_LONG[ARMS[0]]} vs {TRUNK_LONG[ARMS[1]]}", 17, INK, "bold")
T(L, H - 0.85, "5-fold CV x 10 seeds per terminus, scored over every mature anchor of the proteome "
  "(~3M per terminus).  Other annotated ends are left out of the pool.", 10.5, MUTED)
X = [L, L + 2.4, L + 2.85, L + 4.4]                       # set, term, trunk, head
XC = [L + 5.75, L + 7.0, L + 8.3, L + 9.35, L + 10.85, L + 12.1]
y = H - 1.45
for (lab, _, _, _), x in zip(COLS, XC):
    T(x, y, lab, 9.5, MUTED, "bold", "right")
T(X[0], y, "Set", 9.5, MUTED, "bold"); T(X[1], y, "Term", 9.5, MUTED, "bold")
T(X[2], y, "Trunk", 9.5, MUTED, "bold"); T(X[3], y, "Head", 9.5, MUTED, "bold")
y -= 0.22; ax.plot([L, R], [y, y], color=INK, lw=0.8); y -= 0.32
for st in SETS:
    for term in ("C", "N"):
        top = y
        for hd in ("cnn", "xgb_pm"):
            pair = {tr: by.get((st, term, tr, hd)) for tr in ARMS}
            for tr in ARMS:
                r = pair[tr]
                if r is None: continue
                if hd == "cnn" and tr == "flatten":
                    T(X[0], y, st, 10, INK, "bold"); T(X[1], y, term, 10, INK, "bold")
                # a small swatch carries the trunk's identity; the text stays ink
                ax.plot([X[2], X[2] + 0.28], [y, y], color=COL[tr], lw=2.2,
                        ls=DASH[hd], solid_capstyle="butt")
                T(X[2] + 0.38, y, TRUNK[tr].replace(" trunk", ""), 10, BODY)
                T(X[3], y, HEAD[hd].replace(" head", ""), 10, BODY)
                for (lab, key, fmt, better), x in zip(COLS, XC):
                    v = num(r, key)
                    o = pair[ARMS[1] if tr == ARMS[0] else ARMS[0]]
                    win = False
                    if better is not None and o is not None and not math.isnan(v):
                        ov = num(o, key)
                        win = (v - ov) * better > 0 and fmt.format(v) != fmt.format(ov)
                    T(x, y, "-" if math.isnan(v) else fmt.format(v), 10,
                      INK if win else BODY, "bold" if win else "normal", "right")
                y -= ROW
        ax.plot([L, R], [y + ROW / 2 - 0.02, y + ROW / 2 - 0.02], color=RULE, lw=0.6)
        y -= SEP - 0.04
T(L, y + 0.05, "Bold marks the better trunk for that head and row (higher AUC / PR; fewer anchors "
  "to reach the recall).  Normalised PR-AUC is AP re-derived at a common 1% prevalence", 9.5, MUTED)
T(L, y - 0.22, "by reweighting the negatives along the same ranking, so the three sets compare.  "
  "Swatch: colour = trunk, solid = CNN head, dashed = xgb head.", 9.5, MUTED)
p = OUT / f"final_model_metrics_trunks{SUF}.png"; fig.savefig(p, facecolor="white"); plt.close(fig)
print("wrote", p)

# ------------------------------------------------- recall vs rank, 3 x 2 panels
fig, axs = plt.subplots(len(SETS), 2, figsize=(11.6, 10.4), dpi=180, sharex=True, sharey=True)
fig.patch.set_facecolor("white")
for i, st in enumerate(SETS):
    for j, term in enumerate(("C", "N")):
        a = axs[i][j]
        for tr in ARMS:
            for hd in ("cnn", "xgb_pm"):
                pts = [(float(c["N"]), float(c["recall"])) for c in curv
                       if c["set"] == st and c["term"] == term and c["trunk"] == tr and c["head"] == hd]
                if not pts: continue
                xs, ys = zip(*pts)
                a.plot(xs, ys, color=COL[tr], lw=1.9, ls=DASH[hd],
                       label=f"{TRUNK[tr]}, {HEAD[hd]}")
        a.set_xscale("log"); a.set_xlim(1, 3e6); a.set_ylim(0, 1.0)
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True)
        for s in ("top", "right"): a.spines[s].set_visible(False)
        for s in ("left", "bottom"): a.spines[s].set_color(RULE)
        a.tick_params(colors=MUTED, labelsize=9)
        r0 = by.get((st, term, "flatten", "cnn"))
        n = int(float(r0["n"])) if r0 else 0
        a.set_title(f"{term} terminus  |  {st}  (n = {n})", fontsize=10.5, color=INK, loc="left")
        if j == 0: a.set_ylabel("recall", color=BODY, fontsize=10)
        if i == len(SETS) - 1: a.set_xlabel("anchors examined, best-scored first (log)", color=BODY, fontsize=10)
# No in-panel direct labels: the four curves coincide over most of their range
# (both xgb heads almost everywhere), so any placement sits on the data. The
# legend carries both encodings (colour = trunk, dash = head); the metrics table
# is the table view.
h, l = axs[0][0].get_legend_handles_labels()
fig.suptitle("Recall vs rank: how many anchors must be examined to recover the ends",
             x=0.055, y=0.995, ha="left", fontsize=15, color=INK, fontweight="bold")
fig.legend(h, l, loc="upper left", ncol=4, frameon=False, fontsize=9.5,
           bbox_to_anchor=(0.045, 0.968), handlelength=2.6)
fig.tight_layout(rect=(0, 0, 1, 0.935))
p = OUT / f"final_recall_rank_trunks{SUF}.png"; fig.savefig(p, facecolor="white", bbox_inches="tight"); plt.close(fig)
print("wrote", p)
