#!/usr/bin/env python3
"""How close the xgb track's nearest peak lands to the true end."""
import csv, pathlib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
INK, BODY, MUTED, ACCENT, RULE = "#111418", "#3F4A56", "#6B7784", "#197EC0", "#C9D1D9"
S = pathlib.Path(__file__).parent
rows = list(csv.DictReader(open(S/"peak_dist.csv")))
by = {(r["set"], r["term"], r["scorer"]): r for r in rows}
SETS = [("cv", "Cross-validated", "178 trained ends -- IN-SAMPLE for both rules,\nshown only as a ceiling"),
        ("unseen_pep", "Never seen: peptides", "products under 60 aa, no fold ever\ntrained on their gene"),
        ("unseen_big", "Never seen: large products", "cleavage products of 60 aa and over")]
SC = {"all50": "ensemble mean (all 50)", "maxfold": "max of the 5 fold means"}

W = 13.2; L, R = 0.6, W - 0.6
ROW, SEP = 0.40, 0.26
H = 1.75 + 0.60 + 0.60 + len(SETS) * (4 * ROW + SEP) + 0.64 + 0.56
fig = plt.figure(figsize=(W, H), dpi=220); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_visible(False)
ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=11, c=BODY, w="normal", ha="left": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
rule = lambda y, c=RULE, lw=0.9: ax.plot([L, R], [y, y], color=c, lw=lw, solid_capstyle="butt")

T(L, H-0.62, "How close is the nearest peak to the true end?", 20, INK, "bold")
T(L, H-1.05, "xgb window head only.  A residue is a PEAK when its score is at least the "
  "maximum over the surrounding +/-5 residues.", 11, MUTED)
T(L, H-1.35, "Each column is the share of true ends whose nearest peak falls within that "
  "distance; a protein carries a median of 7-14 peaks.", 11, MUTED)

CX =    [L,    4.15,  5.35,  7.95,  8.95,  9.90, 10.85,  R]
ALIGN = ["left","left","left","right","right","right","right","right"]
HEADS = ["Set","Term","Scoring rule","Exact","within 1","within 2","within 5","Median |d|"]
y = H - 1.95
for cx, h, al in zip(CX, HEADS, ALIGN):
    T(cx, y, h.upper(), 9.5, MUTED, "bold", ha=al)
rule(y - 0.19, INK, 1.1)

y -= 0.58
for sk, label, blurb in SETS:
    T(CX[0], y + 0.02, label, 12, INK, "bold")
    T(CX[0], y - 0.36, blurb, 9, MUTED)
    for term in ("N", "C"):
        for sc in ("all50", "maxfold"):
            r = by[(sk, term, sc)]; o = by[(sk, term, "maxfold" if sc == "all50" else "all50")]
            T(CX[1], y, f"{term} terminus", 11, BODY)
            T(CX[2], y, SC[sc], 11, BODY)
            for cx, key, fmt in ((CX[3], "exact", "{:.0%}"), (CX[4], "w1", "{:.0%}"),
                                 (CX[5], "w2", "{:.0%}"), (CX[6], "w5", "{:.0%}"),
                                 (CX[7], "med", "{:.0f}")):
                mine, theirs = fmt.format(float(r[key])), fmt.format(float(o[key]))
                better = float(r[key]) < float(o[key]) if key == "med" else float(r[key]) > float(o[key])
                best = mine != theirs and better
                T(cx, y, mine, 11.5, ACCENT if best else BODY,
                  "bold" if best else "normal", ha="right")
            y -= ROW
    rule(y + (ROW - SEP) / 2); y -= SEP

T(L, y + 0.04, "Blue marks the better rule on that row.  `Exact` means the peak sits ON the "
  "true end; `within 2` is the tolerance the window labels already use.", 10, MUTED)
T(L, y - 0.28, "n = 92 / 86 cross-validated, 45 / 45 never-seen peptides, 294 / 317 large "
  "products (N / C).", 10, MUTED)
out = pathlib.Path.home()/"Desktop/lf_presentation/peak_distance.png"
fig.savefig(out, facecolor="white")
print("wrote", out, f"({out.stat().st_size//1024} KB)")
