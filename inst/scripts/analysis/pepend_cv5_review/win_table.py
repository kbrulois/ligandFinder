#!/usr/bin/env python3
"""Window-model stats: cross-validated, and on ends no fold ever saw."""
import json, pathlib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
INK, BODY, MUTED, ACCENT, RULE = "#111418", "#3F4A56", "#6B7784", "#197EC0", "#C9D1D9"
S = pathlib.Path(__file__).parent
rows = json.load(open(S/"win_table.json"))
by = {(r["set"], r["term"], r["head"]): r for r in rows}
SETS = [("cv", "Cross-validated",
         "178 trained ends, each scored only by\nthe fold that held it out"),
        ("unseen_pep", "Never seen: peptides",
         "products under 60 aa, no fold ever\ntrained on their gene"),
        ("unseen_big", "Never seen: large products",
         "cleavage products of 60 aa and over")]
HEADNM = {"cnn": "CNN", "xgb_pm": "xgb"}

W = 11.4; L, R = 0.6, W - 0.6
ROW, SEP = 0.40, 0.26
H = 1.75 + 0.30 + 0.60 + len(SETS) * (4 * ROW + SEP) + 0.64 + 0.56
fig = plt.figure(figsize=(W, H), dpi=220); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_visible(False)
ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=11, c=BODY, w="normal", ha="left": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
rule = lambda y, c=RULE, lw=0.9: ax.plot([L, R], [y, y], color=c, lw=lw, solid_capstyle="butt")

T(L, H-0.62, "Window models: where a peptide ends", 20, INK, "bold")
T(L, H-1.05, "Every mature residue of all 5,616 secretome precursors is a candidate anchor, "
  "so the pool is ~2.98 M and prevalence is 6-100 per million.", 11, MUTED)
T(L, H-1.35, "Other annotated ends are REMOVED from the pool, never counted as false "
  "positives: they are neither clean positives nor honest negatives.", 11, MUTED)

CX =    [L,    4.05,  5.15,  6.30,  7.95,  9.35,  R]
ALIGN = ["left","left","left","right","right","right","right"]
HEADS = ["Set","Term","Head","Ends","ROC AUC","AP","Norm. PR-AUC"]
y = H - 1.95
for cx, h, al in zip(CX, HEADS, ALIGN):
    T(cx, y, h.upper(), 9.5, MUTED, "bold", ha=al)
rule(y - 0.19, INK, 1.1)

y -= 0.58
for sk, label, blurb in SETS:
    T(CX[0], y + 0.02, label, 12, INK, "bold")
    T(CX[0], y - 0.36, blurb, 9, MUTED)
    for term in ("N", "C"):
        for hd in ("cnn", "xgb_pm"):
            r = by[(sk, term, hd)]; o = by[(sk, term, "xgb_pm" if hd == "cnn" else "cnn")]
            T(CX[1], y, f"{term} terminus", 11, BODY)
            T(CX[2], y, HEADNM[hd], 11, BODY)
            T(CX[3], y, f"{r['npos']:,}", 11, BODY, ha="right")
            auc = f"{r['auc']:.3f}" + (f" ±{r['auc_sd']:.3f}" if r["auc_sd"] else "")
            for cx, val, fmt, key in ((CX[4], auc, None, "auc"),
                                      (CX[5], None, "{:.4f}", "ap"),
                                      (CX[6], None, "{:.3f}", "ap1")):
                mine = val if val is not None else fmt.format(r[key])
                theirs = (f"{o['auc']:.3f}" + (f" ±{o['auc_sd']:.3f}" if o["auc_sd"] else "")) \
                         if key == "auc" else fmt.format(o[key])
                better = r[key] > o[key]
                best = mine != theirs and better
                T(cx, y, mine, 11.5, ACCENT if best else BODY,
                  "bold" if best else "normal", ha="right")
            y -= ROW
    rule(y + (ROW - SEP) / 2); y -= SEP

T(L, y + 0.04, "Blue marks the better head on that row.  Normalised PR-AUC is the average "
  "precision re-derived at a common 1% prevalence, by reweighting the", 10, MUTED)
T(L, y - 0.28, "negatives along the same ranking.  Raw AP cannot be compared across these "
  "sets: their prevalence spans 5.8 to 105 per million, a 17-fold range.", 10, MUTED)
out = pathlib.Path.home()/"Desktop/lf_presentation/window_model_metrics.png"
fig.savefig(out, facecolor="white")
print("wrote", out, f"({out.stat().st_size//1024} KB)")
