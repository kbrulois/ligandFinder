#!/usr/bin/env python3
"""Presentation table: window-free per-residue models on held-out precursors."""
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

for p in ("/System/Library/Fonts/Supplemental/Arial.ttf",
          "/System/Library/Fonts/Supplemental/Arial Bold.ttf"):
    try: fm.fontManager.addfont(p)
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

INK, BODY, MUTED, ACCENT, RULE = "#111418", "#3F4A56", "#6B7784", "#197EC0", "#C9D1D9"
S = pathlib.Path(__file__).parent
rows = json.load(open(S/"resid_metrics.json"))
NICE = {"pep_pocket": "Pocket residue", "pep_other": "Peptide, not pocket",
        "ct_context": "C-terminal cleavage context", "nt_context": "N-terminal cleavage context"}
ORDER = ["pep_pocket", "ct_context", "nt_context", "pep_other"]
by = {(r["target"], r["model"]): r for r in rows}

W, H = 12.6, 7.15
fig = plt.figure(figsize=(W, H), dpi=220)
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
ax.set_xlim(0, W); ax.set_ylim(0, H)

def T(x, y, s, size=11, color=BODY, weight="normal", ha="left"):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va="center")

def rule(y, x0, x1, c=RULE, lw=0.9):
    ax.plot([x0, x1], [y, y], color=c, lw=lw, solid_capstyle="butt")

L, R = 0.6, W - 0.6
T(L, H - 0.62, "Window-free per-residue models, held out", 20, INK, "bold")
T(L, H - 1.05,
  "30,410 residues from precursors the models never trained on \u2014 whole similarity "
  "groups moved to validation, so no homologue straddles the split.", 11, MUTED)

#        target  model  positives  prev   ROC    AP    lift   recall
CX =    [L,      4.05,  6.00,      7.25,  8.45,  9.45, 10.55, R]
ALIGN = ["left", "left", "right",  "right","right","right","right","right"]
HEADS = ["Target", "Model", "Positives", "Prevalence", "ROC AUC", "AP",
         "AP / prev.", "Recall @ top 1%"]
y = H - 1.75
for cx, h, al in zip(CX, HEADS, ALIGN):
    T(cx, y, h.upper(), 9.5, MUTED, "bold", ha=al)
rule(y - 0.19, L, R, INK, 1.1)

ROW, SEP = 0.42, 0.20
y -= 0.60
for t in ORDER:
    m, x = by[(t, "mlp")], by[(t, "xgb")]
    mid = y - ROW / 2                       # the pair's centre line
    T(CX[0], mid, NICE[t], 12, INK, "bold")
    T(CX[2], mid, f"{m['npos']:,}", 11.5, BODY, ha="right")
    T(CX[3], mid, f"{m['prev']:.2f}%", 11.5, BODY, ha="right")
    for mdl, r in (("MLP", m), ("XGBoost", x)):
        other = x if mdl == "MLP" else m
        T(CX[1], y, mdl, 11.5, BODY)
        for cx, key, fmt in ((CX[4], "auc", "{:.3f}"), (CX[5], "ap", "{:.3f}"),
                             (CX[6], "lift", "{:.0f}\u00d7"), (CX[7], "rec1", "{:.0%}")):
            # compare what is DISPLAYED: 0.3770 vs 0.3766 both print 0.377, and
            # highlighting one of them as the winner reads as a mistake
            mine, theirs = fmt.format(r[key]), fmt.format(other[key])
            best = mine != theirs and r[key] > other[key]
            T(cx, y, mine, 12,
              ACCENT if best else BODY, "bold" if best else "normal", ha="right")
        y -= ROW
    rule(y + (ROW - SEP) / 2, L, R)
    y -= SEP

T(L, y + 0.02,
  "Blue marks the better of the two models on that row.  AP is prevalence-dependent, so "
  "AP / prevalence states the lift over chance.", 10, MUTED)
T(L, y - 0.30,
  "The 45 never-seen peptides carry no labels here \u2014 their precursors were never "
  "docked, so pocket and cleavage-context targets are undefined for them.", 10, MUTED)

out = pathlib.Path.home()/"Desktop/lf_presentation/resid_model_metrics.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, facecolor="white", bbox_inches=None)
print("wrote", out, f"({out.stat().st_size//1024} KB)")
