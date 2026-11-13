#!/usr/bin/env python3
"""The three-way split of the window-free per-residue analysis."""
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
rows = json.load(open(S/"split_rows.json"))

SUB = ["166 precursors — the 66 with docked peptides plus 100 random controls — "
       "contribute 75,551 rows, one per mature residue.",
       "The split unit is a sequence-similarity GROUP, so no homologue straddles it: "
       "117 train groups against 46 validation groups, zero overlap."]
NOTES = ["The early-stop set is carved out of TRAIN, never out of validation — stopping on "
         "the set you then score on is how a leak gets reported as a result.",
         "There is no test set here. The 45 never-seen peptides carry no labels, because their "
         "precursors were never docked, so every number reported is a validation number."]

W = 12.6; L, R = 0.6, W - 0.6
ROW, SEP = 0.46, 0.22
H = 1.75 + 0.30 * (len(SUB) - 1) + 0.60 + len(rows) * ROW + SEP + 0.32 * (len(NOTES) - 1) + 0.56
fig = plt.figure(figsize=(W, H), dpi=220); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_visible(False)
ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=11, c=BODY, w="normal", ha="left": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
rule = lambda y, c=RULE, lw=0.9, x0=L: ax.plot([x0, R], [y, y], color=c, lw=lw,
                                               solid_capstyle="butt")

T(L, H - 0.62, "Window-free per-residue models: the data split", 20, INK, "bold")
for i, s in enumerate(SUB):
    T(L, H - 1.05 - i * 0.30, s, 11, MUTED)

CX =    [L,      3.35,    4.60,   5.75,   7.00,    8.40,   9.70,   R]
ALIGN = ["left", "right", "right","right","right", "right","right","right"]
HEADS = ["Split", "Rows", "Share", "Precursors", "Groups", "Pocket", "Body", "Any peptide"]
y = H - 1.75 - 0.30 * (len(SUB) - 1)
for cx, h, al in zip(CX, HEADS, ALIGN):
    T(cx, y, h.upper(), 9.5, MUTED, "bold", ha=al)
rule(y - 0.19, INK, 1.1)

y -= 0.62
for i, r in enumerate(rows):
    top = r["indent"] == 0
    T(CX[0] + 0.34 * r["indent"], y, r["name"], 12 if top else 11.5,
      INK if top else BODY, "bold" if top else "normal")
    vals = [f"{r['rows']:,}", f"{r['share']:.1f}%", f"{r['prec']} ({r['labelled']} lab.)",
            f"{r['groups']}", f"{r['pocket']:,}", f"{r['body']:,}", f"{r['any']:,}"]
    for cx, v in zip(CX[1:], vals):
        T(cx, y, v, 11.5 if top else 11, INK if top else BODY,
          "bold" if top else "normal", ha="right")
    y -= ROW
    if i == len(rows) - 2:                 # between the TRAIN block and VALIDATION
        rule(y + (ROW - SEP) / 2); y -= SEP
rule(y + 0.14)
for i, n in enumerate(NOTES):
    T(L, y - 0.14 - i * 0.32, n, 10, MUTED)

out = pathlib.Path.home()/"Desktop/lf_presentation/resid_data_split.png"
fig.savefig(out, facecolor="white")
print("wrote", out, f"({out.stat().st_size//1024} KB)")
