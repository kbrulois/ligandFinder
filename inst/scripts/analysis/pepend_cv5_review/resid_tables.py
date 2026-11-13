#!/usr/bin/env python3
"""Two presentation tables for the window-free per-residue models."""
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
OUT = pathlib.Path.home()/"Desktop/lf_presentation"
rows = json.load(open(S/"resid_metrics2.json"))
by = {(r["target"], r["model"]): r for r in rows}
NICE = {"pep_any": "Any peptide residue", "pep_pocket": "Pocket residue",
        "pep_other": "Peptide body, not pocket"}
METRICS = (("auc", "ROC AUC", "{:.3f}"), ("ap", "AP", "{:.3f}"),
           ("lift", "AP / prev.", "{:.0f}×"), ("rec1", "Recall @ top 1%", "{:.0%}"))

def render(targets, title, subtitle, notes, out, show_target_col=True):
    nrow = len(targets) * 2
    W = 12.6
    ROW, SEP = 0.42, 0.20
    H = (1.75 + 0.30 * (len(subtitle) - 1)      # title, subtitle lines, header row
         + 0.60                                  # header rule to first data row
         + len(targets) * (2 * ROW + SEP)        # the rows themselves
         + 0.32 * (len(notes) - 1) + 0.52)       # notes and the bottom margin
    fig = plt.figure(figsize=(W, H), dpi=220); fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_visible(False)
    ax.set_xlim(0, W); ax.set_ylim(0, H)
    T = lambda x, y, s, sz=11, c=BODY, w="normal", ha="left": ax.text(
        x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center")
    rule = lambda y, c=RULE, lw=0.9: ax.plot([L, R], [y, y], color=c, lw=lw,
                                             solid_capstyle="butt")
    L, R = 0.6, W - 0.6
    T(L, H - 0.62, title, 20, INK, "bold")
    for i, line in enumerate(subtitle):
        T(L, H - 1.05 - i * 0.30, line, 11, MUTED)

    if show_target_col:
        CX = [L, 4.05, 6.00, 7.25, 8.45, 9.45, 10.55, R]
        HEADS = ["Target", "Model"]
    else:
        CX = [None, L, 4.40, 6.00, 7.80, 9.00, 10.30, R]
        HEADS = [None, "Model"]
    HEADS += [h for _, h, _ in ((None, "Positives", None), (None, "Prevalence", None))]
    HEADS += [h for _, h, _ in METRICS]
    ALIGN = ["left", "left", "right", "right", "right", "right", "right", "right"]
    y = H - 1.75 - 0.30 * (len(subtitle) - 1)
    for cx, h, al in zip(CX, HEADS, ALIGN):
        if cx is not None and h: T(cx, y, h.upper(), 9.5, MUTED, "bold", ha=al)
    rule(y - 0.19, INK, 1.1)

    y -= 0.60
    for t in targets:
        m, x = by[(t, "mlp")], by[(t, "xgb")]
        mid = y - ROW / 2
        if show_target_col: T(CX[0], mid, NICE[t], 12, INK, "bold")
        T(CX[2], mid, f"{m['npos']:,}", 11.5, BODY, ha="right")
        T(CX[3], mid, f"{m['prev']:.2f}%", 11.5, BODY, ha="right")
        for mdl, r in (("MLP", m), ("XGBoost", x)):
            other = x if mdl == "MLP" else m
            T(CX[1], y, mdl, 11.5, BODY)
            for (key, _, fmt), cx in zip(METRICS, CX[4:]):
                mine, theirs = fmt.format(r[key]), fmt.format(other[key])
                best = mine != theirs and r[key] > other[key]
                T(cx, y, mine, 12, ACCENT if best else BODY,
                  "bold" if best else "normal", ha="right")
            y -= ROW
        rule(y + (ROW - SEP) / 2)
        y -= SEP
    for i, n in enumerate(notes):
        T(L, y + 0.02 - i * 0.32, n, 10, MUTED)
    fig.savefig(OUT/out, facecolor="white")
    print("wrote", OUT/out, f"({(OUT/out).stat().st_size//1024} KB)")

SUB = ("30,410 residues from precursors the models never trained on — whole "
       "similarity groups moved to validation, so no homologue straddles the split.")
N_BLUE = ("Blue marks the better of the two models on that row.  AP is "
          "prevalence-dependent, so AP / prevalence states the lift over chance.")
N_LAB  = ("The 45 never-seen peptides carry no labels here — their precursors "
          "were never docked, so these targets are undefined for them.")

render(["pep_any"], "Any peptide residue, held out",
       [SUB, "Pocket and peptide body trained as ONE class."],
       [N_BLUE, N_LAB], "resid_metrics_peptide_residues.png", show_target_col=False)
render(["pep_pocket", "pep_other"], "Pocket vs peptide body, held out",
       [SUB, "The same residues, trained as two separate classes."],
       [N_BLUE, N_LAB], "resid_metrics_by_category.png")
