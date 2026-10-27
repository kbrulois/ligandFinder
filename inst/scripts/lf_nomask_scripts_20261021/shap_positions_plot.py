#!/usr/bin/env python3
"""Per-residue-position attribution figures, from shap_positions.py.

Two figures:

  shap_positions_<arms>.png  the profile along the window -- how much each
      position moves the window score, and which channel family carries it.
      x is the OFFSET from the peptide end (0 = the end residue itself,
      negative = inside the peptide, positive = past the end), so the C and N
      panels are directly readable against the cleavage site.

  shap_letters_<arm>.png     for an arm with the AA one-hot: signed attribution
      per (letter, offset) -- "when position p is R, what does it contribute".
      Diverging, two hues with a neutral grey midpoint, symmetric about zero.

Units are LOGIT contributions, so a value is "how much this position moved the
pre-sigmoid score relative to an average candidate window", and positions sum
(approximately, it is a sampled estimator) to logit(x) - mean logit(background).
The completeness of that sum is printed in the caption rather than assumed.

Colour = arm (the project's two validated categorical slots); text stays in ink.

    python shap_positions_plot.py ~/AF2_analysis/lf_headcmp curated_wAAid curated
"""
import csv, os, pathlib, sys
from collections import defaultdict

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

A2 = pathlib.Path.home() / "AF2_analysis"
D = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else A2 / "lf_headcmp"
ARMS = sys.argv[2:] or ["curated_wAAid", "curated"]
INK, BODY, MUTED, RULE, GRID = "#0b0b0b", "#3F4A56", "#6B7784", "#C9D1D9", "#E6E5E1"
SLOT = ["#2a78d6", "#eb6834", "#4a9c5d", "#9467bd"]
NEG, MID, POS = "#2166ac", "#f2f2f0", "#b2182b"
MINW = int(os.environ.get("MIN_WINDOWS", "8"))   # support floor for a letter x offset cell
LABEL = {"curated": "curated (26 ch)", "curated_wAAid": "curated + AA one-hot (47 ch)",
         "motifind": "curated + motif indicators (29 ch)",
         "esmonly36": "ESM-C only (37 ch)"}


def sfx(arm):
    return "" if arm == "curated" else f"_{arm}"


def load(arm, kind):
    p = D / f"shap_{kind}{sfx(arm)}.csv"
    if not p.exists():
        return []
    return list(csv.DictReader(open(p)))


def mean_by(rows, key, val, where=lambda r: True):
    """mean over folds of `val`, grouped by `key` -- folds are equal-sized by
    construction (one fifth of the knowns each), so an unweighted mean is the
    right one and is what the caption says."""
    acc = defaultdict(list)
    for r in rows:
        if where(r):
            acc[key(r)].append(float(r[val]))
    return {k: float(np.mean(v)) for k, v in acc.items()}


have = [a for a in ARMS if (D / f"shap_position{sfx(a)}.csv").exists()]
if not have:
    sys.exit(f"no shap_position*.csv for {ARMS} in {D}")
print("arms with attribution:", ", ".join(have))

# =========================================================== figure 1: profile
TERMS = ["C", "N"]
fig = plt.figure(figsize=(15.4, 10.2), dpi=150); fig.patch.set_facecolor("white")
fig.text(0.035, 0.963, "Where in the window the score comes from",
         fontsize=17, color=INK, fontweight="bold")
comp = []
for a in have:
    for r in load(a, "completeness"):
        comp.append(f"{LABEL.get(a, a)} r={float(r['completeness_r']):.3f}")
## the caption is two lines, not one: at this figure width a single line runs
## off the right edge and loses the completeness numbers, which are the part
## that says whether the attribution can be believed at all.
cc = [r for a in have for r in load(a, "completeness")]
cstr = "Completeness is not available."
if cc:
    rs = [float(r["completeness_r"]) for r in cc]
    es = [float(r["rel_err"]) for r in cc]
    cstr = (f"Completeness (sum of cells vs logit(x) - mean logit(background)): r = "
            f"{min(rs):.3f}-{max(rs):.3f}, mean relative error {np.mean(es):.2f}, "
            f"over {len(cc)} fold x terminus fits.")
fig.text(0.035, 0.937,
         "Expected gradients on the pre-sigmoid window score, averaged over the 10 ensemble members of each fold and over the 5 folds. "
         "Every positive is explained by members that never trained on it.",
         fontsize=9.6, color=MUTED)
fig.text(0.035, 0.914, cstr, fontsize=9.6, color=MUTED)

GS = fig.add_gridspec(3, 2, left=0.062, right=0.985, top=0.875, bottom=0.075,
                      hspace=0.42, wspace=0.16,
                      height_ratios=[1.0, 0.72, 1.05])
for ti, term in enumerate(TERMS):
    rows = {a: load(a, "position") for a in have}
    # -- signed, positives vs negatives
    ax = fig.add_subplot(GS[0, ti])
    for ai, a in enumerate(have):
        for st, ls, lw in (("positive", "-", 2.0), ("negative", (0, (2, 2)), 1.4)):
            m = mean_by(rows[a], lambda r: int(r["offset"]), "signed",
                        lambda r: r["term"] == term and r["set"] == st)
            if not m: continue
            x = sorted(m)
            ax.plot(x, [m[k] for k in x], ls=ls, lw=lw, color=SLOT[ai],
                    label=f"{LABEL.get(a, a)}" + ("" if st == "positive" else ", negatives"),
                    zorder=3 if st == "positive" else 2)
    ax.axhline(0, color=RULE, lw=1, zorder=1)
    ax.axvline(0, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
    ax.set_title(f"{term} terminus - signed contribution per position",
                 fontsize=11.5, color=INK, loc="left", pad=7)
    ax.set_ylabel("logit contribution", fontsize=9.5, color=BODY)
    ax.legend(frameon=False, fontsize=8.4, labelcolor=BODY, loc="upper left",
              handlelength=1.8)
    # -- total magnitude
    ax2 = fig.add_subplot(GS[1, ti])
    for ai, a in enumerate(have):
        m = mean_by(rows[a], lambda r: int(r["offset"]), "sum_of_abs",
                    lambda r: r["term"] == term and r["set"] == "positive")
        if not m: continue
        x = sorted(m)
        ax2.plot(x, [m[k] for k in x], lw=1.9, color=SLOT[ai])
    ax2.axvline(0, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
    ax2.set_title("total magnitude, sum |contribution| over channels",
                  fontsize=10.2, color=BODY, loc="left", pad=6)
    ax2.set_ylabel("sum |logit|", fontsize=9.5, color=BODY)
    # -- family breakdown, the arm that has the one-hot first
    ax3 = fig.add_subplot(GS[2, ti])
    fa = load(have[0], "pos_family")
    fams = sorted({r["family"] for r in fa},
                  key=lambda f: -sum(float(r["sum_of_abs"]) for r in fa
                                     if r["family"] == f and r["term"] == term
                                     and r["set"] == "positive"))[:5]
    cmap = plt.get_cmap("viridis")
    for fi, fname in enumerate(fams):
        m = mean_by(fa, lambda r: int(r["offset"]), "sum_of_abs",
                    lambda r: (r["term"] == term and r["set"] == "positive"
                               and r["family"] == fname))
        if not m: continue
        x = sorted(m)
        ax3.plot(x, [m[k] for k in x], lw=1.8, color=cmap(0.08 + 0.82 * fi / max(1, len(fams) - 1)),
                 label=fname)
    ax3.axvline(0, color=MUTED, lw=1, ls=(0, (3, 3)), zorder=1)
    ax3.set_title(f"{LABEL.get(have[0], have[0])}: top channel families",
                  fontsize=10.2, color=BODY, loc="left", pad=6)
    ax3.set_xlabel("offset from the peptide end   (0 = the end residue, + = past it)",
                   fontsize=9.5, color=BODY)
    ax3.set_ylabel("sum |logit|", fontsize=9.5, color=BODY)
    ax3.legend(frameon=False, fontsize=8.0, labelcolor=BODY, loc="upper left",
               handlelength=1.8, ncol=1)
    for A in (ax, ax2, ax3):
        A.set_facecolor("white")
        for s in ("top", "right"): A.spines[s].set_visible(False)
        for s in ("left", "bottom"): A.spines[s].set_color(RULE)
        A.tick_params(labelsize=8.6, colors=BODY, length=3)
        A.grid(axis="y", color=GRID, lw=0.8)
        A.set_axisbelow(True)

p1 = D / f"shap_positions_{'_vs_'.join(have)}.png"
fig.savefig(p1, facecolor="white"); plt.close(fig)
print(f"-> {p1}")

# ==================================================== figure 2: letter heatmap
for a in have:
    lt = load(a, "pos_letter")
    if not lt:
        continue
    fig = plt.figure(figsize=(15.0, 7.4), dpi=150); fig.patch.set_facecolor("white")
    fig.text(0.035, 0.955, f"Which amino acid, at which offset - {LABEL.get(a, a)}",
             fontsize=17, color=INK, fontweight="bold")
    fig.text(0.035, 0.918,
             "Signed contribution of the AA one-hot cell, over the held-out positives that actually carry that letter at that offset "
             f"(weighted by how many). Blank = fewer than {MINW} such ends. Red raises the window score, blue lowers it.",
             fontsize=9.6, color=MUTED)
    cm = matplotlib.colors.LinearSegmentedColormap.from_list("div", [NEG, MID, POS])
    GS = fig.add_gridspec(1, 2, left=0.055, right=0.94, top=0.875, bottom=0.11, wspace=0.13)
    ## Pool the five folds BY WINDOW COUNT, and blank any cell carried by fewer
    ## than MINW held-out ends. There are only ~86 C positives over 36 offsets x
    ## 21 letters, so the average cell holds about four and a plain mean lets a
    ## two-window cell shout as loudly as a forty-window one -- which is exactly
    ## what the first version of this figure did: the rare letters (C, H, W) ran
    ## red right across the window and buried the real signal.
    def cell(sub):
        num, den = defaultdict(float), defaultdict(float)
        for r in sub:
            k = (r["letter"], int(r["offset"])); w = float(r["n_windows"])
            num[k] += float(r["signed"]) * w; den[k] += w
        return {k: (num[k] / den[k], den[k]) for k in num}
    vmax = max((v for term in TERMS
                for v, n in cell([r for r in lt if r["term"] == term
                                  and r["set"] == "positive"]).values()
                if n >= MINW), key=abs, default=1.0)
    vmax = abs(vmax)
    for ti, term in enumerate(TERMS):
        sub = [r for r in lt if r["term"] == term and r["set"] == "positive"]
        if not sub: continue
        letters = sorted({r["letter"] for r in sub})
        offs = sorted({int(r["offset"]) for r in sub})
        M = np.full((len(letters), len(offs)), np.nan)
        for (L, o), (v, n) in cell(sub).items():
            if n >= MINW:
                M[letters.index(L), offs.index(o)] = v
        ax = fig.add_subplot(GS[0, ti])
        im = ax.imshow(M, aspect="auto", cmap=cm, vmin=-vmax, vmax=vmax,
                       extent=(offs[0] - 0.5, offs[-1] + 0.5, len(letters) - 0.5, -0.5),
                       interpolation="nearest")
        ax.set_yticks(range(len(letters))); ax.set_yticklabels(letters, fontsize=8.2)
        ax.axvline(0, color=INK, lw=1.2, ls=(0, (3, 2)))
        ax.set_title(f"{term} terminus", fontsize=11.5, color=INK, loc="left", pad=7)
        ax.set_xlabel("offset from the peptide end", fontsize=9.5, color=BODY)
        ax.set_ylabel("residue at that offset", fontsize=9.5, color=BODY)
        ax.tick_params(labelsize=8.6, colors=BODY, length=3)
        for s in ax.spines.values(): s.set_color(RULE)
        if ti == 1:
            cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.045)
            cb.set_label("logit contribution", fontsize=9, color=BODY)
            cb.ax.tick_params(labelsize=8, colors=BODY)
            cb.outline.set_color(RULE)
    p2 = D / f"shap_letters_{a}.png"
    fig.savefig(p2, facecolor="white"); plt.close(fig)
    print(f"-> {p2}")


# ============================ figure 3: the AA property channels per position
## The curated arm has NO identity channel, so the letter heatmap does not
## exist for it. What it does have is the four AA properties, and this is the
## same question asked of them: what does each contribute at each offset. Four
## rows against the one-hot figure's twenty-one is the point -- it is the whole
## difference between the two arms, drawn to the same axes and the same scale.
PROPS = ["AA_hydro", "AA_charge", "AA_mw", "AA_pI"]
for a in have:
    pc = load(a, "pos_channel")
    if not pc:
        continue
    sub0 = [r for r in pc if r["channel"] in PROPS and r["set"] == "positive"]
    if not sub0:
        continue
    fig = plt.figure(figsize=(15.0, 4.3), dpi=150); fig.patch.set_facecolor("white")
    fig.text(0.035, 0.93, f"What the AA properties contribute, by position - {LABEL.get(a, a)}",
             fontsize=17, color=INK, fontweight="bold")
    fig.text(0.035, 0.855,
             "Signed contribution of each AA property channel, averaged over the held-out positives and the five folds. "
             "Red raises the window score, blue lowers it.",
             fontsize=9.6, color=MUTED)
    cm = matplotlib.colors.LinearSegmentedColormap.from_list("div", [NEG, MID, POS])
    GS = fig.add_gridspec(1, 2, left=0.085, right=0.94, top=0.70, bottom=0.21, wspace=0.13)
    vmax = max(abs(float(r["signed"])) for r in sub0)
    for ti, term in enumerate(TERMS):
        sub = [r for r in sub0 if r["term"] == term]
        if not sub: continue
        offs = sorted({int(r["offset"]) for r in sub})
        M = np.full((len(PROPS), len(offs)), np.nan)
        acc = defaultdict(list)
        for r in sub:
            acc[(r["channel"], int(r["offset"]))].append(float(r["signed"]))
        for (c, o), v in acc.items():
            M[PROPS.index(c), offs.index(o)] = np.mean(v)
        ax = fig.add_subplot(GS[0, ti])
        im = ax.imshow(M, aspect="auto", cmap=cm, vmin=-vmax, vmax=vmax,
                       extent=(offs[0] - 0.5, offs[-1] + 0.5, len(PROPS) - 0.5, -0.5),
                       interpolation="nearest")
        ax.set_yticks(range(len(PROPS))); ax.set_yticklabels(PROPS, fontsize=9)
        ax.axvline(0, color=INK, lw=1.2, ls=(0, (3, 2)))
        ax.set_title(f"{term} terminus", fontsize=11.5, color=INK, loc="left", pad=7)
        ax.set_xlabel("offset from the peptide end", fontsize=9.5, color=BODY)
        ax.tick_params(labelsize=8.6, colors=BODY, length=3)
        for sp in ax.spines.values(): sp.set_color(RULE)
        if ti == 1:
            cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.045)
            cb.set_label("logit contribution", fontsize=9, color=BODY)
            cb.ax.tick_params(labelsize=8, colors=BODY); cb.outline.set_color(RULE)
    p3 = D / f"shap_props_{a}.png"
    fig.savefig(p3, facecolor="white"); plt.close(fig)
    print(f"-> {p3}")
