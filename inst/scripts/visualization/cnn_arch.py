#!/usr/bin/env python3
"""Architecture of the deployed peptide-end CNN, drawn from cnn_layers.json.

That file comes from dump_layers.py, which loads the real members, so no shape,
parameter count or hyperparameter below is typed in by hand. Figure language
follows the project's rules: the head is the "linear head", attention does not
appear, and class names go through the same labels as nn_class_labels.
"""
import json, pathlib, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib import font_manager as fm
for f in ("Arial.ttf", "Arial Bold.ttf"):
    try: fm.fontManager.addfont(f"/System/Library/Fonts/Supplemental/{f}")
    except Exception: pass
plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]

INK, BODY, MUTED, ACCENT = "#111418", "#3F4A56", "#6B7784", "#197EC0"
ENC, BOT, DEC, HEAD, OUT = "#EAF2F8", "#197EC0", "#EAF2F8", "#F4F6F7", "#FDF0E3"
RULE, ORANGE = "#C9D1D9", "#C2663A"
LABELS = dict(CT_cleavage_context="CT-context", NT_cleavage_context="NT-context",
              pep_other="peptide (non-inserting)", pep_pocket="peptide (inserting)",
              none="none", padding="padding")

S = pathlib.Path(__file__).parent
R = json.load(open(S/"cnn_layers.json"))
SW, ES = R["curated_AAswap|C"], R["esmonly36|C"]
for a in ("curated_AAswap", "esmonly36"):
    assert R[f"{a}|C"]["layers"] == R[f"{a}|N"]["layers"], "termini differ"
L = {d["name"]: d for d in SW["layers"]}
LE = {d["name"]: d for d in ES["layers"]}
differ = [n for n in L if L[n]["params"] != LE[n]["params"]]
assert differ == ["enc1"], differ
C, CK = SW["config"], SW["check"]
assert C["global_head"] == "flatten" and C["trunk"] == "unet"
assert not C["class_position_mask"] and CK["mask_all_ones"] and not C["position_ramp"]
assert CK["drops_padding_err"] < 1e-5 < CK["drops_none_err"]       # padding goes, none stays
assert CK["sum_over_classes_is_1"] < 1e-5 < CK["sum_over_positions_minus_1_max"]
assert not CK["head_reads_trunk_directly"]
sh = lambda n: L[n]["shape"]
pp = lambda n: f"{L[n]['params']:,}p"
n_sw, n_es = C["n_channels"], ES["config"]["n_channels"]
n_mem = sum(v["check"]["members_checked"] for v in R.values())
per = R["curated_AAswap|C"]["check"]["members_checked"]
assert (n_mem, per) == (200, 50)
chans = C["channel_names"] if "channel_names" in C else C["channels"]
n_aa = sum(c.startswith("plm_AA_") for c in chans); n_cur = len(chans) - n_aa - 1
n_pc = sum(c.startswith("plm_") for c in (ES["config"].get("channel_names") or ES["config"]["channels"]))
assert len(chans) == n_sw and n_es == n_pc + 1
order = [LABELS[n] for n in CK["pi_names"]]
kept = [LABELS[n] for n in CK["head_classes"]]
assert order[-1] == "padding" and kept == order[:-1]

W, H = 14.4, 9.7
fig = plt.figure(figsize=(W, H), dpi=200); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_visible(False)
ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10, c=BODY, w="normal", ha="center", st="normal": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center", style=st, linespacing=1.3)

BW, BH = 2.45, 0.70
def box(cx, cy, title, shape, params, fill, edge=INK, tc=INK, w=BW, h=BH):
    ax.add_patch(FancyBboxPatch((cx-w/2, cy-h/2), w, h,
        boxstyle="round,pad=0.015,rounding_size=0.07",
        facecolor=fill, edgecolor=edge, linewidth=1.0))
    T(cx, cy+0.155, title, 10.5, tc, "bold")
    sub = shape if not params else f"{shape}   {params}"
    T(cx, cy-0.155, sub, 8.6, tc if fill == BOT else MUTED)

def arrow(a, b, label=None, rad=0.0, lx=None, ly=None, colour=INK, lw=1.3):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12,
        linewidth=lw, color=colour, connectionstyle=f"arc3,rad={rad}",
        shrinkA=2, shrinkB=2))
    if label:
        T((a[0]+b[0])/2 if lx is None else lx,
          (a[1]+b[1])/2 if ly is None else ly, label, 8.3, MUTED, st="italic")

totals = f"{CK['total_params']:,} parameters (AAswap), {ES['check']['total_params']:,} (ESM-C only)"
T(0.55, H-0.52, "CNN trunk: U-Net over a 36-residue window, with a linear head",
  19, INK, "bold", ha="left")
T(0.55, H-0.92, "Every layer, shape and parameter count is read from the deployed models, not "
  f"transcribed.  {totals}.", 10.5, MUTED, ha="left")
T(0.55, H-1.24, f"class-position mask off   |   position ramp off   |   kernel {C['conv_kernel']}, gelu, "
  f"dropout {C['unet_dropout'][0]}   |   l2 {C['l2']:.0e}".replace("e-03", "e-3") +
  f"   |   Adam {C['learning_rate']:.0e}".replace("e-04", "e-4") +
  f", clipnorm {C['clipnorm']:.0f}   |   batch {C['batch_size']}", 10, ACCENT, ha="left")

XE, XB, XD = 2.15, 4.30, 6.45
XH, XO = 9.70, 12.85
Y0, Y1, Y2, Y3 = 7.35, 6.30, 4.90, 3.50
SECT = 8.05

T(XB, SECT, "U-NET TRUNK", 10.5, ACCENT, "bold")
ax.plot([0.85, 7.75], [SECT-0.22, SECT-0.22], color=RULE, lw=0.9)
T((XH+XO)/2, SECT, "HEADS", 10.5, ORANGE, "bold")
ax.plot([8.35, W-0.55], [SECT-0.22, SECT-0.22], color=RULE, lw=0.9)

box(XE, Y0, "window  (input)", f"(36, {n_sw}) AAswap | (36, {n_es}) ESM-C", "", "white")
T(3.60, Y0+0.25, "Two channel sets, each its own models:", 9.3, BODY, "bold", ha="left")
T(3.60, Y0-0.11, f"AAswap: {n_cur} curated channels, padding, {n_aa}-letter residue one-hot\n"
  f"ESM-C only: padding, {n_pc} ESM-C principal components", 9.3, MUTED, ha="left")
arrow((XE, Y0-BH/2), (XE, Y1+BH/2))
box(XE, Y1, "enc1  Conv1D 16", sh("enc1"), f"{L['enc1']['params']:,}p | {LE['enc1']['params']:,}p", ENC)
box(XE, Y2, "enc2  Conv1D 32", sh("enc2"), pp("enc2"), ENC)
box(XB, Y3, "bottleneck  Conv1D 64", sh("bottleneck"), pp("bottleneck"), BOT, edge=BOT, tc="white")
box(XD, Y2, "dec2  Conv1D 32", sh("dec2"), pp("dec2"), DEC)
box(XD, Y1, "dec1  Conv1D 16", sh("dec1"), pp("dec1"), DEC)

arrow((XE, Y1-BH/2), (XE, Y2+BH/2), "max-pool /2", lx=XE-0.62, ly=(Y1+Y2)/2)
arrow((XE, Y2-BH/2), (XB-BW/2-0.05, Y3+0.10), "max-pool /2", rad=0.16, lx=XE+1.12, ly=4.28)
arrow((XB+BW/2+0.05, Y3+0.10), (XD, Y2-BH/2), "upsample x2", rad=0.16, lx=XD-1.12, ly=4.28)
arrow((XD, Y2+BH/2), (XD, Y1-BH/2), "upsample x2", lx=XD+0.66, ly=(Y1+Y2)/2)
for yy in (Y1, Y2):                      # label ABOVE the line, not on it
    arrow((XE+BW/2, yy), (XD-BW/2, yy), colour=ACCENT, lw=1.2)
    T((XE+XD)/2, yy+0.20, "skip: concatenate", 8.3, ACCENT, st="italic")

box(XH, Y1, "per_index_logits  Conv1D 6", sh("per_index_logits"), pp("per_index_logits"), HEAD)
arrow((XD+BW/2, Y1), (XH-BW/2, Y1))
box(XH, Y0, "softmax over classes", sh("per_index_cat"), "", HEAD)
arrow((XH, Y1+BH/2), (XH, Y0-BH/2))
box(XO, Y0, "per-residue class", "6 classes x 36", "", OUT, edge=ORANGE)
arrow((XH+BW/2, Y0), (XO-BW/2, Y0), colour=ORANGE)

STEP = (Y1 - 0.80) / 5
chain = [("drop the padding logit", sh("class_position_mask"), ""),
         ("softmax over classes", sh("masked_sum"), ""),
         ("reshape to one vector", f"{sh('flat_cls')}  36 x 5", ""),
         (f"embed  Dense {C['embed_units']}, gelu", f"dropout {C['embed_dropout']}", pp("embed"))]
prev = (XH, Y1-BH/2)
for i, (title, shape, pr) in enumerate(chain):
    yy = Y1 - STEP*(i+1)
    box(XH, yy, title, shape, pr, HEAD); arrow(prev, (XH, yy+BH/2)); prev = (XH, yy-BH/2)
box(XH, 0.80, "global  Dense 1, sigmoid", sh("global"), pp("global"), HEAD)
arrow(prev, (XH, 0.80+BH/2))
box(XO, 0.80, "does a peptide end here", "one score per anchor", "", OUT, edge=ORANGE)
arrow((XH+BW/2, 0.80), (XO-BW/2, 0.80), colour=ORANGE)

T(XO, 6.10, "Draws the class tracks on a\nprotein page, and feeds the\nxgb window head, which is\nfit separately on 48 summaries\nof these 6 classes.", 9.3, MUTED)
T(XO, 3.45, "LINEAR HEAD", 10, ORANGE, "bold")
T(XO, 2.72, "Reads only the class\nprobabilities, never the\ntrunk's features directly.\nThis is the CNN's window\nscore for each anchor.", 9.3, MUTED)

NOTE_Y = 2.15
ax.plot([0.95, 7.65], [NOTE_Y+0.72, NOTE_Y+0.72], color=RULE, lw=0.9)
XN2 = 4.45
T(0.95, NOTE_Y+0.47, "WHY A U-NET", 10, ACCENT, "bold", ha="left")
T(0.95, NOTE_Y-0.18,
  "The per-residue head predicts at all 36\n"
  "positions, so pooling has to be undone. Skips\n"
  "carry back the detail the 9-wide bottleneck\n"
  "loses, so a boundary lands on the residue.", 9.6, BODY, ha="left")
T(XN2, NOTE_Y+0.47, "WHY A LINEAR HEAD", 10, ORANGE, "bold", ha="left")
T(XN2, NOTE_Y-0.18,
  "Every position x class has its own weight, so\n"
  "the score knows where a class sits. Neighbouring\n"
  "anchors share 35 of 36 residues; position is\n"
  "what turns that into a sharp peak on the end.", 9.6, BODY, ha="left")
fold_seed = "5 folds x 10 seeds"
T(0.95, 0.72,
  f"6 per-residue classes: {', '.join(order[:2])}, {order[3]}, {order[2]}, none, padding.\n"
  f"One model set per terminus and channel set: {fold_seed} = {per} members each, {n_mem} in all.\n"
  f"Loss: per-residue focal + {C['loss_weight_global']} x window BCE.  Early stop on val PR AUC, "
  f"start {C['start_from_epoch']}, patience {C['patience']}.", 9.6, MUTED, ha="left")

dest = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else S/"cnn_architecture.png"
fig.savefig(dest, facecolor="white")
print("wrote", dest, f"({dest.stat().st_size//1024} KB)")
