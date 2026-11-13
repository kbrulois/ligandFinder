#!/usr/bin/env python3
"""Architecture of the peptide-end CNN trunk, drawn from the built model."""
import json, pathlib
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
S = pathlib.Path(__file__).parent
L = {d["name"]: d for d in json.load(open(S/"cnn_layers.json"))}
P = lambda *n: sum(L[x]["params"] for x in n)

W, H = 14.4, 9.7
fig = plt.figure(figsize=(W, H), dpi=200); fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_visible(False)
ax.set_xlim(0, W); ax.set_ylim(0, H)
T = lambda x, y, s, sz=10, c=BODY, w="normal", ha="center", st="normal": ax.text(
    x, y, s, fontsize=sz, color=c, fontweight=w, ha=ha, va="center", style=st)

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

TOTAL = sum(d["params"] for d in L.values())
T(0.55, H-0.52, "CNN trunk: U-Net over a 36-residue window, with an attention head",
  19, INK, "bold", ha="left")
T(0.55, H-0.92, "Every layer, shape and parameter count is read from the built model, not "
  f"transcribed.  Total {TOTAL:,} parameters.", 10.5, MUTED, ha="left")
T(0.55, H-1.24, "trunk = unet   |   global head = attn   |   position ramp off   |   "
  "l2 1e-3   |   Adam 3e-4, clipnorm 1   |   batch 32", 10, ACCENT, ha="left")

XE, XB, XD = 2.15, 4.30, 6.45
XH, XO = 9.70, 12.85
Y0, Y1, Y2, Y3 = 7.35, 6.30, 4.90, 3.50
SECT = 8.05

T(XB, SECT, "U-NET TRUNK", 10.5, ACCENT, "bold")
ax.plot([0.85, 7.75], [SECT-0.22, SECT-0.22], color=RULE, lw=0.9)
T((XH+XO)/2, SECT, "HEADS", 10.5, ORANGE, "bold")
ax.plot([8.35, W-0.55], [SECT-0.22, SECT-0.22], color=RULE, lw=0.9)

box(XE, Y0, "window  (input)", "(36, 26)", "26 channels", "white")
arrow((XE, Y0-BH/2), (XE, Y1+BH/2))
box(XE, Y1, "enc1  Conv1D 16", "(36, 16)", f"{L['enc1']['params']:,}p", ENC)
box(XE, Y2, "enc2  Conv1D 32", "(18, 32)", f"{L['enc2']['params']:,}p", ENC)
box(XB, Y3, "bottleneck  Conv1D 64", "(9, 64)", f"{L['bottleneck']['params']:,}p", BOT,
    edge=BOT, tc="white")
box(XD, Y2, "dec2  Conv1D 32", "(18, 32)", f"{L['dec2']['params']:,}p", DEC)
box(XD, Y1, "dec1  Conv1D 16", "(36, 16)", f"{L['dec1']['params']:,}p", DEC)

arrow((XE, Y1-BH/2), (XE, Y2+BH/2), "max-pool /2", lx=XE-0.95, ly=(Y1+Y2)/2)
arrow((XE, Y2-BH/2), (XB-BW/2-0.05, Y3+0.10), "max-pool /2", rad=-0.16, lx=XE+0.05, ly=Y3+0.95)
arrow((XB+BW/2+0.05, Y3+0.10), (XD, Y2-BH/2), "upsample x2", rad=-0.16, lx=XD-0.05, ly=Y3+0.95)
arrow((XD, Y2+BH/2), (XD, Y1-BH/2), "upsample x2", lx=XD+0.98, ly=(Y1+Y2)/2)
for yy in (Y1, Y2):                      # label ABOVE the line, not on it
    arrow((XE+BW/2, yy), (XD-BW/2, yy), colour=ACCENT, lw=1.2)
    T((XE+XD)/2, yy+0.20, "skip: concatenate", 8.3, ACCENT, st="italic")


box(XH, Y1, "per_index_logits  Conv1D 6", "(36, 6)", f"{L['per_index_logits']['params']}p", HEAD)
arrow((XD+BW/2, Y1), (XH-BW/2, Y1))
box(XH, Y0, "softmax", "(36, 6)", "", HEAD)
arrow((XH, Y1+BH/2), (XH, Y0-BH/2))
box(XO, Y0, "per-residue class", "6 classes x 36", "", OUT, edge=ORANGE)
arrow((XH+BW/2, Y0), (XO-BW/2, Y0), colour=ORANGE)

chain = [("class-position mask", "(36, 5)", "", 5.30),
         ("softmax over positions", "(36, 5)", "", 4.40),
         ("MultiHeadAttention", "2 heads, key 8", f"{L['attn']['params']}p", 3.50),
         ("LayerNorm + GlobalAvgPool", "(5,)", f"{L['attn_norm']['params']}p", 2.60),
         ("embed  Dense 16, gelu", "dropout 0.3", f"{L['embed']['params']}p", 1.70)]
prev = (XH, Y1-BH/2)
for title, shape, pr, yy in chain:
    box(XH, yy, title, shape, pr, HEAD); arrow(prev, (XH, yy+BH/2)); prev = (XH, yy-BH/2)
box(XH, 0.80, "global  Dense 1, sigmoid", "(1,)", f"{L['global']['params']}p", HEAD)
arrow(prev, (XH, 0.80+BH/2))
box(XO, 0.80, "does a peptide end here", "one score per anchor", "", OUT, edge=ORANGE)
arrow((XH+BW/2, 0.80), (XO-BW/2, 0.80), colour=ORANGE)

NOTE_Y = 2.15
ax.plot([0.95, 7.65], [NOTE_Y+0.72, NOTE_Y+0.72], color=RULE, lw=0.9)
T(0.95, NOTE_Y+0.45, "WHY A U-NET", 10, ACCENT, "bold", ha="left")
T(0.95, NOTE_Y-0.12,
  "The per-residue head has to emit a prediction at every one of the 36 positions, so the\n"
  "encoder's downsampling must be undone. The skip connections carry back the positional\n"
  "detail the bottleneck discards -- without them the 9-wide bottleneck could not place a\n"
  "boundary to the residue.",
  9.6, BODY, ha="left")
T(0.95, NOTE_Y-1.05,
  "The 6 per-residue classes: NT context, CT context, peptide (pocket), peptide (other),\n"
  "padding, none.  The class-position mask drops `none` before the attention head pools.",
  9.6, MUTED, ha="left")

T(XO, Y1+0.15, "The per-residue head draws\nthe class tracks on a protein\npage.\n\n"
  "The global head is the window\nscore the exhaustive scan\nranks anchors by.", 9.3, MUTED)

out = pathlib.Path.home()/"Desktop/lf_presentation/cnn_architecture.png"
fig.savefig(out, facecolor="white")
print("wrote", out, f"({out.stat().st_size//1024} KB)")
