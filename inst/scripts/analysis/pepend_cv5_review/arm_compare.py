"""Compare the 5-seed arms against a matched 5-seed baseline.

Seeds are assigned 1..n by member index, so members 000-004 of the shipped
20-seed C_near_slack8 run ARE a 5-seed baseline with the same seeds -- no
retraining needed, and the comparison is paired.

Everything is recomputed from the member npz files rather than each run's own
printed summary, so every arm goes through one recipe: mean of the RAW member
scores (ranking and both AUCs are rank-based, so the Platt step is irrelevant),
and mean of the member per-index softmaxes.
"""
import csv, json, os, sys
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

A2 = os.path.expanduser("~/AF2_analysis")
S2 = os.path.dirname(os.path.abspath(__file__))
PI = ["CT_cleavage_context", "NT_cleavage_context", "pep_other", "pep_pocket",
      "none", "padding"]
NONE_I, PAD_I = PI.index("none"), PI.index("padding")
SCORED = [i for i in range(len(PI)) if i not in (NONE_I, PAD_I)]

ARMS = [("baseline (focal CE, MHA)", "lf_pepend_run_C_near_slack8"),
        ("CBAM in the U-Net trunk",  "lf_pepend_run_C_near_slack8_cbamtrunk_slack8"),
        ("CBAM at the global head",  "lf_pepend_run_C_near_slack8_cbamhead_slack8"),
        ("Dice per-index loss",      "lf_pepend_run_C_near_slack8_dice_slack8")]

meta = list(csv.DictReader(open(os.path.join(S2, "arm_all_meta.csv"))))
known = np.array([int(r["known"]) for r in meta])
heldout = np.array([int(r["known"]) == 1 and r["split"] == "val" for r in meta])

arr = np.load(os.path.join(A2, "lf_pepend_run_C_near_slack8_isolated/in/arrays.npz"))
val_y = arr["C__val__y_per_index_cat"]                      # (429, 36, 6)
## locate each val window inside `all` by exact bytes, so per-index can be
## scored on held-out windows even though members only save per-index for `all`.
## Materialise both arrays FIRST: indexing an NpzFile re-decompresses the whole
## array every time, and C__all__x is 157 MB.
all_x, val_x = arr["C__all__x"], arr["C__val__x"]
h_all = {all_x[i].tobytes(): i for i in range(all_x.shape[0])}
vrows, vkeep = [], []
for i in range(val_x.shape[0]):
    j = h_all.get(val_x[i].tobytes())
    if j is not None:
        vrows.append(j); vkeep.append(i)
vrows, vkeep = np.array(vrows), np.array(vkeep)
print(f"val windows located inside `all`: {len(vrows)} of {val_x.shape[0]}", flush=True)
del h_all


def overlap(truth1h, pred_soft):
    """Hard-argmax Dice / IoU per class over the scored positions."""
    keep = (truth1h[:, :, PAD_I] == 0) & (truth1h[:, :, NONE_I] == 0) & \
           (truth1h.sum(-1) > 0)
    t = truth1h.argmax(-1)[keep]
    p = pred_soft.argmax(-1)[keep]
    out = {}
    for c in SCORED:
        tc, pc = t == c, p == c
        inter = int((tc & pc).sum())
        if tc.sum() == 0:
            continue
        out[PI[c]] = dict(dice=2 * inter / (tc.sum() + pc.sum()) if (tc.sum() + pc.sum()) else 0.0,
                          iou=inter / ((tc | pc).sum()) if (tc | pc).sum() else 0.0,
                          n=int(tc.sum()))
    acc = float((t == p).mean())
    return acc, out


rows = []
for label, stem in ARMS:
    d = os.path.join(A2, stem + "_isolated/out/members")
    fs = [os.path.join(d, f"member_{k:03d}_C.npz") for k in range(5)]
    if not all(os.path.exists(f) for f in fs):
        print(f"  !! {label}: only {sum(os.path.exists(f) for f in fs)}/5 members, skipping")
        continue
    g = np.mean([np.load(f)["C__global__x"] for f in fs], axis=0)
    v = np.mean([np.load(f)["C__val__x"] for f in fs], axis=0)
    vl = np.load(fs[0])["C__val_labels__x"]
    pidx = np.mean([np.load(f)["C__per_index__x"] for f in fs], axis=0)

    ## rank each held-out end against the UNKNOWN candidates only (10_8's convention)
    unk = np.sort(g[known == 0])[::-1]
    r = np.searchsorted(-unk, -g[heldout], side="left") + 1
    acc, ov = overlap(val_y[vkeep], pidx[vrows])
    md = float(np.mean([ov[c]["dice"] for c in ov]))
    mi = float(np.mean([ov[c]["iou"] for c in ov]))
    rows.append(dict(label=label, roc=roc_auc_score(vl, v), pr=average_precision_score(vl, v),
                     median_rank=float(np.median(r)), top100=int((r <= 100).sum()),
                     top500=int((r <= 500).sum()), acc=acc, dice=md, iou=mi, per_class=ov,
                     params=json.load(open(fs[0].replace(".npz", ".json")))["params"]["C"]))

print(f"\n{'arm':28s}{'params':>8s}{'val ROC':>9s}{'val PR':>8s}"
      f"{'med rank':>10s}{'top100':>8s}{'top500':>8s}{'pi acc':>8s}{'pi Dice':>9s}{'pi IoU':>8s}")
for r in rows:
    print(f"{r['label']:28s}{r['params']:>8,d}{r['roc']:>9.3f}{r['pr']:>8.3f}"
          f"{r['median_rank']:>10.0f}{r['top100']:>8d}{r['top500']:>8d}"
          f"{r['acc']:>8.3f}{r['dice']:>9.3f}{r['iou']:>8.3f}")

print("\nper-residue Dice by class (hard argmax, val windows):")
cls = [c for c in PI if c in rows[0]["per_class"]] if rows else []
print(f"{'arm':28s}" + "".join(f"{c.replace('_cleavage_context','_ctx'):>22s}" for c in cls))
for r in rows:
    print(f"{r['label']:28s}" +
          "".join(f"{r['per_class'][c]['dice']:>22.3f}" for c in cls))
print("positions per class: " +
      ", ".join(f"{c}={rows[0]['per_class'][c]['n']}" for c in cls))
json.dump([{k: v for k, v in r.items()} for r in rows],
          open(os.path.join(A2, "pepend_arm_compare.json"), "w"), indent=1)
