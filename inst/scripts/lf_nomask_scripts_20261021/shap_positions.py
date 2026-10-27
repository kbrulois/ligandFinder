"""Per-residue-position attribution for a peptide-end arm (Kevin, 2026-10-25).

WHAT THIS COMPUTES, AND WHY IT IS THE SHAPLEY VALUE AND NOT A SALIENCY MAP.

The window is (36 positions x n_channels). The question is how much each
POSITION -- and each channel at that position -- moves the window score. The
estimator is EXPECTED GRADIENTS:

    phi[i,p,c] = E_{b ~ background, alpha ~ U(0,1)}
                 [ (x[i,p,c] - b[p,c]) * d logit(b + alpha (x[i] - b)) / dx[p,c] ]

which is the Aumann-Shapley value of the cell: integrating the gradient along
the straight path from a baseline to x gives the unique attribution that is
exact (sum phi = f(x) - f(baseline)), and averaging over baselines drawn from
the data makes it the Shapley value with respect to the background
distribution rather than to one arbitrary reference. This is the estimator the
`shap` package calls GradientExplainer; it is implemented here directly because
`shap` is in neither venv and NOTHING is ever installed into r-tensorflow
(matplotlib there once upgraded numpy and broke a running render).

THE BASELINE MATTERS. A single all-zero baseline would be wrong in this project:
zero is not "absent", it is a value that already means something (an all-zero
channel block encodes intracellular/TM, and the `padding` channel flags padded
rows). So the background is a random sample of the arm's own candidate windows
-- "what this model usually sees" -- and the attribution reads as "relative to
an average candidate window".

TWO PROPERTIES ARE CHECKED RATHER THAN ASSUMED:
  * COMPLETENESS. sum(phi) must equal logit(x) - mean_b logit(b). With finite
    samples it is an estimate, so the script reports the correlation and the
    relative error per fold. A broken gradient path shows up here immediately.
  * THE LOGIT, NOT THE PROBABILITY. The head ends in a sigmoid, whose gradient
    vanishes wherever the model is confident -- exactly the windows worth
    explaining. The logit is rebuilt as embed %*% W + b from the final Dense
    layer's own weights, so nothing is recomputed by hand and the sigmoid is
    simply never applied. Signs are unchanged: sigmoid is monotonic.

THE EXPLAIN SET IS HELD OUT. For fold k the `val` split is fold k's own
held-out fifth of the knowns (17 of 86 ends for C, 18 of 92 for N, the rest of
`val` being ~400 deliberately hard negatives). So every positive is explained
by ten members that never trained on it, and the five folds together cover
every known end exactly once.

The per-position axis is reported as an OFFSET from the peptide end, not as a
0..35 column index: PEPEND_ANCHOR is 28 for C and 8 for N (1-based), so offset 0
is the end residue itself, negative is inside the peptide and positive is past
the end.

    ARM=curated_wAAid ~/.virtualenvs/r-tensorflow/bin/python shap_positions.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")
## lf_dcnn / lf_winxgb from the checkout this script sits in, or $LF_SRC (lf_paths.py)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths  # noqa: E402

# how channels are grouped in the per-position x family table. Order is the
# order the families are reported in; `family_of` takes the first that matches.
FAMILIES = (
    ("AA identity (one-hot)", lambda c: c.startswith("plm_AA_")),
    ## the motifind arm's indicators are plm_* but are not an embedding; without
    ## this rule they fall into "PLM (ESM-C)" and the family table says ESM-C
    ## for an arm that has none. No other arm has a plm_is_* channel, so this
    ## changes nothing for curated_wAAid or esmonly36.
    ("motif indicators",      lambda c: c.startswith("plm_is_")),
    ("PLM (ESM-C)",           lambda c: c.startswith("plm_") and not c.startswith("plm_AA_")),
    ("AA properties",         lambda c: c.startswith("AA_")),
    ("conservation",          lambda c: c == "cons_rs_n"),
    ("AFM",                   lambda c: c.endswith("_afm")),
    ("relASA",                lambda c: c == "relASA"),
    ("secondary structure",   lambda c: c.startswith("SS_")),
    ("dihedrals",             lambda c: c.startswith(("Phi_", "Psi_"))),
    ("H-bond energies",       lambda c: "->" in c),
    ("padding",               lambda c: c == "padding"),
)


def family_of(c: str) -> str:
    for name, test in FAMILIES:
        if test(c):
            return name
    return "other"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="shap_positions")
    ap.add_argument("--arm", default=os.environ.get("ARM", "curated_wAAid"),
                    help='arm suffix; "" for the curated arm itself')
    ap.add_argument("--terms", default="C,N")
    ap.add_argument("--folds", default="1,2,3,4,5")
    ap.add_argument("--n-background", type=int, default=64,
                    help="background windows x alphas per explained window")
    ap.add_argument("--n-neg", type=int, default=200,
                    help="held-out negatives explained alongside the positives")
    ap.add_argument("--batch", type=int, default=2048)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--out-dir", default=f"{A2}/lf_headcmp")
    a = ap.parse_args(argv)

    lf_paths.add_python()
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    import warnings
    import tensorflow as tf
    import keras

    ## Two warnings fire once per member per fold and say nothing: the trunk's
    ## single input is NAMED ("window"), so handing the sub-model a bare tensor
    ## is correct but noticed; and the last batch of a split is shorter than the
    ## rest, which retraces a tf.function that is cheap to retrace. Neither is a
    ## symptom of anything, and 100 copies of them hide the completeness lines.
    warnings.filterwarnings("ignore", message=r"The structure of `inputs` doesn't match")
    tf.get_logger().setLevel("ERROR")
    from lf_dcnn.config import PEPEND_ANCHOR, Config
    from lf_dcnn.scan import load_member_models

    arm = a.arm
    sfx = f"_{arm}" if arm else ""
    tag = f"_nomask_flatten{sfx}"
    rows_pos, rows_chan, rows_fam, rows_let, rows_chk = [], [], [], [], []
    rows_poschan = []

    for term in a.terms.split(","):
        a0 = PEPEND_ANCHOR[term] - 1              # 0-based column of the end residue
        for fold in [int(f) for f in a.folds.split(",")]:
            iso = f"{A2}/lf_pepend_run_{term}_cv5_f{fold}{tag}_isolated"
            if not os.path.isdir(iso):
                print(f"  {term} f{fold}: no run at {os.path.basename(iso)}, skipping",
                      flush=True)
                continue
            t0 = time.time()
            cfg = Config.from_json(f"{iso}/in/config.json")
            chans = [str(c) for c in cfg.channel_names]
            models = load_member_models(f"{iso}/out", cfg, term)
            z = np.load(f"{iso}/in/arrays.npz", allow_pickle=False)
            xv = np.ascontiguousarray(z[f"{term}__val__x"], dtype="float32")
            yv = z[f"{term}__val__y_global"].reshape(-1)
            xall = np.ascontiguousarray(z[f"{term}__all__x"], dtype="float32")
            yall = z[f"{term}__all__y_global"].reshape(-1)

            rng = np.random.default_rng(a.seed + 97 * fold + (0 if term == "C" else 1))
            ipos = np.flatnonzero(yv == 1)
            ineg = np.flatnonzero(yv == 0)
            if a.n_neg and ineg.size > a.n_neg:
                ineg = rng.choice(ineg, a.n_neg, replace=False)
            X = np.concatenate([xv[ipos], xv[ineg]], axis=0)
            which = np.array(["positive"] * ipos.size + ["negative"] * ineg.size)
            # background: the candidate pool MINUS the knowns, i.e. what the
            # model mostly scores. Drawn once per fold and shared by all members
            # so the completeness reference is the same number for each.
            bgpool = np.flatnonzero(yall == 0)
            bg = xall[rng.choice(bgpool, min(512, bgpool.size), replace=False)]
            print(f"  {term} f{fold}: {len(models)} members, explain "
                  f"{ipos.size} held-out positives + {ineg.size} negatives, "
                  f"background {bg.shape[0]}", flush=True)

            phi = np.zeros((X.shape[0], cfg.seq_len, len(chans)), dtype="float64")
            fx_acc = np.zeros(X.shape[0]); fb_acc = 0.0

            for mi, model in enumerate(models):
                g = model.get_layer("global")
                pen = keras.Model(model.inputs, g.input)
                Wg = tf.convert_to_tensor(g.kernel)
                bg_ = tf.convert_to_tensor(g.bias)

                @tf.function(reduce_retracing=True)
                def grad_logit(xb):
                    with tf.GradientTape() as tape:
                        tape.watch(xb)
                        out = tf.matmul(pen(xb, training=False), Wg) + bg_
                    return tape.gradient(out, xb)

                @tf.function(reduce_retracing=True)
                def logit(xb):
                    return tf.matmul(pen(xb, training=False), Wg) + bg_

                # phi for this member, accumulated over the sampled baselines
                acc = np.zeros_like(phi)
                for s in range(a.n_background):
                    bidx = rng.integers(0, bg.shape[0], X.shape[0])
                    base = bg[bidx]
                    alpha = rng.random((X.shape[0], 1, 1)).astype("float32")
                    pt = base + alpha * (X - base)
                    gr = np.empty_like(X)
                    for i in range(0, X.shape[0], a.batch):
                        sl = slice(i, i + a.batch)
                        gr[sl] = grad_logit(tf.convert_to_tensor(pt[sl])).numpy()
                    acc += (X - base) * gr
                acc /= a.n_background
                phi += acc

                fx_acc += np.concatenate([logit(tf.convert_to_tensor(X[i:i + a.batch])
                                                ).numpy().reshape(-1)
                                          for i in range(0, X.shape[0], a.batch)])
                fb_acc += float(np.mean(
                    np.concatenate([logit(tf.convert_to_tensor(bg[i:i + a.batch])
                                          ).numpy().reshape(-1)
                                    for i in range(0, bg.shape[0], a.batch)])))
            phi /= len(models)
            fx = fx_acc / len(models)
            fb = fb_acc / len(models)

            # ---- COMPLETENESS: sum(phi) == logit(x) - mean_b logit(b) --------
            tot = phi.sum(axis=(1, 2))
            ref = fx - fb
            den = np.maximum(np.abs(ref).mean(), 1e-9)
            r = float(np.corrcoef(tot, ref)[0, 1]) if tot.size > 2 else float("nan")
            rel = float(np.abs(tot - ref).mean() / den)
            rows_chk.append(dict(arm=arm or "curated", term=term, fold=fold,
                                 n=int(X.shape[0]), members=len(models),
                                 n_background=a.n_background, completeness_r=r,
                                 rel_err=rel, mean_ref=float(ref.mean()),
                                 mean_sum_phi=float(tot.mean())))
            print(f"      completeness r={r:.4f} rel.err={rel:.3f}  "
                  f"({time.time() - t0:.0f}s)", flush=True)

            # ---- aggregate ----------------------------------------------------
            fam = np.array([family_of(c) for c in chans])
            for st in ("positive", "negative"):
                m = which == st
                if not m.any():
                    continue
                P = phi[m]                                  # (n, 36, C)
                base = dict(arm=arm or "curated", term=term, fold=fold, set=st,
                            n=int(m.sum()))
                # per position
                psum = P.sum(axis=2)                        # (n, 36)
                pabs = np.abs(P).sum(axis=2)
                for p in range(cfg.seq_len):
                    rows_pos.append(dict(base, offset=p - a0, position=p,
                                         signed=float(psum[:, p].mean()),
                                         abs_of_sum=float(np.abs(psum[:, p]).mean()),
                                         sum_of_abs=float(pabs[:, p].mean())))
                # per channel
                csum = P.sum(axis=1); cabs = np.abs(P).sum(axis=1)
                for ci, c in enumerate(chans):
                    rows_chan.append(dict(base, channel=c, family=fam[ci],
                                          signed=float(csum[:, ci].mean()),
                                          abs_of_sum=float(np.abs(csum[:, ci]).mean()),
                                          sum_of_abs=float(cabs[:, ci].mean())))
                # per position x CHANNEL. The family table sums related
                # channels, which is the right default for ten families but
                # hides what the four AA properties do individually -- and
                # "what does AA_charge contribute at offset +1" is exactly the
                # question an arm without identity channels is asked.
                for ci, c in enumerate(chans):
                    v = P[:, :, ci]
                    for p in range(cfg.seq_len):
                        rows_poschan.append(dict(base, channel=c, family=fam[ci],
                                                 offset=p - a0,
                                                 signed=float(v[:, p].mean()),
                                                 abs=float(np.abs(v[:, p]).mean())))
                # per position x family
                for fname in dict.fromkeys(fam):
                    sel = fam == fname
                    fs = P[:, :, sel].sum(axis=2)
                    fa = np.abs(P[:, :, sel]).sum(axis=2)
                    for p in range(cfg.seq_len):
                        rows_fam.append(dict(base, family=fname, offset=p - a0,
                                             signed=float(fs[:, p].mean()),
                                             abs_of_sum=float(np.abs(fs[:, p]).mean()),
                                             sum_of_abs=float(fa[:, p].mean())))
                # per position x letter -- only meaningful where the one-hot exists.
                # A one-hot cell is nonzero for exactly one letter per position, so
                # the mean is taken over the windows that HAVE that letter there:
                # "when position p is R, what does it contribute".
                oh = [(ci, c[len("plm_AA_"):]) for ci, c in enumerate(chans)
                      if c.startswith("plm_AA_")]
                for ci, letter in oh:
                    v = P[:, :, ci]
                    present = X[m][:, :, ci] > 0.5
                    for p in range(cfg.seq_len):
                        k = present[:, p]
                        if not k.any():
                            continue
                        rows_let.append(dict(base, letter=letter, offset=p - a0,
                                             n_windows=int(k.sum()),
                                             signed=float(v[k, p].mean()),
                                             abs=float(np.abs(v[k, p]).mean())))
            del models
            keras.backend.clear_session()

    os.makedirs(os.path.expanduser(a.out_dir), exist_ok=True)
    import csv as _csv

    def wr(name, rows):
        if not rows:
            print(f"  (no rows for {name})", flush=True)
            return
        p = os.path.join(os.path.expanduser(a.out_dir), name)
        with open(p, "w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"-> {p}  ({len(rows):,} rows)", flush=True)

    wr(f"shap_position{sfx}.csv", rows_pos)
    wr(f"shap_channel{sfx}.csv", rows_chan)
    wr(f"shap_pos_family{sfx}.csv", rows_fam)
    wr(f"shap_pos_channel{sfx}.csv", rows_poschan)
    wr(f"shap_pos_letter{sfx}.csv", rows_let)
    wr(f"shap_completeness{sfx}.csv", rows_chk)


if __name__ == "__main__":
    main()
