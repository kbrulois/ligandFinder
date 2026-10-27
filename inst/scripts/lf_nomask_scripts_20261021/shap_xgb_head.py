"""Exact TreeSHAP for the per-member xgb window head (Kevin, 2026-10-25).

WHY THIS IS A SEPARATE SCRIPT FROM shap_positions.py. The two heads see
different things. The flatten head reads the (36 x n_channels) window directly,
so attribution there is per residue position in the input -- that is
shap_positions.py. The xgb head does NOT see the window: it sees 48 POOLED
summaries of the trunk's per-index surface, eight statistics for each of the six
per-index classes. So "per residue position" is not defined for its inputs, and
pretending otherwise would be the more misleading answer.

What IS defined, and is what this reports: four of the eight statistics are
positional -- `at_anchor`, `near_anchor` (+/- 2), `pre_anchor` (everything
inside the peptide) and `post_anchor` (everything past the end) -- so TreeSHAP
over them says which REGION of the window the xgb head leans on, per class.
Read together with shap_positions.py's profile, that covers both heads honestly.

Unlike expected gradients this is EXACT, not sampled: `pred_contribs=True` runs
TreeSHAP, whose contributions sum to the margin prediction by construction. The
script asserts that identity rather than trusting it -- it is the one place a
feature/column mix-up would show up.

THE FIT IS THE PRODUCTION FIT. The per-member head is not saved to disk; the
full scan fits it on the fly. So this refits it with full_scan.py's own
defaults (inner_frac 0.3, inner_seed 7, near 2, max_depth 4, eta 0.05,
pm_rounds 100) and the same accession-held-out inner split, so the head
explained here is the head that produced the tracks.

    ARM=curated_wAAid ~/.virtualenvs/r-tensorflow/bin/python shap_xgb_head.py
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")
## lf_dcnn / lf_winxgb from the checkout this script sits in, or $LF_SRC (lf_paths.py)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths  # noqa: E402

# the four positional statistics pool_features emits, and what they cover
POSITIONAL = {"at_anchor": "the end residue itself",
              "near_anchor": "within 2 of the end",
              "pre_anchor": "inside the peptide",
              "post_anchor": "past the end",
              "argmax_pos": "where the class peaks"}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="shap_xgb_head")
    ap.add_argument("--arm", default=os.environ.get("ARM", "curated_wAAid"),
                    help='arm suffix; "" or "curated" for the curated arm')
    ap.add_argument("--terms", default="C,N")
    ap.add_argument("--folds", default="1,2,3,4,5")
    ap.add_argument("--inner-frac", type=float, default=0.3)
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--near", type=int, default=2)
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--pm-rounds", type=int, default=100)
    ap.add_argument("--out-dir", default=f"{A2}/lf_headcmp")
    a = ap.parse_args(argv)

    lf_paths.add_python()
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    import warnings
    import xgboost as xgb
    import tensorflow as tf
    from lf_dcnn.config import PEPEND_ANCHOR, Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import load_member_models
    from lf_winxgb.__main__ import pool_features
    warnings.filterwarnings("ignore", message=r"The structure of `inputs` doesn't match")
    tf.get_logger().setLevel("ERROR")

    arm = "" if a.arm in ("", "curated") else a.arm
    sfx = f"_{arm}" if arm else ""
    tag = f"_nomask_flatten{sfx}"
    rows, rows_chk = [], []

    for term in a.terms.split(","):
        anchor = PEPEND_ANCHOR[term]
        for fold in [int(f) for f in a.folds.split(",")]:
            iso = f"{A2}/lf_pepend_run_{term}_cv5_f{fold}{tag}_isolated"
            grp = f"{A2}/lf_pepend_window_groups_{term}_cv5_f{fold}.npz"
            if not os.path.isdir(iso) or not os.path.exists(grp):
                print(f"  {term} f{fold}: missing run or groups, skipping", flush=True)
                continue
            t0 = time.time()
            cfg = Config.from_json(f"{iso}/in/config.json")
            models = load_member_models(f"{iso}/out", cfg, term)
            pi_names = [str(c) for c in cfg.pi_names]
            z = np.load(f"{iso}/in/arrays.npz", allow_pickle=False)
            Xtr = np.ascontiguousarray(z[f"{term}__train__x"], dtype="float32")
            Xva = np.ascontiguousarray(z[f"{term}__val__x"], dtype="float32")
            ytr_all = z[f"{term}__train__y_global"].reshape(-1).astype("float32")
            yva = z[f"{term}__val__y_global"].reshape(-1).astype("float32")

            # full_scan's own inner split: held out BY ACCESSION, same seed
            gz = np.load(os.path.expanduser(grp), allow_pickle=True)
            acc_tr = np.asarray([str(v) for v in gz[f"{term}__train__accession"]])
            rng = np.random.default_rng(a.inner_seed)
            pos_acc = np.unique(acc_tr[ytr_all == 1])
            oth_acc = np.setdiff1d(np.unique(acc_tr), pos_acc)
            hold = set()
            for pool in (pos_acc, oth_acc):
                pool = pool.copy(); rng.shuffle(pool)
                hold.update(pool[: max(1, int(round(a.inner_frac * pool.size)))].tolist())
            itr = np.flatnonzero(~np.array([v in hold for v in acc_tr]))
            ytr = ytr_all[itr]
            spw = float((ytr == 0).sum()) / max(1.0, float((ytr == 1).sum()))
            base = {"objective": "binary:logistic", "eta": a.eta,
                    "max_depth": a.max_depth, "subsample": 0.7,
                    "colsample_bytree": 0.5, "min_child_weight": 1,
                    "lambda": 1.0, "alpha": 0.0, "tree_method": "hist",
                    "eval_metric": ["auc", "aucpr"]}

            names = None
            phi = None
            n_mem = len(models)
            for mi, m in enumerate(models):
                f = np.asarray(predict_all(m, Xtr)["per_index_cat"], dtype="float32")
                pf, names = pool_features(f, anchor, pi_names, near=a.near)
                b = xgb.train(dict(base, seed=mi, scale_pos_weight=spw),
                              xgb.DMatrix(pf[itr], label=ytr, feature_names=names),
                              num_boost_round=a.pm_rounds, verbose_eval=False)
                fv = np.asarray(predict_all(m, Xva)["per_index_cat"], dtype="float32")
                Fva = pool_features(fv, anchor, pi_names, near=a.near)[0]
                dm = xgb.DMatrix(Fva, feature_names=names)
                c = b.predict(dm, pred_contribs=True)       # (n, n_feat + 1) EXACT
                margin = b.predict(dm, output_margin=True)
                err = float(np.abs(c.sum(axis=1) - margin).max())
                rows_chk.append(dict(arm=arm or "curated", term=term, fold=fold,
                                     member=mi, max_abs_err=err,
                                     n_features=len(names)))
                if phi is None:
                    phi = np.zeros((Fva.shape[0], len(names)))
                phi += c[:, :-1]
            phi /= n_mem

            for st, sel in (("positive", yva == 1), ("negative", yva == 0)):
                if not sel.any():
                    continue
                P = phi[sel]
                for fi, nm in enumerate(names):
                    cls, stat = _split(nm, pi_names)
                    rows.append(dict(arm=arm or "curated", term=term, fold=fold,
                                     set=st, n=int(sel.sum()), feature=nm,
                                     pi_class=cls, statistic=stat,
                                     positional=POSITIONAL.get(stat, ""),
                                     signed=float(P[:, fi].mean()),
                                     abs=float(np.abs(P[:, fi]).mean())))
            worst = max(r["max_abs_err"] for r in rows_chk[-n_mem:])
            print(f"  {term} f{fold}: {n_mem} members, {len(names)} pooled features, "
                  f"{int((yva == 1).sum())} held-out positives; TreeSHAP sum vs margin "
                  f"max |err| {worst:.2e}  ({time.time() - t0:.0f}s)", flush=True)
            del models

    os.makedirs(os.path.expanduser(a.out_dir), exist_ok=True)

    def wr(name, rr):
        if not rr:
            print(f"  (no rows for {name})", flush=True)
            return
        p = os.path.join(os.path.expanduser(a.out_dir), name)
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rr[0].keys()))
            w.writeheader(); w.writerows(rr)
        print(f"-> {p}  ({len(rr):,} rows)", flush=True)

    wr(f"shap_xgb_feature{sfx}.csv", rows)
    wr(f"shap_xgb_exactness{sfx}.csv", rows_chk)
    if rows_chk:
        print(f"TreeSHAP exactness over all fits: max |sum(phi) - margin| = "
              f"{max(r['max_abs_err'] for r in rows_chk):.2e}")


def _split(name: str, pi_names: list[str]) -> tuple[str, str]:
    """pool_features names a feature `<class>_<statistic>`, and both sides can
    contain underscores (`CT_cleavage_context_near_anchor`), so split on the
    class list rather than on the last underscore."""
    for c in sorted(pi_names, key=len, reverse=True):
        if name.startswith(c + "_"):
            return c, name[len(c) + 1:]
    return "", name


if __name__ == "__main__":
    main()
