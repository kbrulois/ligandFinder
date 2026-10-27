"""Feature importance for the xgb window head, at two levels.

(1) The head's OWN features: 48 pooled summaries of the CNN per-index class
    surface. Gain, averaged over the 20 per-member boosters.

(2) The INPUT channels it ultimately depends on -- relASA, conservation,
    min_afm and the rest -- which the head never sees directly. Measured by
    permuting one (channel, position) cell across the validation windows and
    re-running the whole pipeline (CNN trunk -> pooling -> head); importance is
    the drop in validation ROC AUC. All 36 positions of a channel go through as
    one batch, so it is one forward pass per member per channel.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

## inst/python from the checkout this script sits in, or $LF_SRC (lf_paths.py).
## The installed package has no lf_winxgb and no lf_dcnn/scan.py, so it cannot
## stand in here; this used to be a hardcoded throwaway git worktree.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths; lf_paths.add_python()
from lf_winxgb.__main__ import _iso, pool_features


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--term", default="C")
    ap.add_argument("--groups", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pm-rounds", type=int, default=100)
    ap.add_argument("--inner-frac", type=float, default=0.3)
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--near", type=int, default=2)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--skip-perm", action="store_true")
    a = ap.parse_args(argv)

    import xgboost as xgb
    from sklearn.metrics import roc_auc_score
    from lf_dcnn.config import Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import PEPEND_ANCHOR, load_member_models

    iso = _iso(a.run)
    cfg = Config.from_json(iso / "in" / "config.json")
    models = load_member_models(iso / "out", cfg, a.term)
    pi_names = [str(c) for c in cfg.pi_names]
    chan = [str(c) for c in cfg.channel_names]
    anchor = PEPEND_ANCHOR[a.term]

    z = np.load(iso / "in" / "arrays.npz")
    X = {s: np.ascontiguousarray(z[f"{a.term}__{s}__x"], dtype="float32") for s in ("train", "val")}
    y = {s: z[f"{a.term}__{s}__y_global"].reshape(-1).astype("float32") for s in ("train", "val")}

    ## the same grouped inner split the heads are always fitted on
    gz = np.load(os.path.expanduser(a.groups), allow_pickle=True)
    acc_tr = np.asarray([str(v) for v in gz[f"{a.term}__train__accession"]])
    rng = np.random.default_rng(a.inner_seed)
    pos_acc = np.unique(acc_tr[y["train"] == 1])
    hold = set()
    for pool in (pos_acc, np.setdiff1d(np.unique(acc_tr), pos_acc)):
        pool = pool.copy(); rng.shuffle(pool)
        hold.update(pool[: max(1, int(round(a.inner_frac * pool.size)))].tolist())
    itr = np.flatnonzero(~np.array([v in hold for v in acc_tr]))

    names = None
    Ftr, Fva = [], []
    for m in models:
        f = np.asarray(predict_all(m, X["train"])["per_index_cat"], dtype="float32")
        pf, names = pool_features(f, anchor, pi_names, near=a.near)
        Ftr.append(pf[itr])
        Fva.append(pool_features(np.asarray(predict_all(m, X["val"])["per_index_cat"],
                                            dtype="float32"), anchor, pi_names, near=a.near)[0])
    ytr = y["train"][itr]
    base = {"objective": "binary:logistic", "eta": 0.05, "max_depth": 4, "subsample": 0.7,
            "colsample_bytree": 0.5, "min_child_weight": 1, "lambda": 1.0, "alpha": 0.0,
            "tree_method": "hist", "eval_metric": ["auc", "aucpr"]}
    spw = float((ytr == 0).sum()) / max(1.0, float((ytr == 1).sum()))
    boosters = [xgb.train(dict(base, seed=mi, scale_pos_weight=spw),
                          xgb.DMatrix(Ftr[mi], label=ytr, feature_names=names),
                          num_boost_round=a.pm_rounds, verbose_eval=False)
                for mi in range(len(models))]

    def head_pred(F):
        return np.mean([boosters[mi].predict(xgb.DMatrix(F[mi], feature_names=names))
                        for mi in range(len(models))], axis=0)

    base_auc = roc_auc_score(y["val"], head_pred(Fva))
    print(f"  baseline val ROC {base_auc:.4f}", flush=True)

    ## ---- (1) the head's own features ----------------------------------------
    gain = {}
    for b in boosters:
        for k, v in b.get_score(importance_type="gain").items():
            gain.setdefault(k, []).append(v)
    feat = sorted(((k, float(np.mean(v)), len(v)) for k, v in gain.items()),
                  key=lambda t: -t[1])
    out = dict(run=os.path.basename(a.run), term=a.term, baseline_roc=float(base_auc),
               feature_names=names, channels=chan,
               head_gain=[dict(feature=k, gain=g, members_using=n) for k, g, n in feat])

    ## ---- (2) input channel x window position --------------------------------
    if not a.skip_perm:
        L, C = cfg.seq_len, len(chan)
        Xv = X["val"]; n = Xv.shape[0]
        drop = np.zeros((C, L), dtype="float64")
        t0 = time.time()
        for ci in range(C):
            acc = np.zeros(L)
            for r in range(a.repeats):
                rr = np.random.default_rng(1000 * r + ci)
                ## one batch holding all 36 single-position permutations
                big = np.repeat(Xv[None, ...], L, axis=0)          # (L, n, 36, C)
                for p in range(L):
                    big[p, :, p, ci] = Xv[rr.permutation(n), p, ci]
                big = big.reshape(L * n, L, C)
                F = []
                for m in models:
                    f = np.asarray(predict_all(m, big)["per_index_cat"], dtype="float32")
                    F.append(pool_features(f, anchor, pi_names, near=a.near)[0])
                pred = np.mean([boosters[mi].predict(xgb.DMatrix(F[mi], feature_names=names))
                                for mi in range(len(models))], axis=0).reshape(L, n)
                acc += [base_auc - roc_auc_score(y["val"], pred[p]) for p in range(L)]
            drop[ci] = acc / a.repeats
            el = time.time() - t0
            print(f"  {chan[ci]:16s} max drop {drop[ci].max():+.4f} at position "
                  f"{int(drop[ci].argmax()) + 1:2d}   [{ci+1}/{C}, {el/60:.1f} min, "
                  f"~{el/(ci+1)*(C-ci-1)/60:.0f} min left]", flush=True)
        out["perm_drop"] = drop.tolist()
        out["repeats"] = a.repeats
    with open(os.path.expanduser(a.out), "w") as fh:
        json.dump(out, fh)
    print(f"-> {a.out}", flush=True)


if __name__ == "__main__":
    main()
