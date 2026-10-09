"""Score one gene's step-1 anchors with the xgboost window heads.

The CNN track on the gallery figures is the global head's ensemble mean over 20
members at every anchor. This produces the same thing for the two xgboost heads,
by the same route and the same averaging, so the three tracks are comparable:

  cnn        the attention-pooled sigmoid, averaged over members
  xgb_pm     one xgboost per member on that member's pooled surface, averaged
  xgb_st     ONE xgboost fitted on all members' surfaces stacked, then applied
             to each member's surface and averaged

Heads are fitted here on the candidate windows (the same inner split, grouped by
precursor, that the comparison used) and then applied to the scan anchors, which
the heads have never seen.

    python -m lf_winxgb.scan_genes --genes NPY,CXCL14,ANO8
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from .__main__ import _iso, _metrics, fit_stacked, pool_features, predict_stacked


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lf_winxgb.scan_genes")
    ap.add_argument("--genes", required=True)
    ap.add_argument("--run", default="~/AF2_analysis/lf_pepend_run_C.rds")
    ap.add_argument("--term", default="C")
    ap.add_argument("--residues-npz", default="~/AF2_analysis/lf_pepend_residues.npz")
    ap.add_argument("--labels-npz", default="~/AF2_analysis/lf_resid_labels.npz",
                    help="only for its accession/gene arrays, which index the residues npz")
    ap.add_argument("--groups", default="~/AF2_analysis/lf_pepend_window_groups_C.npz")
    ap.add_argument("--out", default="~/AF2_analysis/lf_winxgb_scan_C.npz")
    ap.add_argument("--inner-frac", type=float, default=0.3)
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--nrounds", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--near", type=int, default=2)
    ap.add_argument("--stack-seeds", type=int, default=20,
                    help="how many stacked fits to KEEP")
    ap.add_argument("--skip-stacked", action="store_true",
                    help="fit only the per-member heads and omit xgb_st from the "
                         "output. The stacked head needs `stack_seeds` fits above "
                         "`min_rounds`; on a flat surface it cannot get them and "
                         "raises rather than return a head that never learned. "
                         "Use this when only xgb_pm is being drawn")
    ap.add_argument("--pm-rounds", type=int, default=0,
                    help="train the PER-MEMBER heads this many rounds with no "
                         "early stopping (0 = early-stop on the inner split). "
                         "Use when the stop set holds too few positives to rank "
                         "anything -- see the comment above the fit")
    ap.add_argument("--min-rounds", type=int, default=100,
                    help="a stacked fit whose best_iteration is below this is "
                         "discarded and another seed drawn")
    ap.add_argument("--max-attempts", type=int, default=200)
    a = ap.parse_args(argv)

    import xgboost as xgb
    from lf_dcnn.config import Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import PEPEND_ANCHOR, Residues, build_windows, load_member_models

    iso = _iso(a.run)
    cfg = Config.from_json(iso / "in" / "config.json")
    models = load_member_models(iso / "out", cfg, a.term)
    pi_names = [str(c) for c in cfg.pi_names]
    anchor = PEPEND_ANCHOR[a.term]
    pad_col = cfg.channel_names.index("padding")

    z = np.load(iso / "in" / "arrays.npz")
    Xc = {s: np.ascontiguousarray(z[f"{a.term}__{s}__x"], dtype="float32")
          for s in ("train", "val")}
    yc = {s: z[f"{a.term}__{s}__y_global"].reshape(-1).astype("float32")
          for s in ("train", "val")}

    # the same grouped inner split the comparison used
    gz = np.load(os.path.expanduser(a.groups), allow_pickle=True)
    acc_tr = np.asarray([str(v) for v in gz[f"{a.term}__train__accession"]])
    rng = np.random.default_rng(a.inner_seed)
    pos_acc = np.unique(acc_tr[yc["train"] == 1])
    oth_acc = np.setdiff1d(np.unique(acc_tr), pos_acc)
    hold = set()
    for pool in (pos_acc, oth_acc):
        pool = pool.copy(); rng.shuffle(pool)
        hold.update(pool[: max(1, int(round(a.inner_frac * pool.size)))].tolist())
    is_st = np.array([v in hold for v in acc_tr])
    itr, ist = np.flatnonzero(~is_st), np.flatnonzero(is_st)
    print(f"inner split: fit {itr.size} ({int(yc['train'][itr].sum())} pos), "
          f"stop {ist.size} ({int(yc['train'][ist].sum())} pos)", flush=True)

    # ---- per-member pooled features on the candidate windows ------------------
    names = None
    F = {"train": [], "val": []}
    for m in models:
        for s in ("train", "val"):
            f = np.asarray(predict_all(m, Xc[s])["per_index_cat"], dtype="float32")
            pf, names = pool_features(f, anchor, pi_names, near=a.near)
            F[s].append(pf)
    base = {"objective": "binary:logistic", "eta": a.eta, "max_depth": a.max_depth,
            "subsample": 0.7, "colsample_bytree": 0.5, "min_child_weight": 1,
            "lambda": 1.0, "alpha": 0.0, "tree_method": "hist",
            "eval_metric": ["auc", "aucpr"]}

    def _spw(y):
        return float((y == 0).sum()) / max(1.0, float((y == 1).sum()))

    # The stacked head already refuses a fit that stopped in a handful of rounds
    # (`--min-rounds`): the stop set is small and its curve nearly flat, so the
    # minimum is noise, and a head truncated there has landed in the flat part
    # rather than learned the surface. The per-member path had no such guard,
    # and on a trunk whose positives were filtered it collapsed -- 8 positives
    # in the stop set, median best_iteration 4, ~5 trees per member, every gene
    # squeezed into a 0.32-0.68 band. `--pm-rounds N` trains a fixed N rounds
    # with no early stopping at all and predicts with the whole booster, which
    # is the honest option when the stop set is too small to rank anything.
    boosters_pm = []
    for mi in range(len(models)):
        p = dict(base, seed=mi, scale_pos_weight=_spw(yc["train"][itr]))
        d1 = xgb.DMatrix(F["train"][mi][itr], label=yc["train"][itr], feature_names=names)
        d2 = xgb.DMatrix(F["train"][mi][ist], label=yc["train"][ist], feature_names=names)
        if a.pm_rounds:
            b = xgb.train(p, d1, num_boost_round=a.pm_rounds, verbose_eval=False)
        else:
            b = xgb.train(p, d1, num_boost_round=a.nrounds, evals=[(d2, "stop")],
                          early_stopping_rounds=100, verbose_eval=False)
        boosters_pm.append(b)
    if a.pm_rounds:
        print(f"per-member heads: {len(boosters_pm)}, fixed {a.pm_rounds} rounds "
              f"(no early stopping)", flush=True)
    else:
        print(f"per-member heads: {len(boosters_pm)}, median best_iter "
              f"{np.median([b.best_iteration for b in boosters_pm]):.0f}", flush=True)

    def _pm_predict(mi, dm):
        """Score with the trees the member actually earned.

        With early stopping that is 0..best_iteration; with `--pm-rounds` the
        booster has no best_iteration and the whole thing is used.
        """
        b = boosters_pm[mi]
        rng = None if a.pm_rounds else (0, b.best_iteration + 1)
        return b.predict(dm, iteration_range=rng) if rng else b.predict(dm)

    # ---- the stacked head, now an ensemble of its own -------------------------
    # One fit was a single draw, and its stopping point swung 5 -> 182 on the
    # seed alone: the stop set is 20 correlated copies of the same windows, so
    # the early-stopping curve is nearly flat and noise picks the minimum. Draw
    # seeds until `stack_seeds` of them survive a floor on best_iteration -- a
    # fit that stops in a handful of rounds has not learned the surface, it has
    # landed in the flat part.
    Xs = np.vstack([F["train"][mi][itr] for mi in range(len(models))])
    ys = np.tile(yc["train"][itr], len(models))
    ds = xgb.DMatrix(Xs, label=ys, feature_names=names)
    dstop = xgb.DMatrix(np.vstack([F["train"][mi][ist] for mi in range(len(models))]),
                        label=np.tile(yc["train"][ist], len(models)), feature_names=names)
    boosters_st = None
    if not a.skip_stacked:
        boosters_st = fit_stacked(ds, dstop, base, _spw(ys),
                                  n_keep=a.stack_seeds, min_rounds=a.min_rounds,
                                  nrounds=a.nrounds, max_attempts=a.max_attempts)
        print(f"stacked head: {Xs.shape[0]:,} rows, {len(boosters_st)} kept fits", flush=True)
    else:
        print("stacked head: skipped (--skip-stacked)", flush=True)

    def _pred_st(dm):
        return predict_stacked(boosters_st, dm)

    # sanity: the heads reproduce the comparison's val numbers
    pv_cnn = np.mean([np.asarray(predict_all(m, Xc["val"])["global"]).reshape(-1)
                      for m in models], axis=0)
    pv_pm = np.mean([_pm_predict(mi, xgb.DMatrix(F["val"][mi], feature_names=names))
                     for mi in range(len(models))], axis=0)
    heads = [("cnn", pv_cnn), ("xgb_pm", pv_pm)]
    if boosters_st is not None:
        heads.append(("xgb_st", np.mean(
            [_pred_st(xgb.DMatrix(F["val"][mi], feature_names=names))
             for mi in range(len(models))], axis=0)))
    for nm, p in heads:
        print(f"  val {nm:6s} " + json.dumps({k: round(v, 4)
                                              for k, v in _metrics(yc["val"], p).items()}),
              flush=True)

    # ---- now the scan anchors, gene by gene -----------------------------------
    res = Residues.from_npz(os.path.expanduser(a.residues_npz))
    lz = np.load(os.path.expanduser(a.labels_npz), allow_pickle=True)
    genes_all = np.asarray([str(g) for g in lz["gene"]])
    accs_all = np.asarray([str(g) for g in lz["accession"]])

    out = {}
    for g in [x.strip() for x in a.genes.split(",")]:
        hits = np.flatnonzero(genes_all == g)
        if not hits.size:
            print(f"  {g}: not in the secretome -- skipped", flush=True)
            continue
        i = int(hits[0])
        anchors = np.arange(int(res.n_prot[i]), int(res.c_prot[i]) + 1)
        x = build_windows(res, i, anchors, cfg, anchor, pad_col)
        pm, st, cn = [], [], []
        for mi, m in enumerate(models):
            o = predict_all(m, x)
            pf, _ = pool_features(np.asarray(o["per_index_cat"], dtype="float32"),
                                  anchor, pi_names, near=a.near)
            dm = xgb.DMatrix(pf, feature_names=names)
            cn.append(np.asarray(o["global"], dtype="float32").reshape(-1))
            pm.append(_pm_predict(mi, dm))
            if boosters_st is not None:
                st.append(_pred_st(dm))
        out[f"{g}__anchor"] = anchors.astype("int32")
        out[f"{g}__accession"] = np.asarray(accs_all[i])
        for nm, v in ([("cnn", cn), ("xgb_pm", pm)] +
                      ([("xgb_st", st)] if boosters_st is not None else [])):
            out[f"{g}__{nm}"] = np.mean(v, axis=0).astype("float32")
            out[f"{g}__{nm}__sd"] = np.std(v, axis=0, ddof=1).astype("float32")
        line = (f"  {g}: {anchors.size} anchors | best cnn {out[f'{g}__cnn'].max():.3f} "
                f"| xgb_pm {out[f'{g}__xgb_pm'].max():.3f}")
        if boosters_st is not None:
            line += f" | xgb_st {out[f'{g}__xgb_st'].max():.3f}"
        print(line, flush=True)

    np.savez_compressed(os.path.expanduser(a.out), **out,
                        genes=np.asarray([k.split("__")[0] for k in out
                                          if k.endswith("__anchor")]))
    print(f"\n-> {a.out}", flush=True)


if __name__ == "__main__":
    main()
