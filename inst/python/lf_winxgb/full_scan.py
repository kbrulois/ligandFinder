"""Score every step-1 anchor in the secretome with the stacked xgboost head.

Mirrors lf_dcnn.scan: walk the precursors in chunks, build each chunk's windows
once, and score them. The difference is what happens after the trunk -- instead
of keeping the attention head's sigmoid, the per-index surface is pooled to 48
summaries and passed to the 20 kept stacked boosters.

Both heads come back, from the same forward pass, so the two tracks over
2,980,852 anchors are strictly comparable.

Work per chunk is 20 members x 20 boosters, so the boosters see each member's
surface stacked into one matrix rather than 400 separate predict calls.

    python -m lf_winxgb.full_scan --run ~/AF2_analysis/lf_pepend_run_C.rds
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from .__main__ import _iso, _metrics, fit_stacked, pool_features, predict_stacked


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lf_winxgb.full_scan")
    ap.add_argument("--run", default="~/AF2_analysis/lf_pepend_run_C.rds")
    ap.add_argument("--term", default="C")
    ap.add_argument("--residues-npz", default="~/AF2_analysis/lf_pepend_residues.npz")
    ap.add_argument("--groups", default="~/AF2_analysis/lf_pepend_window_groups_C.npz")
    ap.add_argument("--out", default="~/AF2_analysis/lf_winxgb_fullscan_C.npz")
    ap.add_argument("--chunk-windows", type=int, default=50_000)
    ap.add_argument("--batch-size", type=int, default=4096)
    ap.add_argument("--inner-frac", type=float, default=0.3)
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--stack-seeds", type=int, default=20)
    ap.add_argument("--min-rounds", type=int, default=100)
    ap.add_argument("--nrounds", type=int, default=500)
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--near", type=int, default=2)
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

    # ---- fit the head on the candidate windows, exactly as the comparison did --
    z = np.load(iso / "in" / "arrays.npz")
    Xc = {s: np.ascontiguousarray(z[f"{a.term}__{s}__x"], dtype="float32")
          for s in ("train", "val")}
    yc = {s: z[f"{a.term}__{s}__y_global"].reshape(-1).astype("float32")
          for s in ("train", "val")}
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

    names = None
    Ftr, Fst, Fva = [], [], []
    for m in models:
        f = np.asarray(predict_all(m, Xc["train"])["per_index_cat"], dtype="float32")
        pf, names = pool_features(f, anchor, pi_names, near=a.near)
        Ftr.append(pf[itr]); Fst.append(pf[ist])
        f = np.asarray(predict_all(m, Xc["val"])["per_index_cat"], dtype="float32")
        Fva.append(pool_features(f, anchor, pi_names, near=a.near)[0])
    ys = np.tile(yc["train"][itr], len(models))
    base = {"objective": "binary:logistic", "eta": a.eta, "max_depth": a.max_depth,
            "subsample": 0.7, "colsample_bytree": 0.5, "min_child_weight": 1,
            "lambda": 1.0, "alpha": 0.0, "tree_method": "hist",
            "eval_metric": ["auc", "aucpr"]}
    spw = float((ys == 0).sum()) / max(1.0, float((ys == 1).sum()))
    boosters = fit_stacked(
        xgb.DMatrix(np.vstack(Ftr), label=ys, feature_names=names),
        xgb.DMatrix(np.vstack(Fst), label=np.tile(yc["train"][ist], len(models)),
                    feature_names=names),
        base, spw, n_keep=a.stack_seeds, min_rounds=a.min_rounds, nrounds=a.nrounds)
    pv = np.mean([predict_stacked(boosters, xgb.DMatrix(v, feature_names=names))
                  for v in Fva], axis=0)
    print("  val xgb_st " + json.dumps({k: round(v, 4)
                                        for k, v in _metrics(yc["val"], pv).items()}),
          flush=True)
    del Ftr, Fst, Fva

    # ---- the exhaustive scan ---------------------------------------------------
    res = Residues.from_npz(os.path.expanduser(a.residues_npz))
    counts = res.n_anchors
    total = int(counts.sum())
    starts = np.concatenate(([0], np.cumsum(counts)))[:-1]
    cnn_all = np.zeros(total, dtype="float32")
    xgb_all = np.zeros(total, dtype="float32")
    xgb_sd = np.zeros(total, dtype="float32")
    prot_idx = np.zeros(total, dtype="int32")
    anchor_of = np.zeros(total, dtype="int32")
    for i in range(len(res)):
        s, n = int(starts[i]), int(counts[i])
        prot_idx[s:s + n] = i
        anchor_of[s:s + n] = np.arange(int(res.n_prot[i]), int(res.c_prot[i]) + 1)
    print(f"scan: {len(res)} precursors, {total:,} {a.term} anchors, "
          f"{len(models)} members x {len(boosters)} boosters", flush=True)

    t0 = time.time()
    done = 0
    i = 0
    while i < len(res):
        j, k = i, 0
        while j < len(res) and (k == 0 or k + counts[j] <= a.chunk_windows):
            k += int(counts[j]); j += 1
        lo, hi = int(starts[i]), int(starts[i]) + k
        x = np.empty((k, cfg.seq_len, cfg.n_channels), dtype="float32")
        at = 0
        for pi in range(i, j):
            an = anchor_of[int(starts[pi]):int(starts[pi]) + int(counts[pi])]
            x[at:at + an.size] = build_windows(res, pi, an, cfg, anchor, pad_col)
            at += int(an.size)

        G = np.empty((len(models), k), dtype="float32")
        P = np.empty((len(models), k, len(names)), dtype="float32")
        for mi, m in enumerate(models):
            o = predict_all(m, x)
            G[mi] = np.asarray(o["global"], dtype="float32").reshape(-1)
            P[mi] = pool_features(np.asarray(o["per_index_cat"], dtype="float32"),
                                  anchor, pi_names, near=a.near)[0]
        # one matrix of every member's view, so each booster predicts once
        dm = xgb.DMatrix(P.reshape(-1, len(names)), feature_names=names)
        S = np.mean([b.predict(dm, iteration_range=(0, b.best_iteration + 1))
                     for b in boosters], axis=0).reshape(len(models), k)
        cnn_all[lo:hi] = G.mean(0)
        xgb_all[lo:hi] = S.mean(0)
        xgb_sd[lo:hi] = S.std(0, ddof=1)
        done += k
        el = time.time() - t0
        print(f"  {done:,}/{total:,} ({100*done/total:.1f}%) "
              f"{el/60:.1f} min elapsed, ~{el/done*(total-done)/60:.1f} min left",
              flush=True)
        i = j

    np.savez_compressed(os.path.expanduser(a.out), cnn=cnn_all, xgb=xgb_all,
                        xgb_sd=xgb_sd, prot_idx=prot_idx, anchor=anchor_of,
                        term=np.asarray(a.term))
    print(f"\n-> {a.out}", flush=True)


if __name__ == "__main__":
    main()
