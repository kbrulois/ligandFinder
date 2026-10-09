"""Replace the CNN's attention-pooled global head with xgboost over its logits.

The window model ends in `attn` -> GlobalAveragePooling -> Dense(16) -> Dense(1,
sigmoid): the 36x6 per-index surface is pooled by learned attention into one
number. This asks whether trees pool it better. Same trunk, same members, same
windows -- only the thing that turns (36, 6) into a score changes, so the
comparison is clean.

Features are the per-index head's (36, 6) output, flattened to 216 columns:
position x class, with position preserved rather than averaged away. That is
the point -- attention pooling is permutation-invariant over positions once its
weights are learned, while a tree can say "class pep_pocket at position 27".

On "logits": the head emits a softmax, and this uses it as-is. Trees are
invariant to any monotone transform of a single feature, so log-odds and
probabilities give identical splits -- converting would cost a pass over the
data and change nothing.

    python -m lf_winxgb --run ~/AF2_analysis/lf_pepend_run_C.rds

Writes an npz with the xgboost head's score for every `all` window beside the
CNN head's, so retrieval can be compared downstream.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np


def _iso(run: str) -> Path:
    return Path(os.path.expanduser(run)).with_suffix("").parent / (
        Path(os.path.expanduser(run)).stem + "_isolated")


def pool_features(f, anchor, pi_names, near=2, skip=()):
    """Summarise a (n, 36, K) per-index surface into a few numbers per class.

    216 raw cells against 41 positives is a wider problem than the data can
    support. These keep the thing attention pooling cannot express -- WHERE on
    the window a class sits -- while collapsing the 36 positions to 8 numbers:

      max, mean, sd      how strong the class is anywhere in the window
      argmax_pos         where its peak is, as a fraction of the window
      at_anchor          its value ON the anchor (the peptide's own last residue)
      near_anchor        its mean over anchor +/- `near`
      pre_anchor         its mean over the peptide body, positions < anchor
      post_anchor        its mean over the cleavage context, positions > anchor

    `pre`/`post`/`at` are the insertion geometry written down directly: a pocket
    that reaches the anchor scores on at/near, one that sits away from it scores
    on pre.

    `skip` drops whole classes. `none` is the obvious candidate: it leads the
    gain table, but it is the vocabulary's bookkeeping column -- "this position
    is not any real class" -- so a head leaning on it may be reading how
    confidently the per-index head declined rather than anything about peptides.
    """
    n, L, K = f.shape
    a = int(anchor) - 1                      # cfg anchors are 1-based
    lo, hi = max(0, a - near), min(L, a + near + 1)
    skip = set(skip)
    out, names = [], []
    for k in range(K):
        if pi_names[k] in skip:
            continue
        v = f[:, :, k]
        cols = {
            "max": v.max(axis=1),
            "mean": v.mean(axis=1),
            "sd": v.std(axis=1),
            "argmax_pos": v.argmax(axis=1).astype("float32") / float(L - 1),
            "at_anchor": v[:, a],
            "near_anchor": v[:, lo:hi].mean(axis=1),
            "pre_anchor": v[:, :a].mean(axis=1) if a > 0 else np.zeros(n, "float32"),
            "post_anchor": (v[:, a + 1:].mean(axis=1) if a + 1 < L
                            else np.zeros(n, "float32")),
        }
        for nm, col in cols.items():
            out.append(col.astype("float32")); names.append(f"{pi_names[k]}_{nm}")
    return np.column_stack(out), names


def fit_stacked(ds, dstop, base, spw, *, n_keep=20, min_rounds=100, nrounds=500,
                early_stopping_rounds=100, max_attempts=200, verbose=False):
    """Fit the stacked head `n_keep` times, keeping only fits that boost far enough.

    A single fit was a single draw, and its stopping point swung 5 -> 182 on the
    seed alone: the stop set is 20 correlated copies of the same windows, so the
    early-stopping curve is nearly flat and noise picks the minimum. Drawing
    until `n_keep` fits clear a floor on best_iteration discards the draws that
    landed in the flat part without having learned anything.

    Lives here, and is imported by scan_genes, because two copies of this is
    exactly how the 5-vs-182 disagreement happened in the first place.
    """
    import xgboost as xgb
    kept, iters, rejected, attempt = [], [], [], 0
    while len(kept) < n_keep and attempt < max_attempts:
        b = xgb.train(dict(base, seed=1000 + attempt, scale_pos_weight=spw), ds,
                      num_boost_round=nrounds, evals=[(dstop, "stop")],
                      early_stopping_rounds=early_stopping_rounds,
                      verbose_eval=verbose)
        bi = int(b.best_iteration)
        if bi >= min_rounds:
            kept.append(b); iters.append(bi)
        else:
            rejected.append(bi)
        attempt += 1
    if len(kept) < n_keep:
        raise SystemExit(f"only {len(kept)} of {n_keep} stacked fits reached "
                         f"{min_rounds} rounds in {attempt} attempts")
    print(f"  stacked head: kept {len(kept)} of {attempt} fits "
          f"(best_iter >= {min_rounds}); kept median {np.median(iters):.0f} "
          f"(range {min(iters)}-{max(iters)})"
          + (f"; rejected {len(rejected)} at median {np.median(rejected):.0f}"
             if rejected else ""), flush=True)
    return kept


def predict_stacked(boosters, dm):
    """Mean over the kept stacked fits."""
    return np.mean([b.predict(dm, iteration_range=(0, b.best_iteration + 1))
                    for b in boosters], axis=0)


def load_plm(path, accs, anchors, how="anchor", near=2):
    """Per-residue PLM channels for each window, keyed on its anchor residue.

    `how`:
      anchor       the 32 channels AT the anchor -- the residue whose identity
                   defines the cleavage
      anchor+near  those, plus the mean over anchor +/- `near`

    Note what these are: your own evaluation found the per-residue embedding
    space is dominated by amino-acid identity, and largely redundant with
    relASA and AlphaMissense, which the CNN already has as channels. So a gain
    here is informative, and no gain is the expected outcome rather than a bug.
    """
    import pyarrow.dataset as ds
    want = {(str(a), int(r)) for a, r in zip(accs, anchors)}
    accs_needed = {a for a, _ in want}
    tbl = ds.dataset(os.path.expanduser(path), format="parquet").to_table(
        filter=ds.field("accession").isin(list(accs_needed)))
    cols = [c for c in tbl.column_names if c.startswith("plm_")]
    import pandas as pd
    df = tbl.to_pandas()
    df["index"] = df["index"].astype("int64")
    piv = {(a, i): v for a, i, v in zip(df["accession"], df["index"],
                                        df[cols].to_numpy(dtype="float32"))}
    n = len(accs)
    out = [np.zeros((n, len(cols)), dtype="float32")]
    names = [f"esm_{c[4:]}" for c in cols]
    miss = 0
    for j, (a, r) in enumerate(zip(accs, anchors)):
        v = piv.get((str(a), int(r)))
        if v is None:
            miss += 1
        else:
            out[0][j] = v
    if how == "anchor+near":
        nb = np.zeros((n, len(cols)), dtype="float32")
        for j, (a, r) in enumerate(zip(accs, anchors)):
            vs = [piv[(str(a), int(r) + d)] for d in range(-near, near + 1)
                  if (str(a), int(r) + d) in piv]
            if vs:
                nb[j] = np.mean(vs, axis=0)
        out.append(nb)
        names = names + [f"{x}_near" for x in names]
    if miss:
        print(f"  PLM: {miss} of {n} anchors not in the table (left as zeros)",
              flush=True)
    return np.concatenate(out, axis=1), names


def _metrics(y, p):
    from sklearn.metrics import average_precision_score, roc_auc_score
    y = np.asarray(y).reshape(-1)
    return {"roc_auc": float(roc_auc_score(y, p)),
            "pr_auc": float(average_precision_score(y, p))}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="lf_winxgb")
    ap.add_argument("--run", default="~/AF2_analysis/lf_pepend_run_C.rds")
    ap.add_argument("--term", default="C")
    ap.add_argument("--out", default=None)
    ap.add_argument("--nrounds", type=int, default=500)
    ap.add_argument("--early-stopping-rounds", type=int, default=100)
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--features", choices=("pooled", "raw", "plm"), default="pooled",
                    help="pooled: 8 summaries per class (48 cols); "
                         "raw: every position x class cell (216 cols); "
                         "plm: the PLM channels ONLY, with nothing from the CNN")
    # One xgboost per CNN member is the default. The stacked head scored better
    # on the exhaustive scan, but it is one tree model fitted on 20 correlated
    # views of the same 41 positives -- its stopping point swung 5 -> 182 on the
    # seed alone, and it needed a best_iter floor plus 20 draws to be stable at
    # all. The per-member head keeps one model per member, which is how every
    # other ensemble in this project is built.
    ap.add_argument("--stack", action="store_true",
                    help="ONE xgboost on all members' surfaces stacked as rows, "
                         "instead of the default one-per-member")
    ap.set_defaults(stack=False)
    ap.add_argument("--stack-seeds", type=int, default=20)
    ap.add_argument("--min-rounds", type=int, default=100)
    ap.add_argument("--plm", default=None,
                    help="parquet of per-residue PLM channels to append, "
                         "e.g. ~/AF2_analysis/lf_plm/esm_c_pca32.parquet")
    ap.add_argument("--plm-how", choices=("anchor", "anchor+near"), default="anchor")
    ap.add_argument("--drop-classes", default="",
                    help="comma-separated per-index classes to leave out of the "
                         "pooled features, e.g. none")
    ap.add_argument("--near", type=int, default=2,
                    help="half-width of the near_anchor window, in residues")
    ap.add_argument("--with-global", action="store_true",
                    help="also give xgboost the CNN's own global score as a feature")
    ap.add_argument("--members", type=int, default=0, help="0 = all")
    ap.add_argument("--groups", default="~/AF2_analysis/lf_pepend_window_groups_C.npz",
                    help="window accessions per split, for the grouped inner split")
    ap.add_argument("--inner-frac", type=float, default=0.3,
                    help="fraction of TRAIN precursors held out for early stopping")
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--no-inner", action="store_true",
                    help="early-stop on val instead (leaks val into model selection)")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args(argv)

    import xgboost as xgb
    from lf_dcnn.config import Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import load_member_models

    iso = _iso(a.run)
    cfg = Config.from_json(iso / "in" / "config.json")
    models = load_member_models(iso / "out", cfg, a.term)
    if a.members:
        models = models[: a.members]
    z = np.load(iso / "in" / "arrays.npz")
    X = {s: np.ascontiguousarray(z[f"{a.term}__{s}__x"], dtype="float32")
         for s in ("train", "val", "all")}
    y = {s: z[f"{a.term}__{s}__y_global"].reshape(-1).astype("float32")
         for s in ("train", "val")}
    pi_names = [str(c) for c in cfg.pi_names]
    from lf_dcnn.scan import PEPEND_ANCHOR
    anchor = PEPEND_ANCHOR[a.term]
    print(f"{a.term}: train {X['train'].shape[0]} ({int(y['train'].sum())} pos), "
          f"val {X['val'].shape[0]} ({int(y['val'].sum())} pos), "
          f"all {X['all'].shape[0]} | {len(models)} members", flush=True)

    # ---- optional PLM channels, one row per window ----------------------------
    PLM = {}
    plm_names = []
    if a.plm:
        gz0 = np.load(os.path.expanduser(a.groups), allow_pickle=True)
        for s_ in ("train", "val", "all"):
            acc = [str(v) for v in gz0[f"{a.term}__{s_}__accession"]]
            anc = [int(v) for v in gz0[f"{a.term}__{s_}__anchor"]]
            if len(acc) != X[s_].shape[0]:
                raise SystemExit(f"groups have {len(acc)} {s_} rows, arrays have "
                                 f"{X[s_].shape[0]}")
            PLM[s_], plm_names = load_plm(a.plm, acc, anc, how=a.plm_how)
        print(f"  PLM: {PLM['train'].shape[1]} channels from "
              f"{os.path.basename(os.path.expanduser(a.plm))} ({a.plm_how})", flush=True)

    # ---- inner split, carved out of TRAIN ------------------------------------
    # Early stopping has to score something, and scoring it on `val` makes val a
    # model-selection set: 500 nested models get ranked on 29 positives and the
    # best is kept. That is why the first run of this looked better than the CNN
    # on val and worse on held-out retrieval. So the stopping set comes out of
    # train instead, grouped by PRECURSOR -- several train windows share an
    # accession, and splitting on rows would put near-identical windows of one
    # protein on both sides.
    if a.no_inner:
        inner_tr = np.arange(X["train"].shape[0])
        inner_st = None
        print("  inner split: OFF -- early stopping on val (val is no longer held out)",
              flush=True)
    else:
        gz = np.load(os.path.expanduser(a.groups), allow_pickle=True)
        acc = np.asarray([str(v) for v in gz[f"{a.term}__train__accession"]])
        if acc.size != X["train"].shape[0]:
            raise SystemExit(f"groups have {acc.size} train rows, arrays have "
                             f"{X['train'].shape[0]} -- different nn_input")
        rng = np.random.default_rng(a.inner_seed)
        # stratify the PRECURSORS by whether they carry a positive, so the
        # stopping set cannot come out with no positives to stop on
        pos_acc = np.unique(acc[y["train"] == 1])
        oth_acc = np.setdiff1d(np.unique(acc), pos_acc)
        hold = set()
        for pool in (pos_acc, oth_acc):
            pool = pool.copy(); rng.shuffle(pool)
            hold.update(pool[: max(1, int(round(a.inner_frac * pool.size)))].tolist())
        is_st = np.array([v in hold for v in acc])
        inner_tr = np.flatnonzero(~is_st)
        inner_st = np.flatnonzero(is_st)
        print(f"  inner split (grouped by precursor, seed {a.inner_seed}): "
              f"fit {inner_tr.size} ({int(y['train'][inner_tr].sum())} pos), "
              f"stop {inner_st.size} ({int(y['train'][inner_st].sum())} pos) "
              f"-- val stays held out", flush=True)
    # 457 training windows against 216 features is a wide, short problem -- hence
    # the shallower trees and slower eta than the residue model uses.
    drop_cls = [c.strip() for c in a.drop_classes.split(",") if c.strip()]
    if drop_cls:
        unknown = set(drop_cls) - set(pi_names)
        if unknown:
            raise SystemExit(f"unknown class(es) {sorted(unknown)}; "
                             f"the vocabulary is {pi_names}")
        print(f"  dropping class(es) from the features: {', '.join(drop_cls)}", flush=True)
    feat_names = ([f"p{p:02d}_{c}" for p in range(cfg.seq_len) for c in pi_names
                   if c not in drop_cls] if a.features == "raw" else None)
    print(f"  features: {a.features}"
          + (f" (anchor at window position {anchor})" if a.features == "pooled" else ""),
          flush=True)

    cnn_val, cnn_all, xgb_val, xgb_all = [], [], [], []
    best_iters = []
    stack = {"tr": [], "st": [], "val": [], "all": []} if a.stack else None
    for mi, m in enumerate(models):
        feats, gl = {}, {}
        for s in ("train", "val", "all"):
            o = predict_all(m, X[s])
            f = np.asarray(o["per_index_cat"], dtype="float32")
            if a.features == "plm":
                # no CNN features at all: the head sees only the embedding at the
                # anchor, so this asks what the PLM alone knows about a cleavage
                feats[s] = np.zeros((f.shape[0], 0), dtype="float32")
                pooled_names = []
            elif a.features == "raw":
                feats[s] = f.reshape(f.shape[0], -1)
                pooled_names = feat_names
            else:
                feats[s], pooled_names = pool_features(f, anchor, pi_names,
                                                       near=a.near, skip=drop_cls)
            gl[s] = np.asarray(o["global"], dtype="float32").reshape(-1)
        names = list(pooled_names)
        if PLM:
            for s_ in ("train", "val", "all"):
                feats[s_] = np.column_stack([feats[s_], PLM[s_]])
            names = names + plm_names
        if a.with_global:
            for s in ("train", "val", "all"):
                feats[s] = np.column_stack([feats[s], gl[s]])
            names = names + ["cnn_global"]
        cnn_val.append(gl["val"]); cnn_all.append(gl["all"])

        dtr = xgb.DMatrix(feats["train"][inner_tr], label=y["train"][inner_tr],
                          feature_names=names)
        dstop = (xgb.DMatrix(feats["train"][inner_st], label=y["train"][inner_st],
                             feature_names=names)
                 if inner_st is not None else None)
        dva = xgb.DMatrix(feats["val"], label=y["val"], feature_names=names)
        base_params = {"objective": "binary:logistic", "eta": a.eta,
                       "max_depth": a.max_depth, "subsample": 0.7,
                       "colsample_bytree": 0.5, "min_child_weight": 1,
                       "lambda": 1.0, "alpha": 0.0, "tree_method": "hist",
                       "eval_metric": ["auc", "aucpr"]}
        params = {"objective": "binary:logistic", "eta": a.eta,
                  "max_depth": a.max_depth, "subsample": 0.7, "colsample_bytree": 0.5,
                  "min_child_weight": 1, "lambda": 1.0, "alpha": 0.0,
                  "tree_method": "hist", "eval_metric": ["auc", "aucpr"],
                  "seed": mi,
                  "scale_pos_weight": float((y["train"][inner_tr] == 0).sum())
                                      / max(1.0, float((y["train"][inner_tr] == 1).sum()))}
        watch = [(dtr, "train"), (dstop if dstop is not None else dva, "stop")]
        if a.stack:
            # Hold this member's view; one model is fitted on all of them below.
            stack["tr"].append(feats["train"][inner_tr])
            if inner_st is not None:
                stack["st"].append(feats["train"][inner_st])
            stack["val"].append(feats["val"])
            stack["all"].append(feats["all"])
            print(f"  member {mi:02d}: features held for stacking | "
                  f"cnn {json.dumps({k: round(v,4) for k,v in _metrics(y['val'], gl['val']).items()})}",
                  flush=True)
            continue

        bst = xgb.train(params, dtr, num_boost_round=a.nrounds, evals=watch,
                        early_stopping_rounds=a.early_stopping_rounds,
                        verbose_eval=bool(a.verbose))
        bi = getattr(bst, "best_iteration", a.nrounds - 1)
        best_iters.append(int(bi))
        rng = {"iteration_range": (0, bi + 1)}
        xgb_val.append(bst.predict(dva, **rng))
        xgb_all.append(bst.predict(
            xgb.DMatrix(feats["all"], feature_names=names), **rng))
        print(f"  member {mi:02d}: best_iter {bi:3d} | "
              f"cnn {json.dumps({k: round(v,4) for k,v in _metrics(y['val'], gl['val']).items()})}"
              f" | xgb {json.dumps({k: round(v,4) for k,v in _metrics(y['val'], xgb_val[-1]).items()})}",
              flush=True)

    if a.stack:
        # One model over every member's view of the same windows. This is NOT 20x
        # the labels -- the same 41 positives recur, one row per member -- so it
        # buys averaging over member disagreement, not independent data.
        Xtr = np.vstack(stack["tr"])
        ytr = np.tile(y["train"][inner_tr], len(models))
        dtr = xgb.DMatrix(Xtr, label=ytr, feature_names=names)
        if stack["st"]:
            dstop = xgb.DMatrix(np.vstack(stack["st"]),
                                label=np.tile(y["train"][inner_st], len(models)),
                                feature_names=names)
        else:
            dstop = xgb.DMatrix(np.vstack(stack["val"]),
                                label=np.tile(y["val"], len(models)),
                                feature_names=names)
        spw = float((ytr == 0).sum()) / max(1.0, float((ytr == 1).sum()))
        print(f"\n  stacked fit: {Xtr.shape[0]:,} rows x {Xtr.shape[1]} features "
              f"({int(ytr.sum())} positives, {len(models)} views of "
              f"{inner_tr.size} windows)", flush=True)
        boosters = fit_stacked(dtr, dstop, base_params, spw,
                               n_keep=a.stack_seeds, min_rounds=a.min_rounds,
                               nrounds=a.nrounds,
                               early_stopping_rounds=a.early_stopping_rounds,
                               verbose=bool(a.verbose))
        best_iters.extend(int(b.best_iteration) for b in boosters)
        # the ensemble applied to each member's view, then averaged -- same
        # ensembling as the per-member path, so the comparison is like for like
        for v in stack["val"]:
            xgb_val.append(predict_stacked(boosters, xgb.DMatrix(v, feature_names=names)))
        for v in stack["all"]:
            xgb_all.append(predict_stacked(boosters, xgb.DMatrix(v, feature_names=names)))
        imp = sorted(boosters[0].get_score(importance_type="gain").items(),
                     key=lambda kv: -kv[1])[:12]
        print("  top gain (first kept fit): "
              + ", ".join(f"{k} {v:.0f}" for k, v in imp), flush=True)

    rep = {"cnn_head": _metrics(y["val"], np.mean(cnn_val, axis=0)),
           "xgb_head": _metrics(y["val"], np.mean(xgb_val, axis=0)),
           "n_members": len(models), "median_best_iter": float(np.median(best_iters)),
           "with_global": bool(a.with_global),
           "inner_split": (not a.no_inner), "inner_frac": a.inner_frac,
           "stack_seeds": a.stack_seeds if a.stack else None,
           "min_rounds": a.min_rounds if a.stack else None,
           "features": a.features, "n_features": len(names), "stacked": bool(a.stack),
           "dropped_classes": drop_cls, "plm": a.plm, "plm_how": a.plm_how}
    print("\nENSEMBLE on val (" + str(len(models)) + " members)")
    for k in ("cnn_head", "xgb_head"):
        print(f"  {k:9s} roc {rep[k]['roc_auc']:.4f}  pr {rep[k]['pr_auc']:.4f}")

    out = a.out or os.path.expanduser(
        str(Path(os.path.expanduser(a.run)).with_name(
            Path(os.path.expanduser(a.run)).stem + "_xgbhead.npz")))
    np.savez_compressed(out,
                        xgb_all=np.mean(xgb_all, axis=0).astype("float32"),
                        xgb_all_sd=np.std(xgb_all, axis=0, ddof=1).astype("float32"),
                        cnn_all=np.mean(cnn_all, axis=0).astype("float32"),
                        xgb_val=np.mean(xgb_val, axis=0).astype("float32"),
                        cnn_val=np.mean(cnn_val, axis=0).astype("float32"),
                        y_val=y["val"], report=np.asarray(json.dumps(rep)))
    print(f"\n-> {out}", flush=True)


if __name__ == "__main__":
    main()
