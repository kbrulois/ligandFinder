"""Grouped permutation importance for the xgb window head, scored on RETRIEVAL.

Two changes from the per-cell AUC version, both deliberate:

  * GROUPED. A whole channel is permuted across all 36 positions at once, and
    correlated families (the two AFM channels, the SS one-hots, the dihedrals,
    the H-bond energies, the AA properties) are permuted together. Permuting a
    single cell lets the trunk recover the information from a neighbouring
    channel, so per-cell drops are close to a lower bound.
    The permutation shuffles WINDOW IDENTITY for that block -- each window's
    profile for those channels stays intact but is attached to the wrong window
    -- so the marginal distribution and the within-channel positional structure
    survive while the association with the label does not.

  * RETRIEVAL, not AUC. Importance is what happens to the held-out known ends'
    rank among the candidate windows -- the currency the project actually acts
    on -- rather than a rank statistic over a 400-window validation sample that
    is deliberately enriched for hard negatives. Ranks follow 10_8: a known is
    ranked against the UNKNOWN candidate windows only.
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


def groups_for(chan):
    """Single channels, then the correlated families they sit in."""
    g = {}
    for c in ("relASA", "cons_rs_n", "min_afm", "mean_afm", "padding"):
        if c in chan: g[c] = [c]
    fam = {
        "AFM (min+mean)":      [c for c in chan if c.endswith("_afm")],
        "secondary structure": [c for c in chan if c.startswith("SS_")],
        "dihedrals":           [c for c in chan if c.startswith(("Phi_", "Psi_"))],
        "H-bond energies":     [c for c in chan if "->" in c],
        "AA properties":       [c for c in chan if c.startswith("AA_")],
        "structure (relASA+SS+dihedral+Hbond)":
            [c for c in chan if c == "relASA" or c.startswith(("SS_", "Phi_", "Psi_")) or "->" in c],
        "everything but padding": [c for c in chan if c != "padding"],
    }
    g.update({k: v for k, v in fam.items() if len(v) > 1})
    return g


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--term", default="C")
    ap.add_argument("--groups", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--pm-rounds", type=int, default=100)
    ap.add_argument("--inner-frac", type=float, default=0.3)
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--near", type=int, default=2); ap.add_argument("--repeats", type=int, default=3)
    a = ap.parse_args(argv)

    import xgboost as xgb
    from lf_dcnn.config import Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import PEPEND_ANCHOR, load_member_models

    iso = _iso(a.run)
    cfg = Config.from_json(iso / "in" / "config.json")
    models = load_member_models(iso / "out", cfg, a.term)
    pi_names = [str(c) for c in cfg.pi_names]; chan = [str(c) for c in cfg.channel_names]
    anchor = PEPEND_ANCHOR[a.term]
    z = np.load(iso / "in" / "arrays.npz")
    X = {s: np.ascontiguousarray(z[f"{a.term}__{s}__x"], dtype="float32") for s in ("train", "all")}
    ytr_all = z[f"{a.term}__train__y_global"].reshape(-1).astype("float32")

    gz = np.load(os.path.expanduser(a.groups), allow_pickle=True)
    acc_tr = np.asarray([str(v) for v in gz[f"{a.term}__train__accession"]])
    rng = np.random.default_rng(a.inner_seed)
    pos_acc = np.unique(acc_tr[ytr_all == 1])
    hold = set()
    for pool in (pos_acc, np.setdiff1d(np.unique(acc_tr), pos_acc)):
        pool = pool.copy(); rng.shuffle(pool)
        hold.update(pool[: max(1, int(round(a.inner_frac * pool.size)))].tolist())
    itr = np.flatnonzero(~np.array([v in hold for v in acc_tr]))

    ## which rows of `all` are the HELD-OUT known ends
    all_peps = np.asarray([str(v) for v in gz[f"{a.term}__all__peps"]])
    all_known = np.asarray(gz[f"{a.term}__all__known"]).astype(int)
    val_peps = set(str(v) for v, k in zip(gz[f"{a.term}__val__peps"], gz[f"{a.term}__val__known"]) if int(k) == 1)
    is_val_pos = np.array([p in val_peps for p in all_peps])
    is_unknown = all_known == 0
    print(f"  {int(is_val_pos.sum())} held-out ends among {int(is_unknown.sum()):,} unknown candidates", flush=True)

    names = None; Ftr = []
    for m in models:
        f = np.asarray(predict_all(m, X["train"])["per_index_cat"], dtype="float32")
        pf, names = pool_features(f, anchor, pi_names, near=a.near)
        Ftr.append(pf[itr])
    ytr = ytr_all[itr]
    base = {"objective": "binary:logistic", "eta": 0.05, "max_depth": 4, "subsample": 0.7,
            "colsample_bytree": 0.5, "min_child_weight": 1, "lambda": 1.0, "alpha": 0.0,
            "tree_method": "hist", "eval_metric": ["auc", "aucpr"]}
    spw = float((ytr == 0).sum()) / max(1.0, float((ytr == 1).sum()))
    boosters = [xgb.train(dict(base, seed=mi, scale_pos_weight=spw),
                          xgb.DMatrix(Ftr[mi], label=ytr, feature_names=names),
                          num_boost_round=a.pm_rounds, verbose_eval=False)
                for mi in range(len(models))]

    def retrieval(Xall):
        """Score every candidate, then rank the held-out ends against the unknowns."""
        P = []
        for mi, m in enumerate(models):
            f = np.asarray(predict_all(m, Xall)["per_index_cat"], dtype="float32")
            pf = pool_features(f, anchor, pi_names, near=a.near)[0]
            P.append(boosters[mi].predict(xgb.DMatrix(pf, feature_names=names)))
        p = np.mean(P, axis=0)
        unk = np.sort(p[is_unknown])[::-1]
        r = np.searchsorted(-unk, -p[is_val_pos], side="left") + 1
        return dict(median_rank=float(np.median(r)), mean_log_rank=float(np.mean(np.log(r))),
                    top100=int((r <= 100).sum()), top500=int((r <= 500).sum()),
                    top1000=int((r <= 1000).sum()))

    t0 = time.time()
    b = retrieval(X["all"])
    print(f"  baseline  median rank {b['median_rank']:.0f}  top100 {b['top100']}  "
          f"top500 {b['top500']}  top1000 {b['top1000']}  [{time.time()-t0:.0f}s]", flush=True)

    G = groups_for(chan)
    out = dict(run=os.path.basename(a.run), term=a.term, baseline=b, repeats=a.repeats,
               n_val_pos=int(is_val_pos.sum()), n_unknown=int(is_unknown.sum()), groups={})
    n = X["all"].shape[0]
    for gi, (gname, cols) in enumerate(G.items()):
        idx = [chan.index(c) for c in cols]
        reps = []
        for r in range(a.repeats):
            rr = np.random.default_rng(97 * r + gi)
            Xp = X["all"].copy()
            ## one shuffle of window identity for the whole block, so the block
            ## keeps its internal structure and loses only its link to the window
            Xp[:, :, idx] = X["all"][rr.permutation(n)][:, :, idx]
            reps.append(retrieval(Xp))
        agg = {k: float(np.mean([x[k] for x in reps])) for k in reps[0]}
        agg["median_rank_sd"] = float(np.std([x["median_rank"] for x in reps], ddof=1))
        out["groups"][gname] = dict(channels=cols, **agg)
        el = (time.time() - t0) / 60
        print(f"  {gname:40s} median rank {b['median_rank']:5.0f} -> {agg['median_rank']:6.0f} "
              f"(+-{agg['median_rank_sd']:4.0f})   top500 {b['top500']:2d} -> {agg['top500']:4.1f}"
              f"   [{gi+1}/{len(G)}, {el:.0f} min]", flush=True)
    with open(os.path.expanduser(a.out), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"-> {a.out}", flush=True)


if __name__ == "__main__":
    main()
