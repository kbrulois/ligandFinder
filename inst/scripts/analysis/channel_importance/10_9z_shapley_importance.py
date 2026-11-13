"""Grouped Shapley values for the xgb window head, on retrieval.

Permutation importance gives a group zero credit whenever its information
survives in a correlated neighbour -- and here every structural channel is
71-97% linearly predictable from the others, so single-channel numbers are
near-meaningless and even grouped ones under-credit whatever overlaps.

Shapley fixes the credit assignment: a group's value is its marginal
contribution averaged over every subset of the other groups, so redundant
groups SPLIT the shared signal instead of each showing nothing.

  groups      a PARTITION of the 26 channels (Shapley needs disjoint parts)
  v(S)        retrieval when only the groups in S survive; everything else is
              permuted. Value is -mean(log rank) of the held-out ends, so
              higher is better and the scale is multiplicative in rank.
  exhaustive  all 2^8 = 256 subsets, so the values are exact for this v, not
              sampled. They sum to v(all) - v(none), which is checked.

One permutation vector is drawn per group and REUSED in every subset, so
v(S u {j}) - v(S) differs only by group j and not by the draw.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from itertools import combinations
from math import comb, factorial
import numpy as np
## inst/python from the checkout this script sits in, or $LF_SRC (lf_paths.py).
## The installed package has no lf_winxgb and no lf_dcnn/scan.py, so it cannot
## stand in here; this used to be a hardcoded throwaway git worktree.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lf_paths; lf_paths.add_python()
from lf_winxgb.__main__ import _iso, pool_features

PARTITION = [
    ("relASA",              lambda c: c == "relASA"),
    ("conservation",        lambda c: c == "cons_rs_n"),
    ("AFM (min+mean)",      lambda c: c.endswith("_afm")),
    ("secondary structure", lambda c: c.startswith("SS_")),
    ("dihedrals",           lambda c: c.startswith(("Phi_", "Psi_"))),
    ("H-bond energies",     lambda c: "->" in c),
    ("AA properties",       lambda c: c.startswith("AA_")),
    ("padding",             lambda c: c == "padding"),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--term", default="C")
    ap.add_argument("--groups", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--pm-rounds", type=int, default=100)
    ap.add_argument("--inner-frac", type=float, default=0.3)
    ap.add_argument("--inner-seed", type=int, default=7)
    ap.add_argument("--near", type=int, default=2)
    ap.add_argument("--n-neg", type=int, default=12000,
                    help="unknown candidates kept, to make 256 subsets affordable; "
                         "ranks are reported scaled back to the full candidate set")
    ap.add_argument("--perm-seed", type=int, default=11)
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
    Xtr = np.ascontiguousarray(z[f"{a.term}__train__x"], dtype="float32")
    ytr_all = z[f"{a.term}__train__y_global"].reshape(-1).astype("float32")
    Xall = np.ascontiguousarray(z[f"{a.term}__all__x"], dtype="float32")

    gz = np.load(os.path.expanduser(a.groups), allow_pickle=True)
    acc_tr = np.asarray([str(v) for v in gz[f"{a.term}__train__accession"]])
    rng = np.random.default_rng(a.inner_seed)
    pos_acc = np.unique(acc_tr[ytr_all == 1])
    hold = set()
    for pool in (pos_acc, np.setdiff1d(np.unique(acc_tr), pos_acc)):
        pool = pool.copy(); rng.shuffle(pool)
        hold.update(pool[: max(1, int(round(a.inner_frac * pool.size)))].tolist())
    itr = np.flatnonzero(~np.array([v in hold for v in acc_tr]))

    all_peps = np.asarray([str(v) for v in gz[f"{a.term}__all__peps"]])
    all_known = np.asarray(gz[f"{a.term}__all__known"]).astype(int)
    val_peps = {str(v) for v, k in zip(gz[f"{a.term}__val__peps"], gz[f"{a.term}__val__known"]) if int(k) == 1}
    is_val_pos = np.array([p in val_peps for p in all_peps])
    unk_idx = np.flatnonzero(all_known == 0)
    ## a fixed subsample of the unknowns, identical for every subset
    sub = np.random.default_rng(5).choice(unk_idx, min(a.n_neg, unk_idx.size), replace=False)
    keep = np.concatenate([np.flatnonzero(is_val_pos), np.sort(sub)])
    Xk = np.ascontiguousarray(Xall[keep]); pos_mask = is_val_pos[keep]
    scale = unk_idx.size / sub.size
    print(f"  {int(pos_mask.sum())} held-out ends vs {sub.size:,} of {unk_idx.size:,} unknowns "
          f"(ranks scaled x{scale:.2f})", flush=True)

    names = None; Ftr = []
    for m in models:
        f = np.asarray(predict_all(m, Xtr)["per_index_cat"], dtype="float32")
        pf, names = pool_features(f, anchor, pi_names, near=a.near); Ftr.append(pf[itr])
    ytr = ytr_all[itr]
    base = {"objective": "binary:logistic", "eta": 0.05, "max_depth": 4, "subsample": 0.7,
            "colsample_bytree": 0.5, "min_child_weight": 1, "lambda": 1.0, "alpha": 0.0,
            "tree_method": "hist", "eval_metric": ["auc", "aucpr"]}
    spw = float((ytr == 0).sum()) / max(1.0, float((ytr == 1).sum()))
    boosters = [xgb.train(dict(base, seed=mi, scale_pos_weight=spw),
                          xgb.DMatrix(Ftr[mi], label=ytr, feature_names=names),
                          num_boost_round=a.pm_rounds, verbose_eval=False)
                for mi in range(len(models))]

    gidx = []
    for nm, pred in PARTITION:
        cols = [i for i, c in enumerate(chan) if pred(c)]
        if cols: gidx.append((nm, cols))
    assert sorted(i for _, c in gidx for i in c) == list(range(len(chan))), "not a partition"
    G = len(gidx)
    prng = np.random.default_rng(a.perm_seed)
    perms = [prng.permutation(Xk.shape[0]) for _ in range(G)]   # one draw per group, reused

    def value(S):
        """-mean log rank with only the groups in S surviving."""
        X = Xk.copy()
        for gi, (_, cols) in enumerate(gidx):
            if gi in S: continue
            X[:, :, cols] = Xk[perms[gi]][:, :, cols]
        P = []
        for mi, m in enumerate(models):
            f = np.asarray(predict_all(m, X)["per_index_cat"], dtype="float32")
            P.append(boosters[mi].predict(
                xgb.DMatrix(pool_features(f, anchor, pi_names, near=a.near)[0], feature_names=names)))
        p = np.mean(P, axis=0)
        unk = np.sort(p[~pos_mask])[::-1]
        r = (np.searchsorted(-unk, -p[pos_mask], side="left") + 1) * scale
        return -float(np.mean(np.log(r))), float(np.median(r))

    t0 = time.time(); V = {}
    for k in range(G + 1):
        for S in combinations(range(G), k):
            V[S] = value(set(S))
        print(f"  subsets of size {k}: {comb(G, k):3d} done  [{(time.time()-t0)/60:.0f} min]", flush=True)
    phi = {}
    for gi, (nm, _) in enumerate(gidx):
        others = [j for j in range(G) if j != gi]
        tot = 0.0
        for k in range(G):
            for S in combinations(others, k):
                w = factorial(k) * factorial(G - k - 1) / factorial(G)
                tot += w * (V[tuple(sorted(S + (gi,)))][0] - V[tuple(S)][0])
        phi[nm] = tot
    v_none, v_all = V[()][0], V[tuple(range(G))][0]
    print(f"\n  v(none) {v_none:+.3f} (median rank {V[()][1]:,.0f})   "
          f"v(all) {v_all:+.3f} (median rank {V[tuple(range(G))][1]:,.0f})", flush=True)
    print(f"  sum of Shapley values {sum(phi.values()):+.4f} vs v(all)-v(none) {v_all-v_none:+.4f}", flush=True)
    for nm, v in sorted(phi.items(), key=lambda t: -t[1]):
        print(f"    {nm:22s} {v:+.3f}   ({np.exp(v):.2f}x better ranking)", flush=True)
    json.dump(dict(run=os.path.basename(a.run), term=a.term, groups=[n for n, _ in gidx],
                   shapley=phi, v_none=v_none, v_all=v_all,
                   median_none=V[()][1], median_all=V[tuple(range(G))][1],
                   n_pos=int(pos_mask.sum()), scale=scale),
              open(os.path.expanduser(a.out), "w"), indent=1)
    print(f"-> {a.out}", flush=True)


if __name__ == "__main__":
    main()
