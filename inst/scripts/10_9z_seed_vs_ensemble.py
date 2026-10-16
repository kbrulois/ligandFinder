"""Is the 20-seed ensemble beaten by its best single seed?

The answer depends entirely on HOW the best seed is chosen, so all three ways
are reported:

  post hoc     the seed with the best retrieval OF THE HELD-OUT ENDS. This is
               selection on the evaluation set: an upper bound nobody can
               achieve in practice, shown only to bound the question.
  on val       the seed with the best window-level PR AUC on the validation
               windows -- the quantity training already early-stops on.
               Recomputed from the weights, never from the saved histories,
               which report a masked metric roughly 14x low.
  on train     the seed with the best retrieval of the TRAIN ends, evaluated
               on the held-out ends. No leakage; this is what a person could
               actually do.
"""
from __future__ import annotations
import csv, json, os, sys
import numpy as np
sys.path.insert(0, "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8/inst/python")

S2 = sys.argv[1]
A = os.path.expanduser("~/AF2_analysis/")


def ends(term):
    prec = list(csv.DictReader(open(f"{S2}/prec_order.csv")))
    idx = {r["accession"]: i for i, r in enumerate(prec)}
    val, train = set(), set()
    for r in csv.DictReader(open(f"{S2}/knowns_near.csv")):
        if r["term"] != term: continue
        i = idx.get(r["accession"])
        if i is None: continue
        (val if r["split"] == "val" else train).add((i, int(r["anchor"])))
    return val, train


def ranks(score, keys, targets, known):
    """rank of each target end among anchors that are not any known end"""
    is_t = np.array([k in targets for k in keys])
    is_k = np.array([k in known for k in keys])
    pos = score[is_t]; neg = np.sort(score[~is_k])
    r = neg.size - np.searchsorted(neg, pos, side="right") + 1
    return r


def summarize(r, n):
    return dict(median=int(np.median(r)), top100=int((r <= 100).sum()),
                top1k=int((r <= 1000).sum()), top10k=int((r <= 10000).sum()), n=int(n))


out = {}
for term, tag, run in (("C", "C_near", "lf_pepend_run_C_near"), ("N", "N_near", "lf_pepend_run_N_near")):
    f = f"{A}lf_pepend_scanmem_{tag}.npz"
    if not os.path.exists(f): print("missing", f); continue
    z = np.load(f)
    M = z["members"]                       # (n_seeds, n_anchors)
    keys = list(zip(z["prot_idx"].tolist(), z["anchor"].tolist()))
    val, train = ends(term); known = val | train
    print(f"\n=== {tag}: {M.shape[0]} seeds, {M.shape[1]:,} anchors, "
          f"{len(val)} held-out ends, {len(train)} train ends", flush=True)

    per_seed = [summarize(ranks(M[s], keys, val, known), len(val)) for s in range(M.shape[0])]
    per_seed_tr = [summarize(ranks(M[s], keys, train, known), len(train)) for s in range(M.shape[0])]
    ens = summarize(ranks(z["score"], keys, val, known), len(val))

    ## selection rule 2: window-level val PR AUC, recomputed from the weights
    from sklearn.metrics import average_precision_score
    from lf_dcnn.config import Config
    from lf_dcnn.pipeline import predict_all
    from lf_dcnn.scan import load_member_models
    iso = f"{A}{run}_isolated"
    cfg = Config.from_json(f"{iso}/in/config.json")
    models = load_member_models(f"{iso}/out", cfg, term)
    arr = np.load(f"{iso}/in/arrays.npz")
    Xv = np.ascontiguousarray(arr[f"{term}__val__x"], dtype="float32")
    yv = arr[f"{term}__val__y_global"].reshape(-1).astype(int)
    val_ap = [float(average_precision_score(yv, np.asarray(predict_all(m, Xv)["global"]).reshape(-1)))
              for m in models]

    med = np.array([p["median"] for p in per_seed])
    best_post = int(np.argmin(med))
    best_val = int(np.argmax(val_ap))
    best_train = int(np.argmin([p["median"] for p in per_seed_tr]))
    out[tag] = dict(term=term, n_seeds=int(M.shape[0]), ensemble=ens, per_seed=per_seed,
                    val_ap=val_ap, best_post=best_post, best_val=best_val, best_train=best_train,
                    seed_median_spread=[int(med.min()), int(np.median(med)), int(med.max())])
    print(f"  ensemble mean of {M.shape[0]}      median rank {ens['median']:7,d}  "
          f"top1k {ens['top1k']:2d}  top10k {ens['top10k']:2d} of {len(val)}")
    print(f"  single seeds, range          median rank {med.min():7,d} - {med.max():,d}  "
          f"(median seed {int(np.median(med)):,d})")
    for nm, s_ in (("best seed, POST HOC", best_post), ("best seed, picked on val windows", best_val),
                   ("best seed, picked on train ends", best_train)):
        p = per_seed[s_]
        print(f"  {nm:34s} seed {s_:2d}  median rank {p['median']:7,d}  "
              f"top1k {p['top1k']:2d}  top10k {p['top10k']:2d}")
    better = int((med < ens["median"]).sum())
    print(f"  seeds beating the ensemble on median rank: {better}/{M.shape[0]}")
json.dump(out, open(f"{S2}/seed_vs_ensemble.json", "w"), indent=1)
print(f"\n-> {S2}/seed_vs_ensemble.json")
