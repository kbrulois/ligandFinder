"""Grouped Shapley for the WINDOW-FREE xgb, the same v and method as the window head.

  groups   a PARTITION of the per-residue channels. Each channel contributes 5
           columns (itself and +/-2 of context), and all 5 move together.
  v(S)     retrieval when only the groups in S survive, everything else
           permuted. Value is -mean(log rank) of the held-out positives among
           all held-out rows, so higher is better and the scale is
           multiplicative in rank -- identical in form to the window version.
  exact    all 2^7 = 128 subsets. `padding` is constant once the signal peptide
           is excluded, so lf_resid drops it and there are 7 groups here, not 8.
"""
import json, sys
from itertools import combinations
from math import factorial
import numpy as np, xgboost as xgb
sys.path.insert(0, "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8/inst/python")
from lf_resid import data
from lf_resid.models import XGB_PARAMS

OUT = sys.argv[1]
PARTITION = [
    ("relASA",              lambda c: c == "relASA"),
    ("conservation",        lambda c: c == "cons_rs_n"),
    ("AFM (min+mean)",      lambda c: c.endswith("_afm")),
    ("secondary structure", lambda c: c.startswith("SS_")),
    ("dihedrals",           lambda c: c.startswith(("Phi_", "Psi_"))),
    ("H-bond energies",     lambda c: "->" in c),
    ("AA properties",       lambda c: c.startswith("AA_")),
]
d = data.load("~/AF2_analysis/lf_pepend_residues.npz", "~/AF2_analysis/lf_resid_labels.npz")
rows, lab_prot, ctrl = data.training_rows(d, n_ctrl=100, seed=42)
tr, va = data.split_by_precursor(d, rows, val_frac=0.4, seed=42, groups=True)
Xtr, Xva = d.x[tr], d.x[va]
base = [n.split("_lag")[0].split("_lead")[0] for n in d.feat_names]
cols = {g: np.array([i for i, b in enumerate(base) if f(b)]) for g, f in PARTITION}
assert sum(len(v) for v in cols.values()) == Xva.shape[1], "partition must cover every column"
print({g: len(v) for g, v in cols.items()}, flush=True)

rng = np.random.default_rng(7)
perm = {g: rng.permutation(Xva.shape[0]) for g in cols}       # one draw per group, reused
names = [g for g, _ in PARTITION]; n = len(names)
out = {}
for target in ("pep_pocket", "pep_any"):
    y_tr, y_va = d.y[target][tr], d.y[target][va]
    p = dict(XGB_PARAMS["reg"]); p.update(objective="binary:logistic", eval_metric="aucpr")
    m = xgb.train(p, xgb.DMatrix(Xtr, label=y_tr), num_boost_round=300)
    pos = np.flatnonzero(y_va == 1)
    def v(S):
        X = Xva.copy()
        for g in names:
            if g not in S:
                X[:, cols[g]] = Xva[perm[g]][:, cols[g]]
        s = m.predict(xgb.DMatrix(X))
        r = np.empty(len(s), dtype=np.int64)
        r[np.argsort(-s, kind="stable")] = np.arange(1, len(s) + 1)
        return -float(np.mean(np.log(r[pos])))
    cache = {}
    for k in range(n + 1):
        for S in combinations(names, k):
            cache[frozenset(S)] = v(set(S))
    phi = {}
    for j in names:
        rest = [g for g in names if g != j]
        tot = 0.0
        for k in range(len(rest) + 1):
            w = factorial(k) * factorial(n - k - 1) / factorial(n)
            for S in combinations(rest, k):
                tot += w * (cache[frozenset(S) | {j}] - cache[frozenset(S)])
        phi[j] = tot
    v_all, v_none = cache[frozenset(names)], cache[frozenset()]
    assert abs(sum(phi.values()) - (v_all - v_none)) < 1e-9, "Shapley must sum to v(all)-v(none)"
    out[target] = dict(shapley=phi, median_none=float(np.exp(-v_none)),
                       median_all=float(np.exp(-v_all)), n_val_pos=int(len(pos)),
                       n_val=int(Xva.shape[0]))
    print(f"{target}: geo-mean rank {np.exp(-v_none):,.0f} -> {np.exp(-v_all):,.0f} "
          f"({np.exp(v_all - v_none):.1f}x), sum check OK", flush=True)
json.dump(out, open(OUT, "w"), indent=1)
print("wrote", OUT)
