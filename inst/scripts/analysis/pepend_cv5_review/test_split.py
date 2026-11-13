import sys, numpy as np
sys.path.insert(0, "inst/python")
from lf_resid import data as D

# 40 precursors, 10 rows each; groups of 4 precursors share a group id
n_prot, per = 40, 10
prot_idx = np.repeat(np.arange(n_prot), per).astype("int32")
prec_group = (np.arange(n_prot) // 4).astype("int32")
d = D.Residues(x=np.zeros((n_prot*per, 3), "float32"), feat_names=["a","b","c"], y={},
               prot_idx=prot_idx, resno=np.tile(np.arange(1, per+1), n_prot).astype("int32"),
               mature=np.ones(n_prot*per, "int8"), labelled=np.ones(n_prot, "int8"),
               accession=np.array([f"P{i}" for i in range(n_prot)]),
               gene=np.array([f"G{i}" for i in range(n_prot)]), prec_group=prec_group)
rows = np.arange(n_prot*per)

tr, va = D.split_by_precursor(d, rows, val_frac=0.4, seed=1)
gtr, gva = set(prec_group[prot_idx[tr]]), set(prec_group[prot_idx[va]])
ptr, pva = set(prot_idx[tr]), set(prot_idx[va])
print(f"groups: train {len(gtr)} val {len(gva)} overlap {len(gtr & gva)}")
print(f"precursors: train {len(ptr)} val {len(pva)} overlap {len(ptr & pva)}")
print(f"val row fraction {len(va)/len(rows):.3f} (target 0.40)")
assert not (gtr & gva), "a similarity group straddles the split"
assert not (ptr & pva), "a precursor straddles the split"

# swap is the exact complement
tr2, va2 = D.split_by_precursor(d, rows, val_frac=0.5, seed=1)
tr3, va3 = D.split_by_precursor(d, rows, val_frac=0.5, seed=1, swap=True)
assert set(tr2) == set(va3) and set(va2) == set(tr3), "swap is not the complement"
print("swap is the exact complement: ok")

# groups=False reverts to per-precursor and is allowed to straddle groups
tr4, va4 = D.split_by_precursor(d, rows, val_frac=0.4, seed=1, groups=False)
g4 = set(prec_group[prot_idx[tr4]]) & set(prec_group[prot_idx[va4]])
print(f"groups=False: {len(g4)} groups straddle (expected > 0), precursor overlap "
      f"{len(set(prot_idx[tr4]) & set(prot_idx[va4]))}")
assert not (set(prot_idx[tr4]) & set(prot_idx[va4]))

# missing prec_group warns and still splits
d.prec_group = None
tr5, va5 = D.split_by_precursor(d, rows, val_frac=0.4, seed=1)
assert len(va5) > 0 and not (set(prot_idx[tr5]) & set(prot_idx[va5]))
print("fallback without prec_group: ok")
print("\nALL SPLIT TESTS PASSED")
