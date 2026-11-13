"""Precision-recall over LOCI, for each trunk, each ensemble statistic.

Scoring every anchor is misleading -- neighbours share 35 of 36 residues, so one
signal becomes a run of near-identical windows. Two ways to reduce a profile to
loci:

  nms     greedy non-maximum suppression: walk the ranked anchors, keep one
          only if nothing within `--gap` is already kept.
  peaks   scipy.signal.find_peaks with a minimum separation and a prominence
          floor. Narrow by construction, but it REQUIRES a prominent local
          maximum, and about half the known ends sit on a plateau and have
          none -- so it cannot return them at any threshold.

A known end scores as the best locus within `--tol`, or, if it has none, as the
window score at its own anchor, so recall is not capped by the locus caller.
Negatives are loci that are not within tolerance of any known end. Loci near a
real end from OUTSIDE this split -- train ends, ligand_list, UniProt
peptide/propeptide -- are reported both as mistakes and held out, since calling
them mistakes punishes the model for being right about an unlabelled peptide.
"""
from __future__ import annotations
import argparse, csv, json, os
import numpy as np
from scipy.signal import find_peaks


def load_ends(S2, term):
    prec = list(csv.DictReader(open(f"{S2}/prec_order.csv")))
    idx = {r["accession"]: i for i, r in enumerate(prec)}
    val, train, other = set(), set(), set()
    for r in csv.DictReader(open(f"{S2}/knowns_near.csv")):
        if r["term"] != term: continue
        i = idx.get(r["accession"])
        if i is None: continue
        (val if r["split"] == "val" else train).add((i, int(r["anchor"])))
    for r in csv.DictReader(open(f"{S2}/ligand_list.csv")):
        i = idx.get(r["accession"])
        if i is not None: other.add((i, int(r["start"]) if term == "N" else int(r["end"])))
    for r in csv.DictReader(open(f"{S2}/uniprot_peptides.csv")):
        if r["type"] not in ("peptide", "propeptide"): continue
        i = idx.get(r["accession"])
        if i is not None: other.add((i, int(r["start"]) if term == "N" else int(r["end"])))
    return val, (other | train) - val


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True); ap.add_argument("--term", required=True)
    ap.add_argument("--label", required=True); ap.add_argument("--score-key", default="score")
    ap.add_argument("--method", choices=("nms", "peaks"), default="nms")
    ap.add_argument("--scratch", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--gap", type=int, default=5); ap.add_argument("--tol", type=int, default=2)
    ap.add_argument("--distance", type=int, default=4); ap.add_argument("--prominence", type=float, default=0.005)
    a = ap.parse_args()

    z = np.load(os.path.expanduser(a.npz))
    s = np.asarray(z[a.score_key], dtype="float64")
    prot, anc = np.asarray(z["prot_idx"]), np.asarray(z["anchor"])
    val, other = load_ends(a.scratch, a.term)

    if a.method == "nms":
        order = np.argsort(-s)
        taken, L_p, L_c, L_h = set(), [], [], []
        for i in order:
            p, c = int(prot[i]), int(anc[i])
            if any((p, c + d) in taken for d in range(-a.gap, a.gap + 1)): continue
            taken.add((p, c)); L_p.append(p); L_c.append(c); L_h.append(s[i])
        L_p, L_c, L_h = np.array(L_p), np.array(L_c), np.array(L_h)
    else:
        bounds = np.flatnonzero(np.diff(prot)) + 1
        P, C, H = [], [], []
        for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, prot.size]):
            x = s[lo:hi]
            if x.size < 3: continue
            pk, _ = find_peaks(x, distance=a.distance, prominence=a.prominence)
            if pk.size: P.append(np.full(pk.size, prot[lo])); C.append(anc[lo + pk]); H.append(x[pk])
        L_p, L_c, L_h = np.concatenate(P), np.concatenate(C), np.concatenate(H)

    anchor_score = {}
    bounds = np.flatnonzero(np.diff(prot)) + 1
    for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, prot.size]):
        anchor_score[int(prot[lo])] = (anc[lo:hi], s[lo:hi])

    from collections import defaultdict
    by_prot = defaultdict(list)
    for i in range(L_h.size): by_prot[int(L_p[i])].append(i)

    def best_locus(p, c):
        best, bh = None, -1.0
        for i in by_prot.get(p, ()):
            if abs(int(L_c[i]) - c) <= a.tol and L_h[i] > bh: best, bh = i, L_h[i]
        return best

    pos, via = [], []
    for p, c in sorted(val):
        i = best_locus(p, c)
        if i is not None: pos.append(float(L_h[i])); via.append(True)
        else:
            ancs, sc = anchor_score.get(p, (np.array([]), np.array([])))
            j = np.flatnonzero(ancs == c)
            pos.append(float(sc[j[0]]) if j.size else 0.0); via.append(False)
    pos, via = np.asarray(pos), np.asarray(via)

    near_val = np.zeros(L_h.size, bool); near_oth = np.zeros(L_h.size, bool)
    for i in range(L_h.size):
        p, c = int(L_p[i]), int(L_c[i])
        if any((p, c + d) in val for d in range(-a.tol, a.tol + 1)): near_val[i] = True
        elif any((p, c + d) in other for d in range(-a.tol, a.tol + 1)): near_oth[i] = True

    def pr(neg):
        sc = np.r_[pos, neg]; lab = np.r_[np.ones(pos.size), np.zeros(neg.size)]
        o = np.argsort(-sc); lab = lab[o]; sc = sc[o]
        tp = np.cumsum(lab); fp = np.cumsum(1 - lab)
        p_ = tp / np.maximum(tp + fp, 1); r_ = tp / pos.size
        return p_, r_, sc, float(np.sum(p_ * np.diff(np.r_[0.0, r_])))
    pA, rA, sA, apA = pr(L_h[~near_val])
    pB, rB, sB, apB = pr(L_h[~near_val & ~near_oth])

    ## recall at a budget of loci, which is comparable across statistics
    o = np.argsort(-L_h)
    budget = {}
    for N in (100, 500, 1000, 5000, 10000, 50000):
        found = set()
        for i in o[:N]:
            p, c = int(L_p[i]), int(L_c[i])
            for d in range(-a.tol, a.tol + 1):
                if (p, c + d) in val: found.add((p, c + d))
        budget[N] = len(found)
    out = dict(label=a.label, term=a.term, score_key=a.score_key, method=a.method,
               n_loci=int(L_h.size), n_val=int(pos.size), pct_via_locus=float(via.mean()),
               n_neg_all=int((~near_val).sum()), n_neg_excl=int((~near_val & ~near_oth).sum()),
               ap_all=apA, ap_excl=apB, budget=budget)
    json.dump(out, open(f"{a.out}.json", "w"), indent=1)
    with open(f"{a.out}.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["label","score_key","method","arm","precision","recall","score"])
        for nm, (p_, r_, s_) in (("all loci counted", (pA, rA, sA)),
                                 ("unlabelled peptides excluded", (pB, rB, sB))):
            step = max(1, p_.size // 1200)
            for i in range(0, p_.size, step):
                w.writerow([a.label, a.score_key, a.method, nm, round(float(p_[i]),5),
                            round(float(r_[i]),5), round(float(s_[i]),5)])
    print(f"  {a.label:22s} {a.score_key:10s} {a.method:5s} loci {out['n_loci']:7,d}  "
          f"{100*out['pct_via_locus']:3.0f}% via locus  AP {apA:.3f} -> {apB:.3f} excl  "
          f"| ends@1k {budget[1000]:2d} @10k {budget[10000]:2d} of {out['n_val']}", flush=True)


if __name__ == "__main__":
    main()
