# Per-gene xgb window-head scan files for ONE precursor the all-secretome scan
# skipped (Kevin, 2026-10-28). allsec_scan.sh keys its *_all.npz by gene name,
# and the residue cache holds two accessions under COL18A1 (and OTOA), so those
# names were left out -- and a page drawn without them has no curated xgb
# window-head track, hence no click-to-open windows.
#
# The full-proteome scans (lf_winxgb_fullscan_<T>_cv5_f<k>_nomask_flatten.npz,
# every anchor of every precursor, indexed by prot_idx = row in the residue
# cache) carry the same xgb_pm / xgb_pm_sd as the per-gene scans -- identical on
# SCG2, all 10 fold x terminus pairs -- so this cuts one accession out of them
# into <T>_cv5_f<k> files, and pools the all-50 file with allsec_scan.sh's
# formula. Only the xgb keys: the page computes the CNN series live, and
# nothing reads a scan's cnn keys with PROV_N=0.
#
#   python own_scan_npz.py build COL18A1 P39060 _col18a1
#   python own_scan_npz.py check SCG2 P13521        # rebuild a gene _all HAS; must match exactly
#
# Render with NPZ_SFX=<sfx> (candidate_pages.sh does it for genes in OWN_NPZ).
# ORDER: the residue cache's accession order (prot_idx), one per line.
import numpy as np, os, sys
A2 = os.path.expanduser("~/AF2_analysis"); TAG = "_nomask_flatten"
ORDER = os.environ.get("ORDER", f"{A2}/lf_nomask_scripts_20261021/residues_accession_order.txt")
n, K = 10, 5   # seeds per fold, folds

def prot_idx(acc):
    order = [l.strip() for l in open(ORDER)]
    hits = [i for i, a in enumerate(order) if a == acc]
    if len(hits) != 1: sys.exit(f"{acc}: {len(hits)} rows in {ORDER}")
    return hits[0]

def folds(T, idx):
    out = []
    for k in range(1, K + 1):
        fs = np.load(f"{A2}/lf_winxgb_fullscan_{T}_cv5_f{k}{TAG}.npz", allow_pickle=True)
        w = fs["prot_idx"] == idx
        out.append((fs["anchor"][w], fs["xgb_pm"][w], fs["xgb_pm_sd"][w]))
    an = out[0][0]
    if not len(an): sys.exit(f"prot_idx {idx}: no anchors in the {T} fullscan")
    if not all(np.array_equal(an, o[0]) for o in out): sys.exit(f"prot_idx {idx}: {T} anchors differ across folds")
    if not np.array_equal(an, np.arange(an[0], an[-1] + 1)): sys.exit(f"prot_idx {idx}: {T} anchors not contiguous")
    return out

def pool(F):   # allsec_scan.sh: mean over 50 members, sd from each fold's 10-member mean and sd
    M = np.stack([f[1] for f in F]).astype("f8"); Sd = np.stack([f[2] for f in F]).astype("f8")
    m = M.mean(0)
    var = ((n - 1) * (Sd ** 2).sum(0) + n * ((M - m) ** 2).sum(0)) / (n * K - 1)
    return m.astype("f4"), np.sqrt(var).astype("f4")

mode, gene, acc = sys.argv[1], sys.argv[2], sys.argv[3]
idx = prot_idx(acc)
if mode == "check":
    for T in ("N", "C"):
        F = folds(T, idx)
        for k, (a, x, s) in enumerate(F, 1):
            z = np.load(f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}_all.npz", allow_pickle=True)
            ok = (np.array_equal(a, z[f"{gene}__anchor"]) and np.array_equal(x, z[f"{gene}__xgb_pm"])
                  and np.array_equal(s, z[f"{gene}__xgb_pm__sd"]))
            print(f"{gene} {T} f{k}: {'identical' if ok else 'DIFFERS'}")
        m, s = pool(F)
        z = np.load(f"{A2}/lf_winxgb_scan_{T}_cv5all50{TAG}_all.npz", allow_pickle=True)
        ok = np.array_equal(m, z[f"{gene}__xgb_pm"]) and np.array_equal(s, z[f"{gene}__xgb_pm__sd"])
        print(f"{gene} {T} all-50: {'identical' if ok else 'DIFFERS'}")
elif mode == "build":
    sfx = sys.argv[4]
    if sfx in ("", "_all"): sys.exit("give the files their own suffix, never _all")
    for T in ("N", "C"):
        F = folds(T, idx)
        paths = [f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}{sfx}.npz" for k in range(1, K + 1)] + \
                [f"{A2}/lf_winxgb_scan_{T}_cv5all50{TAG}{sfx}.npz"]
        if any(os.path.exists(p) for p in paths): sys.exit(f"refusing to overwrite: {[p for p in paths if os.path.exists(p)]}")
        for (a, x, s), p in zip(F, paths):
            np.savez_compressed(p, **{f"{gene}__anchor": a, f"{gene}__xgb_pm": x, f"{gene}__xgb_pm__sd": s})
        m, s = pool(F)
        np.savez_compressed(paths[-1], **{f"{gene}__anchor": F[0][0], f"{gene}__xgb_pm": m, f"{gene}__xgb_pm__sd": s})
        print(f"{gene} ({acc}, prot_idx {idx}) {T}: {len(F[0][0])} anchors {F[0][0][0]}-{F[0][0][-1]}, "
              f"5 folds + all-50 -> *{sfx}.npz")
else:
    sys.exit("mode: build | check")
