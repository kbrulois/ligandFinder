#!/bin/zsh
# Per-gene xgb/CNN window scans for the WHOLE secretome (Kevin, 2026-10-22:
# pages for all secretome proteins). Same call as cv5_post.sh step 3, but with
# all_genes.txt (5,607 genes) instead of the 104 page genes, and its own `_all`
# output so the verified 104-gene npz set is never clobbered.
#
# The CNN per-member sd is why this has to run at all: lf_winxgb_fullscan_*.npz
# already covers every anchor, but carries xgb_pm_sd only -- the pages' spread
# band needs cnn__sd, which only scan_genes emits.
#
# Resumable: a job whose .npz already exists is skipped, so re-running after an
# interruption picks up where it stopped.
#   allsec_scan.sh            -> the 10 fold/terminus jobs, 4 at a time, then pool
#   allsec_scan.sh one K T    -> one job (what xargs calls)
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
TAG=_nomask_flatten
SFX=_all

if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  OUT=$A2/lf_winxgb_scan_${T}_cv5_f${k}${TAG}${SFX}.npz
  if [[ -f $OUT ]]; then echo "scan $T f$k cached"; exit 0; fi
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  echo "### $T f$k start $(date +%H:%M:%S)" > $S/allscan_${T}_f$k.log
  python -m lf_winxgb.scan_genes --genes $(paste -sd, $S/all_genes.txt) --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked \
    --out $OUT >> $S/allscan_${T}_f$k.log 2>&1
  rc=$?
  echo "### $T f$k rc=$rc $(date +%H:%M:%S)" >> $S/allscan_${T}_f$k.log
  echo "scan $T f$k rc=$rc"
  exit 0
fi

echo "ALL-SECRETOME SCAN start $(date '+%F %H:%M:%S')  ($(wc -l < $S/all_genes.txt) genes)"
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P 4 -n 2 ${0:A} one

# every job must have produced its npz before pooling
missing=0
for T in C N; do for k in 1 2 3 4 5; do
  [[ -f $A2/lf_winxgb_scan_${T}_cv5_f${k}${TAG}${SFX}.npz ]] || { echo "MISSING $T f$k"; missing=1; }
done; done
[[ $missing -eq 0 ]] || { echo "SCAN INCOMPLETE -- rerun allsec_scan.sh to resume"; exit 1; }

# pool the 5 folds into the all-50 ensemble, exactly as cv5_post.sh does:
# mean over all 50 members, sd recovered from each fold's 10-member mean and sd
source ~/.virtualenvs/r-tensorflow/bin/activate
TAG=$TAG SFX=$SFX GENES_FILE=$S/all_genes.txt python -I - <<'EOF'
import numpy as np, os
A2 = os.path.expanduser("~/AF2_analysis"); TAG = os.environ["TAG"]; SFX = os.environ["SFX"]
GENES = [g.strip() for g in open(os.environ["GENES_FILE"]) if g.strip()]
n, K = 10, 5
for T in ("C", "N"):
    z = [np.load(f"{A2}/lf_winxgb_scan_{T}_cv5_f{k}{TAG}{SFX}.npz", allow_pickle=True) for k in range(1, K + 1)]
    out, miss = {}, []
    for g in GENES:
        if f"{g}__anchor" not in z[0].files: miss.append(g); continue
        an = z[0][f"{g}__anchor"]
        assert all(np.array_equal(an, zz[f"{g}__anchor"]) for zz in z), g
        out[f"{g}__anchor"] = an
        for h in ("cnn", "xgb_pm"):
            M = np.stack([zz[f"{g}__{h}"].astype("f8") for zz in z])
            Sd = np.stack([zz[f"{g}__{h}__sd"].astype("f8") for zz in z])
            m = M.mean(0)
            var = ((n - 1) * (Sd ** 2).sum(0) + n * ((M - m) ** 2).sum(0)) / (n * K - 1)
            out[f"{g}__{h}"] = m.astype("f4"); out[f"{g}__{h}__sd"] = np.sqrt(var).astype("f4")
    np.savez_compressed(f"{A2}/lf_winxgb_scan_{T}_cv5all50{TAG}{SFX}.npz", **out)
    print(f"{T}: {len(out)//5} genes pooled over 5 folds x 10 members" +
          (f"; {len(miss)} not in the scan" + (f": {', '.join(miss[:8])}..." if miss else "") if miss else ""))
EOF
echo "ALL-SECRETOME SCAN DONE $(date '+%F %H:%M:%S')"
