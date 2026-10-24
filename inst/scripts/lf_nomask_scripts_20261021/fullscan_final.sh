#!/bin/zsh
# Full-proteome scan (every mature anchor, ~3M per terminus) with each 5-fold
# flatten run: CNN flatten head + per-member xgb head, exactly as cv5_fold.sh
# scanned the attention folds -- so ROC / normalised PR / recall-vs-rank compare
# like for like. Ten (fold, terminus) scans, four at a time.
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  OUT=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}_nomask_flatten.npz
  [[ -f $OUT ]] && { echo "$T f$k cached"; exit 0; }
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  echo "### $T f$k start $(date +%H:%M:%S)" > $S/fullscan_${T}_f$k.log
  python -m lf_winxgb.full_scan --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten.rds \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked --out $OUT >> $S/fullscan_${T}_f$k.log 2>&1
  echo "### $T f$k rc=$? $(date +%H:%M:%S)" >> $S/fullscan_${T}_f$k.log
  exit 0
fi
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P 4 -n 2 ${0:A} one
echo "ALL SCANS DONE $(date +%H:%M:%S)"
