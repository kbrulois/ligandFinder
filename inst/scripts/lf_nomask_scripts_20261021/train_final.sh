#!/bin/zsh
# THE model going forward (Kevin, 2026-10-21): the position-aware flatten window
# head at the trunk's learning rate, no class-position mask, on the 5-fold CV
# sets (lf_pepend_nn_input_cv5_f1..5), 10 seeds per fold per terminus = 100
# models. Ten (fold, terminus) jobs, four at a time (xargs -P 4).
#   train_final.sh            -> runs the queue
#   train_final.sh one K T    -> one job (what xargs calls)
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  {
    echo "### $T f$k start $(date +%H:%M:%S)"
    cd $W && LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_8_pepend_train.R \
      --term $T --n-seeds 10 --no-mask --global-head flatten \
      --input $A2/lf_pepend_nn_input_cv5_f$k.rds \
      --out $A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten.rds
    echo "### $T f$k rc=$? $(date +%H:%M:%S)"
  } > $S/final_${T}_f$k.log 2>&1
  exit 0
fi
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done |
  xargs -P 4 -n 2 ${0:A} one
echo "ALL DONE $(date +%H:%M:%S)"
