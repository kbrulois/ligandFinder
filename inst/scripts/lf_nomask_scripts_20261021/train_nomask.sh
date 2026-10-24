#!/bin/zsh
# No-mask arm: 2-fold CV (foldA / foldB, the same cluster partition as
# --n-folds 2), 5 seeds, both termini -- four trainings in parallel.
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
cd $W
for F in A B; do
  for T in C N; do
    ( echo "### $T fold$F start $(date +%H:%M:%S)"
      LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_8_pepend_train.R \
        --term $T --n-seeds 5 --no-mask \
        --input $A2/lf_pepend_nn_input_fold${F}.rds \
        --out $A2/lf_pepend_run_${T}_fold${F}_nomask.rds
      echo "### $T fold$F rc=$? $(date +%H:%M:%S)" ) > $S/train_${T}_fold${F}.log 2>&1 &
  done
done
wait
echo "ALL DONE $(date +%H:%M:%S)"
