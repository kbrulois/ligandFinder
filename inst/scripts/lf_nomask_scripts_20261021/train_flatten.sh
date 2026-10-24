#!/bin/zsh
# Position-aware window head: global_head = flatten, at the trunk's learning
# rate and at 10x it. 2-fold CV (foldA / foldB) x 5 seeds x both termini, no
# class-position mask -- so lf_pepend_run_*_fold*_nomask (attn head, same folds,
# same seeds) is the exact baseline. One arm per wave, four trainings at a time.
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
cd $W
for M in 1 10; do
  tag=_nomask_flatten; [[ $M -ne 1 ]] && tag=${tag}_hlr$M
  for F in A B; do
    for T in C N; do
      ( echo "### $T fold$F x$M start $(date +%H:%M:%S)"
        LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_8_pepend_train.R \
          --term $T --n-seeds 5 --no-mask --global-head flatten --head-lr-mult $M \
          --input $A2/lf_pepend_nn_input_fold${F}.rds \
          --out $A2/lf_pepend_run_${T}_fold${F}${tag}.rds
        echo "### $T fold$F x$M rc=$? $(date +%H:%M:%S)" ) > $S/trainflat_${T}_fold${F}_x$M.log 2>&1 &
    done
  done
  wait
  echo "wave x$M done $(date +%H:%M:%S)"
done
echo "ALL DONE $(date +%H:%M:%S)"
