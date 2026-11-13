#!/bin/zsh
# One CV fold, end to end. Idempotent: trunk members and scans cache on disk,
# so re-running resumes instead of restarting.
set -u
k=$1
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8
S2=/private/tmp/claude-501/-Users-kbrulois-R-projects-ligandFinder--claude-worktrees-modest-satoshi-87b2c8/92c0469f-235a-4310-b3a9-55060cf86b72/scratchpad
A2=~/AF2_analysis
IN=$A2/lf_pepend_nn_input_cv5_f${k}.rds
cd $W
echo "### fold $k start $(date +%H:%M:%S)"
for T in C N; do
  [[ -f $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz ]] || \
    LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/modeling/10_9r_export_window_groups.R \
      --term $T --input $IN --out $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz >/dev/null 2>&1
  echo "### fold $k $T train $(date +%H:%M:%S)"
  LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/modeling/10_8_pepend_train.R \
    --term $T --n-seeds 10 --input $IN \
    --out $A2/lf_pepend_run_${T}_cv5_f${k}.rds >> $S2/cv5_f${k}_${T}_train.log 2>&1
  echo "### fold $k $T train rc=$? $(date +%H:%M:%S)"
done
source ~/.virtualenvs/r-tensorflow/bin/activate
cd $W/inst/python
for T in C N; do
  OUT=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}.npz
  if [[ -f $OUT ]]; then echo "### fold $k $T scan cached"; continue; fi
  echo "### fold $k $T fullscan $(date +%H:%M:%S)"
  python -m lf_winxgb.full_scan --term $T \
    --run $A2/lf_pepend_run_${T}_cv5_f${k}_slack8.rds \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked --out $OUT >> $S2/cv5_f${k}_${T}_scan.log 2>&1
  echo "### fold $k $T fullscan rc=$? $(date +%H:%M:%S)"
done
echo "### fold $k DONE $(date +%H:%M:%S)"
