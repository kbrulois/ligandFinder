#!/bin/zsh
set -u
export LC_ALL=en_US.UTF-8
ROOT=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8
A2=~/AF2_analysis
IN=$A2/lf_pepend_nn_input_near.rds
## every arm needs its OWN --out: lf_dcnn_run_isolated reuses any member already
## sitting in <out>_isolated/out/members/ on presence alone
run () {
  tag=$1; cfg=$2
  echo "=== $tag ($cfg) $(date +%H:%M:%S)"
  /usr/bin/time -p Rscript $ROOT/inst/scripts/modeling/10_8_pepend_train.R \
      --term C --input $IN --n-seeds 5 --cfg "$cfg" \
      --out $A2/lf_pepend_run_C_near_slack8_$tag.rds 2>&1
  echo "=== $tag done $(date +%H:%M:%S)"
}
run cbamtrunk "cbam=trunk"
run cbamhead  "cbam=head"
run dice      "per_index_loss=dice"
echo "ALL ARMS DONE $(date +%H:%M:%S)"
