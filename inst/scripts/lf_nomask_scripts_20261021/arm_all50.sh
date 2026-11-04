#!/bin/zsh
# An arm's all-50 ensemble: the five folds' 10 members each, linked into one
# run directory, the step cv5_post_esmonly36.sh does for ESM-C. What a page
# scores the arm's all-50 per-index surface from (final_pages.R ARMS_ONLY=1).
#
#   ./arm_all50.sh curated_AAswap esmonly36
#
# A directory already holding 50 members is left alone.
set -u
A2=~/AF2_analysis
for arm in $@; do
  TAG=_nomask_flatten_${arm}
  for T in C N; do
    E=$A2/lf_pepend_run_${T}_cv5all50${TAG}_isolated
    if (( $(print -l $E/out/members/member_*_${T}.weights.h5(N) | wc -l) == 50 )); then
      echo "$arm $T: 50 members already in $E"; continue
    fi
    mkdir -p $E/in $E/out/members
    rm -f $E/out/members/member_*_${T}.weights.h5(N)
    cp $A2/lf_pepend_run_${T}_cv5_f1${TAG}_isolated/in/config.json $E/in/
    i=0
    for k in 1 2 3 4 5; do
      for f in $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}_isolated/out/members/member_*_${T}.weights.h5; do
        ln -s $f $E/out/members/$(printf "member_%03d_%s.weights.h5" $i $T); i=$((i+1))
      done
    done
    (( i == 50 )) || { echo "$arm $T: $i members, want 50"; exit 1; }
    echo "$arm $T: $i members -> $E"
  done
done
