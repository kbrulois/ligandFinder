#!/bin/zsh
# Train one ARM's 10 trunks: 5 CV folds x {N, C}, 10 seeds each, flatten head,
# no class-position mask (Kevin, 2026-10-25). The arm is a PARAMETER -- its
# channel set lives in arms/$ARM.sh -- so a new channel set is a new 20-line
# file, not a new copy of this script.
#
#   ARM=curated_wAAid ./train_arm.sh        -> the queue (resumable)
#   ARM=curated_wAAid ./train_arm.sh one 3 C   -> one job (what xargs calls)
#   PROCS=4 ARM=... ./train_arm.sh          -> raise the parallelism
#
# WHY THE EXPLICIT --out MATTERS, and this is the trap that would silently ruin
# a whole arm: 10_8 derives its cache name from --plm but NOT from
# --drop-channels. Two arms that differ only in dropped channels would share a
# cache, and lf_dcnn_run_isolated reuses any member_*.npz already sitting in
# <cache>_isolated/out/members/ on its presence alone -- so one would train once
# and report as both. The ${ARM} in the name is what keeps them apart.
#
# PROCS defaults to 4, matching train_final.sh. This machine is memory-bound
# rather than CPU-bound; drop it if something large is rendering at the time.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
R="LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript"
PY_PLM=~/.virtualenvs/lf-plm/bin/python
PY_TF=~/.virtualenvs/r-tensorflow/bin/python
## export, not just set: xargs re-enters this script as `$0 one k T`, and a
## child that fell back to the default would train the WRONG ARM into the right
## file name.
export ARM=${ARM:-curated_wAAid}
PROCS=${PROCS:-4}
TAG=_nomask_flatten_${ARM}

[[ -r $S/arms/$ARM.sh ]] || { echo "no arm definition at arms/$ARM.sh"; exit 1; }
source $S/arms/$ARM.sh

## ZSH TRAP, and it cost a whole mis-trained arm: zsh does NOT word-split an
## unquoted parameter expansion, so `${DROP:+--drop-channels $DROP}` is ONE
## argument ("--drop-channels AA_hydro,...") and 10_8's .opt() never finds the
## flag -- it drops nothing and trains the wrong channel set under the right
## file name. `${DROP:+--drop-channels} $DROP` is two arguments when DROP is set
## and none when it is empty. Verified both ways with argc, not by eye.
##
## This is why the channel-count assertion below is not optional, and why it is
## now ALSO checked on the first trunk's config rather than only after all ten.
if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  out=$A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds
  {
    echo "### $T f$k start $(date '+%F %H:%M:%S')"
    cd $W && LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_8_pepend_train.R \
      --term $T --n-seeds 10 --no-mask --global-head flatten \
      --plm $PLM ${DROP:+--drop-channels} $DROP \
      --input $A2/lf_pepend_nn_input_cv5_f$k.rds \
      --out $out
    echo "### $T f$k rc=$? $(date '+%F %H:%M:%S')"
  } >> $S/train_${ARM}_${T}_f$k.log 2>&1
  exit 0
fi

[[ -s $PLM ]] || { echo "MISSING $PLM -- run stage 0 (mk_plm) first"; exit 1; }
echo "TRAINING arm=$ARM ($ARM_LABEL) start $(date '+%F %H:%M:%S')"
echo "  plm      : $(basename $PLM) ($(du -h $PLM | cut -f1))"
if [[ -n $DROP ]]; then
  echo "  channels : expect $N_CHAN_EXPECT; $(echo $DROP | tr ',' '\n' | wc -l | tr -d ' ') curated dropped"
else
  echo "  channels : expect $N_CHAN_EXPECT; all 26 curated kept, nothing dropped"
fi
echo "  jobs     : 5 folds x {N,C} = 10, $PROCS at a time, 10 seeds each"

# resumable: a job whose .rds already exists is skipped
jobs=()
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds ]] || jobs+=("$k $T")
done; done
if (( ${#jobs} == 0 )); then
  echo "  nothing to train"
else
  echo "  to train  : ${#jobs}"
  ## EARLY WIRING CHECK. 10_8 writes in/config.json at the start of a run, long
  ## before any trunk finishes, so the channel count is knowable within minutes.
  ## Getting it wrong means every trunk in the queue is the wrong channel set
  ## under the right file name -- which has happened. This reports; it does NOT
  ## kill, because the only pattern broad enough to stop the queue would also
  ## stop another session's training.
  (
    for _ in {1..90}; do
      ## (N) is zsh's null-glob qualifier: without it a glob that matches
      ## nothing is an ERROR, which is what this check printed on its first run.
      cj=$(print -rl -- $A2/lf_pepend_run_[NC]_cv5_f*${TAG}_isolated/in/config.json(N) | head -1)
      if [[ -n $cj ]]; then
        n=$($PY_PLM -c "import json,sys;print(len(json.load(open(sys.argv[1]))['channel_names']))" $cj 2>/dev/null)
        if [[ -n $n && $n -ne $N_CHAN_EXPECT ]]; then
          echo ""
          echo "!!! WIRING ERROR: trunks are training with $n channels, arms/$ARM.sh"
          echo "!!! expects $N_CHAN_EXPECT. Every job in this queue is the wrong channel set."
          echo "!!! Stop it, fix arms/$ARM.sh, delete the *${TAG}_isolated dirs and re-run."
        elif [[ -n $n ]]; then
          echo "  early check: $n channels, as expected"
        fi
        break
      fi
      sleep 10
    done
  ) &
  printf '%s\n' $jobs | xargs -P $PROCS -n 2 ${0:A} one
  wait
fi
echo "TRAINING DONE arm=$ARM $(date '+%F %H:%M:%S')"

## The channel count is the one thing worth asserting after the fact: a wrong
## parquet or a stale DROP list trains happily and only shows up as nonsense ten
## hours later, at the scan. Read it from a trunk's own config.
for k in 1 2 3 4 5; do for T in C N; do
  f=$A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds
  [[ -s $f ]] && echo "  $T f$k  $(du -h $f | cut -f1)" || echo "  $T f$k  MISSING"
done; done
cj=$A2/lf_pepend_run_C_cv5_f1${TAG}_isolated/in/config.json
if [[ -s $cj ]]; then
  n=$($PY_PLM -c "import json,sys; print(len(json.load(open(sys.argv[1]))['channel_names']))" $cj)
  if [[ $n -eq $N_CHAN_EXPECT ]]; then echo "  channels : $n, as expected"
  else echo "  CHANNELS : $n but arms/$ARM.sh expects $N_CHAN_EXPECT -- do not scan these"; exit 1
  fi
fi
