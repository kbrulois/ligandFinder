#!/bin/zsh
## ---- hyperparameter search for the peptide-end window CNN ---------------------
## Learning rate, U-Net width (filters per level), kernel size and early stopping
## (start epoch, patience, epoch budget, monitor), searched over the five cv5
## folds of the model train_final.sh trains -- the flatten head, no class-position
## mask, one terminus. The method, and what is held out from what, is in the
## docstring of inst/python/lf_dcnn/tune.py; read it before reading the numbers.
##
##   inst/scripts/10_8t_pepend_tune.sh                     # everything, terminus C
##   STAGES=run,report inst/scripts/10_8t_pepend_tune.sh   # resume
##   STAGES=report inst/scripts/10_8t_pepend_tune.sh       # re-aggregate only
##   EXTEND="5:2,3" STAGES=extend,run,report ...           # seeds 2,3 for the top 5
##   N_CONFIGS=32 WORKERS=6 THREADS=2 ROOT=~/AF2_analysis/lf_tune_C2 ...
##   EXTRA="--plm ~/AF2_analysis/lf_plm/aa_onehot.parquet" ROOT=... ...   # another arm
##
## Stages, each resumable:
##   export  10_8_pepend_train.R --export-only, once per fold (skipped if done)
##   plan    python -m lf_dcnn tune plan            (skipped if plan.json exists)
##   extend  add seeds to the top configs ($EXTEND = "<top>:<seeds>")
##   run     python -m lf_dcnn tune run             (finished trials are skipped)
##   report  python -m lf_dcnn tune report          (always re-runs; cheap)
##
## A trial is one model (config x fold x seed) in its own process; each answers
## every early-stopping setting at once. ~4-7 min per trial with 4 in parallel on
## the M2 Max, so the default 25 configs x 5 folds x 1 seed is ~3-4 hours.
## Stop it with Ctrl-C (or kill the `tune run` process and its `tune trial`
## children); re-running picks up where it stopped.
##
## The winner trains through the normal pipeline with the flags the report
## prints, e.g. `10_8_pepend_train.R ... --learning-rate 1e-3 --conv-kernel 5`.
## ------------------------------------------------------------------------------
set -u
W=${0:A:h:h:h}                                   # the checkout this script sits in
A2=${A2:-$HOME/AF2_analysis}
T=${T:-C}                                        # not TERM: that is the terminal's
ROOT=${ROOT:-$A2/lf_tune_$T}
## the production model's flags, and any arm-specific channel flags on top
BASE_FLAGS=${BASE_FLAGS:---no-mask --global-head flatten}
EXTRA=${EXTRA:-}
FOLD_INPUT=${FOLD_INPUT:-$A2/lf_pepend_nn_input_cv5_f%d.rds}
FOLDS=${FOLDS:-1 2 3 4 5}
N_CONFIGS=${N_CONFIGS:-24}
SEEDS=${SEEDS:-1}
PARTS=${PARTS:-4}
WORKERS=${WORKERS:-4}
THREADS=${THREADS:-}
OTHER_ENDS=${OTHER_ENDS:-$A2/pepend_locus}
EXTEND=${EXTEND:-5:2,3}
STAGES=${STAGES:-export,plan,run,report}
PY=${PY:-$HOME/.virtualenvs/r-tensorflow/bin/python}
RSCRIPT=${RSCRIPT:-/usr/local/bin/Rscript}

[[ -x $PY ]] || { print -u2 "no python at $PY"; exit 1 }
mkdir -p $ROOT
## the checkout's lf_dcnn, never an installed one (lf-dcnn-path-shadowing)
lf() { PYTHONPATH=$W/inst/python $PY -m lf_dcnn tune "$@" }
has() { [[ ,$STAGES, == *,$1,* ]] }

echo "tune: terminus $T -> $ROOT  (stages $STAGES)"

if has export; then
  for k in ${=FOLDS}; do
    d=$ROOT/in/f$k
    if [[ -s $d/arrays.npz && -s $d/val_meta.parquet ]]; then
      echo "  export f$k: done"
      continue
    fi
    in=$(printf $FOLD_INPUT $k)
    echo "  export f$k <- $in"
    ## ${=VAR}: zsh does not word-split a plain $VAR, and a flag string passed
    ## as ONE argument is silently ignored by 10_8's .opt() (see train_arm.sh)
    (cd $W && LC_ALL=en_US.UTF-8 $RSCRIPT inst/scripts/10_8_pepend_train.R \
        --term $T ${=BASE_FLAGS} ${=EXTRA} --input $in --export-only $d) \
      > $ROOT/export_f$k.log 2>&1 || { print -u2 "export f$k failed; see $ROOT/export_f$k.log"; exit 1 }
    tail -1 $ROOT/export_f$k.log
  done
fi

if has plan; then
  if [[ -s $ROOT/plan.json ]]; then
    echo "  plan: exists"
  else
    oe=()
    [[ -s $OTHER_ENDS/ligand_list.csv ]] && oe=(--other-ends $OTHER_ENDS)
    lf plan --root $ROOT --inputs $ROOT/in --n-configs $N_CONFIGS --seeds $SEEDS \
      --parts $PARTS $oe || exit 1
  fi
fi

if has extend; then
  lf extend --root $ROOT --top ${EXTEND%%:*} --seeds ${EXTEND#*:} || exit 1
fi

if has run; then
  th=()
  [[ -n $THREADS ]] && th=(--threads $THREADS)
  lf run --root $ROOT --workers $WORKERS $th
fi

if has report; then
  lf report --root $ROOT
fi
