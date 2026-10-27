#!/bin/zsh
# The ESM-C-only arm with the UNCOVERED PRECURSORS EXCLUDED (Kevin, 2026-10-24).
#
# Identical to train_esmonly36.sh in every respect -- same 36-channel PCA, same
# 25 dropped curated channels, --no-mask --global-head flatten, 10 seeds, the
# same 5 folds -- except for ONE thing: --input points at the _esmcov inputs,
# whose train and val splits hold only the 5,182 precursors that actually have
# an ESM C embedding. nn_input_esmcov.R made those and reports what it removed.
#
# WHY: in this arm a precursor with no embedding has an all-zero 36-channel
# block and therefore carries no information whatever, so those windows are
# label noise. The question this arm answers is whether the never-seen-gene
# advantage of _esmonly36 survives their removal.
#
# WHAT IT IS NOT: a 7.8% cut. 7.8% (434/5,616) is the share of the SCANNED
# PROTEOME with no embedding. The share of FITTING windows is 0.76-2.15% per
# split, because training windows sit on known peptide precursors, which are
# well-studied and nearly all embedded. See nn_input_esmcov_report.csv for the
# exact per-fold counts. Expect a small effect; measure it anyway.
#
#   train_esmonly36x.sh            -> runs the queue
#   train_esmonly36x.sh one K T    -> one job (what xargs calls)
#   PROCS=4 train_esmonly36x.sh    -> raise the parallelism
#
# THE --out TRAP, restated because it would silently ruin this arm: 10_8 derives
# its cache name from --plm but NOT from --drop-channels and NOT from --input,
# and lf_dcnn_run_isolated reuses any member_*.npz already in
# <cache>_isolated/out/members/ on presence alone. This arm differs from
# _esmonly36 ONLY in --input, which the cache name cannot see at all -- so
# without the distinct _esmonly36x in --out, this script would find
# _esmonly36's members, train nothing, and report the OLD model as the new arm.
# The suffix is the only thing keeping them apart. Do not remove it.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
PLM=$A2/lf_plm/esm_c_pca36.parquet
PROCS=${PROCS:-2}

# the 25 curated channels, i.e. all of all_params3 except `padding`
DROP='cons_rs_n,min_afm,mean_afm,relASA,SS_P,SS_S,SS_E,SS_-,SS_T,SS_G,SS_B,SS_H,SS_I,Phi_cos,Psi_cos,Phi_sin,Psi_sin,NH->O_1_energy,O->NH_1_energy,NH->O_2_energy,O->NH_2_energy,AA_hydro,AA_charge,AA_mw,AA_pI'

if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  out=$A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36x.rds
  {
    echo "### $T f$k start $(date '+%F %H:%M:%S')"
    cd $W && LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_8_pepend_train.R \
      --term $T --n-seeds 10 --no-mask --global-head flatten \
      --plm $PLM --drop-channels $DROP \
      --input $A2/lf_pepend_nn_input_cv5_f${k}_esmcov.rds \
      --out $out
    echo "### $T f$k rc=$? $(date '+%F %H:%M:%S')"
  } >> $S/trainesm36x_${T}_f$k.log 2>&1
  exit 0
fi

[[ -s $PLM ]] || { echo "MISSING $PLM"; exit 1; }
for k in 1 2 3 4 5; do
  [[ -s $A2/lf_pepend_nn_input_cv5_f${k}_esmcov.rds ]] || { echo "MISSING _esmcov input f$k -- run nn_input_esmcov.R first"; exit 1; }
done
echo "ESM-C-ONLY (COVERED PRECURSORS ONLY) TRAINING start $(date '+%F %H:%M:%S')"
echo "  plm      : $(basename $PLM)"
echo "  input    : lf_pepend_nn_input_cv5_f<k>_esmcov.rds (train/val restricted to embedded precursors)"
echo "  jobs     : 5 folds x {N,C} = 10, $PROCS at a time, 10 seeds each"

jobs=()
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36x.rds ]] || jobs+=("$k $T")
done; done
if (( ${#jobs} == 0 )); then
  echo "  nothing to train"
else
  echo "  to train  : ${#jobs}"
  printf '%s\n' $jobs | xargs -P $PROCS -n 2 ${0:A} one
fi
echo "ESM-C-ONLY (COVERED) TRAINING DONE $(date '+%F %H:%M:%S')"
for k in 1 2 3 4 5; do for T in C N; do
  f=$A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36x.rds
  [[ -s $f ]] && echo "  $T f$k  $(du -h $f | cut -f1)" || echo "  $T f$k  MISSING"
done; done
