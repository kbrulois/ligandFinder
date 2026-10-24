#!/bin/zsh
# The ESM-C-only arm (Kevin, 2026-10-24): the flatten window head trained on the
# 36-component PCA of the ESM C per-residue embedding and NOTHING ELSE. All 25
# curated channels are dropped; only `padding` survives beside the 36 ESM
# channels, because padding is the mask channel and 10_8 refuses to drop it.
# So the trunk sees 37 channels, of which 36 carry information.
#
# Mirrors train_final.sh in every other respect -- the same 5-fold CV inputs,
# 10 seeds per fold per terminus, --no-mask --global-head flatten -- so that the
# end-segment comparison at THRESH=0.5 / TOL=4 is like-for-like against
# lf_pepend_run_<T>_cv5_f<k>_nomask_flatten.rds.
#
#   train_esmonly36.sh            -> runs the queue
#   train_esmonly36.sh one K T    -> one job (what xargs calls)
#   PROCS=4 train_esmonly36.sh    -> raise the parallelism
#
# WHY THE EXPLICIT --out MATTERS, and this is the trap that would silently ruin
# the whole arm: 10_8 derives its cache name from --plm but NOT from
# --drop-channels. An ESM-only run and an ESM+curated run would therefore share
# a cache, and lf_dcnn_run_isolated reuses any member_*.npz already sitting in
# <cache>_isolated/out/members/ on its presence alone -- so two arms sharing an
# --out train once and report it as both. The _esmonly36 in the name is what
# keeps these separate. Do not remove it.
#
# PROCS defaults to 2, not the 4 of train_final.sh, because the all-secretome
# page render is using five R workers on this machine. Raise it once that ends.
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
PLM=$A2/lf_plm/esm_c_pca36.parquet
PROCS=${PROCS:-2}

# the 25 curated channels, i.e. all of all_params3 except `padding`
DROP='cons_rs_n,min_afm,mean_afm,relASA,SS_P,SS_S,SS_E,SS_-,SS_T,SS_G,SS_B,SS_H,SS_I,Phi_cos,Psi_cos,Phi_sin,Psi_sin,NH->O_1_energy,O->NH_1_energy,NH->O_2_energy,O->NH_2_energy,AA_hydro,AA_charge,AA_mw,AA_pI'

if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  out=$A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36.rds
  {
    echo "### $T f$k start $(date '+%F %H:%M:%S')"
    cd $W && LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript inst/scripts/10_8_pepend_train.R \
      --term $T --n-seeds 10 --no-mask --global-head flatten \
      --plm $PLM --drop-channels $DROP \
      --input $A2/lf_pepend_nn_input_cv5_f$k.rds \
      --out $out
    echo "### $T f$k rc=$? $(date '+%F %H:%M:%S')"
  } >> $S/trainesm36_${T}_f$k.log 2>&1
  exit 0
fi

[[ -s $PLM ]] || { echo "MISSING $PLM -- run lf_plm reduce --k 36 first"; exit 1; }
echo "ESM-C-ONLY TRAINING start $(date '+%F %H:%M:%S')"
echo "  plm      : $(basename $PLM) ($(du -h $PLM | cut -f1))"
echo "  channels : 36 ESM + padding; $(echo $DROP | tr ',' '\n' | wc -l | tr -d ' ') curated dropped"
echo "  jobs     : 5 folds x {N,C} = 10, $PROCS at a time, 10 seeds each"

# resumable: a job whose .rds already exists is skipped
jobs=()
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36.rds ]] || jobs+=("$k $T")
done; done
if (( ${#jobs} == 0 )); then
  echo "  nothing to train"
else
  echo "  to train  : ${#jobs}"
  printf '%s\n' $jobs | xargs -P $PROCS -n 2 ${0:A} one
fi
echo "ESM-C-ONLY TRAINING DONE $(date '+%F %H:%M:%S')"
for k in 1 2 3 4 5; do for T in C N; do
  f=$A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36.rds
  [[ -s $f ]] && echo "  $T f$k  $(du -h $f | cut -f1)" || echo "  $T f$k  MISSING"
done; done
