#!/bin/zsh
# Full-proteome scan of the ESM-C-only COVERED-ONLY trunks (Kevin, 2026-10-24).
# The exact counterpart of fullscan_esmonly36.sh; only TAG and the log prefix
# differ.
#
# THE RESIDUE TABLE IS REUSED, and that is correct, not a shortcut.
# lf_pepend_residues_esmonly36.npz holds `padding` plus plm_01..plm_36 for every
# residue of every precursor. Those columns are a property of the EMBEDDING and
# the channel order, not of which precursors were in a training split -- and
# this arm differs from _esmonly36 only in its training split. The two arms have
# byte-identical channel_names, which esmonly36x_driver.sh asserts from the
# trained config.json before it scans, so build_windows' positional copy is
# safe. Rebuilding the table would produce the same 3,056,424 x 37 array.
#
# THE --groups NPZ, BY CONTRAST, IS *NOT* REUSABLE, and assuming it was cost a
# full scan queue. fullscan_esmonly36.sh reuses the curated
# lf_pepend_window_groups_<T>_cv5_f<k>.npz on the grounds that they "carry the
# fold/member structure, not the channels". True for an arm that differs only in
# channels; false here, because this arm differs in its TRAINING SET.
# full_scan.py line 85 does `acc_tr[yc["train"] == 1]`, indexing the groups
# table's accessions with the MODEL's own y_global, so the two must agree row
# for row -- and they did not:
#   IndexError: boolean index did not match indexed array along dimension 0;
#   dimension is 469 but corresponding boolean dimension is 461
# groups_esmcov.sh rebuilds them from the same _esmcov inputs used for training
# and asserts every split's length against the trained arrays.npz.
#
# One honest consequence of reusing the residue table: the 434 precursors with no embedding are
# still SCANNED, on all-zero input, by trunks that never saw an all-zero window
# in training. Their scores are out-of-distribution and should not be read as
# predictions. That is why the evaluation is also run restricted to the 5,182
# covered precursors (KEEP=, see segments_roc.R).
#
#   ./fullscan_esmonly36x.sh          -> the queue
#   ./fullscan_esmonly36x.sh one K T  -> one scan
#   PROCS=4 ./fullscan_esmonly36x.sh  -> raise the parallelism
set -u
W=/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/class-position-mask-diagram-b398c7
A2=~/AF2_analysis
S=${0:A:h}
TAG=_nomask_flatten_esmonly36x
PROCS=${PROCS:-2}
RESID=$A2/lf_pepend_residues_esmonly36.npz

if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  OUT=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}${TAG}.npz
  [[ -f $OUT ]] && { echo "$T f$k cached"; exit 0; }
  RUN=$A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds
  [[ -s $RUN ]] || { echo "$T f$k SKIPPED: no trunk at $(basename $RUN)"; exit 0; }
  [[ -s $RESID ]] || { echo "$T f$k SKIPPED: no $(basename $RESID)"; exit 1; }
  GRP=$A2/lf_pepend_window_groups_${T}_cv5_f${k}_esmcov.npz
  [[ -s $GRP ]] || { echo "$T f$k SKIPPED: no $(basename $GRP) -- run groups_esmcov.sh"; exit 1; }
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  echo "### $T f$k start $(date '+%F %H:%M:%S')" > $S/fsesm36x_${T}_f$k.log
  python -m lf_winxgb.full_scan --term $T \
    --run $RUN \
    --residues-npz $RESID \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}_esmcov.npz \
    --pm-rounds 100 --skip-stacked --out $OUT >> $S/fsesm36x_${T}_f$k.log 2>&1
  rc=$?
  echo "### $T f$k rc=$rc $(date '+%F %H:%M:%S')" >> $S/fsesm36x_${T}_f$k.log
  echo "scan $T f$k rc=$rc"
  exit 0
fi

missing=0
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds ]] || { echo "MISSING trunk $T f$k"; missing=1; }
done; done
[[ $missing -eq 0 ]] || { echo "train_esmonly36x.sh has not finished -- not scanning"; exit 1; }
[[ -s $RESID ]] || { echo "MISSING $(basename $RESID)"; exit 1; }
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_window_groups_${T}_cv5_f${k}_esmcov.npz ]] || { echo "MISSING groups $T f$k -- run groups_esmcov.sh"; missing=1; }
done; done
[[ $missing -eq 0 ]] || exit 1

echo "ESM-ONLY (COVERED) FULL SCAN start $(date '+%F %H:%M:%S')  ($PROCS at a time)"
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P $PROCS -n 2 ${0:A} one
echo "ESM-ONLY (COVERED) FULL SCAN DONE $(date '+%F %H:%M:%S')"
for k in 1 2 3 4 5; do for T in C N; do
  f=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}${TAG}.npz
  [[ -f $f ]] && echo "  $T f$k  $(du -h $f | cut -f1)" || echo "  $T f$k  MISSING"
done; done
