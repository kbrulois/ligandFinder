#!/bin/zsh
# Full-proteome scan of the ESM-C-only trunks (Kevin, 2026-10-24), the exact
# counterpart of fullscan_final.sh: every mature anchor, ~3M per terminus, CNN
# flatten head + per-member xgb head in one pass. Both heads land in the same
# npz (`cnn` and `xgb_pm`), which is why "flat and xgb heads" needs no separate
# training step -- ten scans give all four models.
#
#   ./fullscan_esmonly36.sh          -> the queue
#   ./fullscan_esmonly36.sh one K T  -> one scan (what xargs calls)
#   PROCS=4 ./fullscan_esmonly36.sh  -> raise the parallelism
#
# The --groups npz are REUSED from the curated run: they carry the fold/member
# structure, not the channels, so the same ten files serve both arms.
#
# PROCS defaults to 2, not the 4 of fullscan_final.sh, because this machine is
# MEMORY-bound, not CPU-bound: with the all-secretome render going, free RAM sits
# near zero and swap near its allocation, and adding load halves the render's
# throughput rather than sharing it. Raise this once the render is done.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
ARM=${ARM:-esmonly36}
TAG=_nomask_flatten_${ARM}
PROCS=${PROCS:-2}
# The residue table MUST be this arm's own 37-channel one. build_windows copies
# res.feat straight into a (k, seq_len, n_channels) array positionally, so the
# default 26-channel table does not merely mismatch -- a wrong-width table
# raises, and a right-width table in the wrong order would scan silently
# against the wrong channels. residues_esmonly36.py asserts the order against
# the trained model's channel_names.
RESID=${RESID:-$A2/lf_pepend_residues_${ARM}.npz}

if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  OUT=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}${TAG}.npz
  [[ -f $OUT ]] && { echo "$T f$k cached"; exit 0; }
  RUN=$A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds
  [[ -s $RUN ]] || { echo "$T f$k SKIPPED: no trunk at $(basename $RUN)"; exit 0; }
  [[ -s $RESID ]] || { echo "$T f$k SKIPPED: no $(basename $RESID), run residues_esmonly36.py"; exit 1; }
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  echo "### $T f$k start $(date '+%F %H:%M:%S')" > $S/fsesm36_${T}_f$k.log
  python -m lf_winxgb.full_scan --term $T \
    --run $RUN \
    --residues-npz $RESID \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked --out $OUT >> $S/fsesm36_${T}_f$k.log 2>&1
  rc=$?
  echo "### $T f$k rc=$rc $(date '+%F %H:%M:%S')" >> $S/fsesm36_${T}_f$k.log
  echo "scan $T f$k rc=$rc"
  exit 0
fi

# every trunk must exist before scanning, or the queue silently half-fills
missing=0
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds ]] || { echo "MISSING trunk $T f$k"; missing=1; }
done; done
[[ $missing -eq 0 ]] || { echo "train_esmonly36.sh has not finished -- not scanning"; exit 1; }
[[ -s $RESID ]] || { echo "MISSING $(basename $RESID) -- run residues_esmonly36.py first"; exit 1; }

echo "ESM-ONLY FULL SCAN start $(date '+%F %H:%M:%S')  ($PROCS at a time)"
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P $PROCS -n 2 ${0:A} one
echo "ESM-ONLY FULL SCAN DONE $(date '+%F %H:%M:%S')"
for k in 1 2 3 4 5; do for T in C N; do
  f=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}${TAG}.npz
  [[ -f $f ]] && echo "  $T f$k  $(du -h $f | cut -f1)" || echo "  $T f$k  MISSING"
done; done
