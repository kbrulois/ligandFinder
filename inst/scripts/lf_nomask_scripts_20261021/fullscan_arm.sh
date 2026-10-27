#!/bin/zsh
# Full-proteome scan of one ARM's trunks (Kevin, 2026-10-25): every mature
# anchor, ~3M per terminus, CNN flatten head + per-member xgb head in ONE pass.
# Both heads land in the same npz (`cnn` and `xgb_pm`), which is why "flat and
# xgb heads" needs no separate training step -- ten scans give all four models.
#
#   ARM=curated_wAAid ./fullscan_arm.sh         -> the queue (resumable)
#   ARM=curated_wAAid ./fullscan_arm.sh one 3 C -> one scan (what xargs calls)
#   PROCS=4 ARM=... ./fullscan_arm.sh           -> raise the parallelism
#
# The --groups npz are REUSED from the curated run: they carry the fold/member
# structure, not the channels, so the same ten files serve every arm.
#
# THE RESIDUE TABLE MUST BE THIS ARM'S OWN. build_windows copies res.feat
# straight into a (k, seq_len, n_channels) array POSITIONALLY: a wrong-width
# table raises, but a right-width table in the wrong order scans silently
# against the wrong channels and yields a full proteome of plausible nonsense.
# check_resid_order.py is run before the queue starts for exactly that reason.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
PY_PLM=~/.virtualenvs/lf-plm/bin/python
PY_TF=~/.virtualenvs/r-tensorflow/bin/python
export ARM=${ARM:-curated_wAAid}      # exported: xargs re-enters this script
PROCS=${PROCS:-4}
TAG=_nomask_flatten_${ARM}

[[ -r $S/arms/$ARM.sh ]] || { echo "no arm definition at arms/$ARM.sh"; exit 1; }
source $S/arms/$ARM.sh
RESID=${RESID:-$A2/lf_pepend_residues_${ARM}.npz}

if [[ ${1:-} == one ]]; then
  k=$2 T=$3
  OUT=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}${TAG}.npz
  [[ -f $OUT ]] && { echo "$T f$k cached"; exit 0; }
  RUN=$A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds
  [[ -s $RUN ]] || { echo "$T f$k SKIPPED: no trunk at $(basename $RUN)"; exit 0; }
  source ~/.virtualenvs/r-tensorflow/bin/activate
  cd $W/inst/python
  echo "### $T f$k start $(date '+%F %H:%M:%S')" > $S/fs_${ARM}_${T}_f$k.log
  python -m lf_winxgb.full_scan --term $T \
    --run $RUN \
    --residues-npz $RESID \
    --groups $A2/lf_pepend_window_groups_${T}_cv5_f${k}.npz \
    --pm-rounds 100 --skip-stacked --out $OUT >> $S/fs_${ARM}_${T}_f$k.log 2>&1
  rc=$?
  echo "### $T f$k rc=$rc $(date '+%F %H:%M:%S')" >> $S/fs_${ARM}_${T}_f$k.log
  echo "scan $T f$k rc=$rc"
  exit 0
fi

# every trunk must exist before scanning, or the queue silently half-fills
missing=0
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}${TAG}.rds ]] || { echo "MISSING trunk $T f$k"; missing=1; }
done; done
[[ $missing -eq 0 ]] || { echo "training has not finished -- not scanning"; exit 1; }
[[ -s $RESID ]] || { echo "MISSING $(basename $RESID) -- run stage 2 first"; exit 1; }

## the check that stands between a mis-ordered table and a proteome of nonsense
$PY_PLM $S/check_resid_order.py --arm $ARM --plm $PLM --residues-npz $RESID \
  || { echo "residue table does not match the trunks -- not scanning"; exit 1; }

echo "FULL SCAN arm=$ARM start $(date '+%F %H:%M:%S')  ($PROCS at a time)"
for k in 1 2 3 4 5; do for T in C N; do echo "$k $T"; done; done | xargs -P $PROCS -n 2 ${0:A} one
echo "FULL SCAN DONE arm=$ARM $(date '+%F %H:%M:%S')"
for k in 1 2 3 4 5; do for T in C N; do
  f=$A2/lf_winxgb_fullscan_${T}_cv5_f${k}${TAG}.npz
  [[ -f $f ]] && echo "  $T f$k  $(du -h $f | cut -f1)" || echo "  $T f$k  MISSING"
done; done
