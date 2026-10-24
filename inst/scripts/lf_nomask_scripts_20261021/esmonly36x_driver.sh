#!/bin/zsh
# Unattended driver for the ESM-C-only COVERED-ONLY arm (Kevin, 2026-10-24).
#
#   nohup ./esmonly36x_driver.sh > esmonly36x_driver.log 2>&1 &
#
# Mirrors esmonly36_driver.sh's shape -- gate, one cheap proof job, then the
# rest -- and adds the scan/tracks/evaluation tail so the whole chain runs
# unattended. Every stage is gated on the previous one having actually worked.
#
#   0  wait for the all-secretome page render to finish. This Mac is MEMORY-
#      bound, not CPU-bound: free RAM sits near zero and swap near its
#      allocation, so training alongside five render workers costs the render
#      ~35% of its throughput and risks thrashing both. Waiting is faster.
#   1  one job (C fold 1) to prove the path executes, and assert the train
#      window count actually DROPPED relative to the _esmonly36 run -- the one
#      cheap check that proves --input took effect, because a cache collision
#      would otherwise hand back the old model silently (see train_esmonly36x.sh)
#   2  the remaining nine
#   3  assert this arm's channel_names equal _esmonly36's, which is what makes
#      reusing lf_pepend_residues_esmonly36.npz sound
#   4  ten full scans
#   5  tracks
#   6  the evaluations, including the covered-only ones for BOTH existing arms
#      and a byte-identity check that the KEEP parameterisation left the default
#      path alone
set -u
S=${0:A:h}
A2=~/AF2_analysis
HC=$A2/lf_headcmp
KEEPF=$S/zero_esm_covered.csv
PROCS=${PROCS:-2}
say() { echo "[$(date '+%F %H:%M:%S')] $*"; }

say "=== ESM-C-ONLY (COVERED PRECURSORS ONLY) ARM START ==="

# ---- stage 0: wait for the render ------------------------------------------
if [[ ${SKIP_WAIT:-0} != 1 ]]; then
  say "stage 0: waiting for final_pages.R to finish"
  for i in {1..300}; do                    # up to 5 h, checked each minute
    pgrep -f 'final_pages.R' > /dev/null || break
    [[ $((i % 15)) -eq 0 ]] && say "  still rendering ($(pgrep -f 'final_pages.R' | wc -l | tr -d ' ') workers, ${i} min)"
    sleep 60
  done
  if pgrep -f 'final_pages.R' > /dev/null; then
    say "stage 0: render STILL running after 5 h -- proceeding anyway at PROCS=1"
    PROCS=1
  else
    say "stage 0: render done"
  fi
fi
say "stage 0: load$(uptime | sed 's/.*load averages//')"

# ---- stage 1: one job, and prove --input took effect ------------------------
one=$A2/lf_pepend_run_C_cv5_f1_nomask_flatten_esmonly36x.rds
if [[ -s $one ]]; then
  say "stage 1: $(basename $one) already present, skipping the proof job"
else
  say "stage 1: one job (C fold 1) to prove the path executes"
  $S/train_esmonly36x.sh one 1 C
  [[ -s $one ]] || { say "stage 1 FAILED: no rds. tail:"; tail -25 $S/trainesm36x_C_f1.log; exit 1; }
  say "stage 1: ok, $(du -h $one | cut -f1)"
fi
# the decisive check: 10_8 prints the train/val counts it actually fitted on.
# _esmonly36 fitted C f1 on train 469 / val 417; this arm must be 461 / 413.
nw=$(grep -h '^terminus C:' $S/trainesm36x_C_f1.log | tail -1)
say "stage 1: $nw"
if [[ "$nw" != *"train 461"* || "$nw" != *"val 413"* ]]; then
  say "stage 1 FAILED: expected 'train 461, val 413' (see nn_input_esmcov_report.csv)."
  say "  Got: $nw"
  say "  This means --input did NOT take effect, or a cache collision returned the"
  say "  _esmonly36 model. Not continuing -- the whole arm would be the old one."
  exit 1
fi
say "stage 1: train/val counts confirm the covered-only input took effect"

# ---- stage 2: the rest ------------------------------------------------------
say "stage 2: the remaining nine training jobs (PROCS=$PROCS)"
PROCS=$PROCS $S/train_esmonly36x.sh
for k in 1 2 3 4 5; do for T in C N; do
  [[ -s $A2/lf_pepend_run_${T}_cv5_f${k}_nomask_flatten_esmonly36x.rds ]] || { say "stage 2 FAILED: $T f$k missing"; exit 1; }
done; done
say "stage 2: all 10 trunks present"

# ---- stage 3: the residue table is reusable only if the channels match ------
say "stage 3: asserting channel_names match _esmonly36"
cmp_rc=$(LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript -e '
  a <- jsonlite::fromJSON("~/AF2_analysis/lf_pepend_run_C_cv5_f1_nomask_flatten_esmonly36_isolated/in/config.json")$channel_names
  b <- jsonlite::fromJSON("~/AF2_analysis/lf_pepend_run_C_cv5_f1_nomask_flatten_esmonly36x_isolated/in/config.json")$channel_names
  cat(if (identical(a, b)) sprintf("MATCH %d", length(b)) else "MISMATCH", "\n")' 2>/dev/null | tail -1)
say "stage 3: $cmp_rc"
[[ $cmp_rc == MATCH* ]] || { say "stage 3 FAILED: channel order differs -- do NOT reuse the residue npz; write residues_esmonly36x.py"; exit 1; }

# ---- stage 4: the scans -----------------------------------------------------
say "stage 4: ten full scans (PROCS=$PROCS)"
PROCS=$PROCS $S/fullscan_esmonly36x.sh
for k in 1 2 3 4 5; do for T in C N; do
  [[ -f $A2/lf_winxgb_fullscan_${T}_cv5_f${k}_nomask_flatten_esmonly36x.npz ]] || { say "stage 4 FAILED: scan $T f$k missing"; exit 1; }
done; done
say "stage 4: all 10 scans present"

# ---- stage 5: tracks --------------------------------------------------------
say "stage 5: tracks"
cd $S && LC_ALL=en_US.UTF-8 TAG=_nomask_flatten_esmonly36x /usr/local/bin/Rscript tracks_esmonly36.R
for T in C N; do
  [[ -s $A2/lf_pepend_cv5_tracks_${T}_flatten_esmonly36x.rds ]] || { say "stage 5 FAILED: tracks $T missing"; exit 1; }
done
say "stage 5: tracks built"

# ---- stage 6: the evaluations ----------------------------------------------
# Back up the curated TOL-4 files on record BEFORE the byte-identity run, which
# rewrites them in place. If they come back different, the originals are here.
mkdir -p $S/segroc_backup
cp -f $HC/segment_confusion_tol4.csv $HC/segment_roc_auc_tol4.csv $S/segroc_backup/ 2>/dev/null

run_eval() {  # label, ARM, KEEP
  local label=$1 arm=$2 keep=$3
  say "stage 6: $label"
  cd $S && LC_ALL=en_US.UTF-8 TOL=4 ARM=$arm KEEP=$keep \
    /usr/local/bin/Rscript segments_roc.R > $S/segroc_${label}.log 2>&1
  local rc=$?
  say "stage 6: $label rc=$rc"
}

run_eval esmonly36x        _esmonly36x ""
run_eval esmonly36x_cov    _esmonly36x $KEEPF
run_eval esmonly36_cov     _esmonly36  $KEEPF
run_eval curated_cov       ""          $KEEPF
run_eval curated_identity  ""          ""

say "stage 6: byte-identity check of the KEEP parameterisation's default path"
for f in segment_confusion_tol4.csv segment_roc_auc_tol4.csv; do
  a=$(md5 -q $S/segroc_backup/$f 2>/dev/null); b=$(md5 -q $HC/$f 2>/dev/null)
  if [[ "$a" == "$b" ]]; then say "  $f IDENTICAL ($b)"
  else say "  $f DIFFERS -- backup $a, new $b -- restoring the original"; cp -f $S/segroc_backup/$f $HC/$f; fi
done

say "=== ESM-C-ONLY (COVERED) ARM DONE ==="
ls -la $HC/segment_*_esmonly36x*_tol4*.csv $HC/segment_*_cov.csv 2>/dev/null
