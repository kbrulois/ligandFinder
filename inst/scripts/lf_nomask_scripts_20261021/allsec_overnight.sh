#!/bin/zsh
# The overnight run (Kevin, 2026-10-22: pages for all secretome proteins, in
# their current state). Chains the two stages and leaves a readable trail:
#
#   1. allsec_scan.sh   per-gene CNN + xgb window scans, 5 folds x 2 termini,
#                       then the all-50 pool            (~3.5 h, 4 at a time)
#   2. allsec_pages.sh  one HTML page per gene, shortest first, then the index
#
# Stage 2 cannot start before stage 1's pooled npz exists, so if a scan is
# already running this waits for it rather than starting a second one.
#
# BOTH STAGES ARE RESUMABLE. If this dies, is killed, or the machine sleeps,
# just run allsec_overnight.sh again: completed scans and existing pages are
# skipped and it continues from where it stopped. Nothing is recomputed.
#
#   nohup ./allsec_overnight.sh > allsec_overnight.log 2>&1 &
set -u
A2=~/AF2_analysis
S=${0:A:h}
L=$S/allsec_overnight.log
say() { echo "[$(date '+%F %H:%M:%S')] $*" }

say "=== OVERNIGHT RUN START ==="
say "disk free before: $(df -h ~ | tail -1 | awk '{print $4}')"

# ---- stage 1: scans -----------------------------------------------------
POOLED_C=$A2/lf_winxgb_scan_C_cv5all50_nomask_flatten_all.npz
POOLED_N=$A2/lf_winxgb_scan_N_cv5all50_nomask_flatten_all.npz
if [[ -f $POOLED_C && -f $POOLED_N ]]; then
  say "stage 1: pooled scans already present, skipping"
elif pgrep -f "lf_winxgb.scan_genes" > /dev/null; then
  # Wait on the running scan's OWN completion marker, not on process presence:
  # xargs leaves a gap between jobs, and starting a second allsec_scan.sh in
  # that gap would duplicate an in-flight job (its .npz does not exist yet, so
  # the skip-if-cached guard would not catch it).
  say "stage 1: a scan is already running -- waiting for its DONE marker"
  while ! grep -q "ALL-SECRETOME SCAN DONE" $S/allsec_scan_driver.log 2>/dev/null; do
    sleep 60
    pgrep -f "lf_winxgb.scan_genes" > /dev/null && continue
    # no scan process and no DONE marker: it died, so take over and resume
    grep -q "ALL-SECRETOME SCAN DONE" $S/allsec_scan_driver.log 2>/dev/null && break
    say "stage 1: the running scan stopped without finishing -- resuming it here"
    $S/allsec_scan.sh && break || { say "STAGE 1 FAILED -- rerun to resume"; exit 1; }
  done
  say "stage 1: scan finished"
else
  say "stage 1: starting scans"
  $S/allsec_scan.sh || { say "STAGE 1 FAILED -- rerun allsec_overnight.sh to resume"; exit 1; }
fi
[[ -f $POOLED_C && -f $POOLED_N ]] || { say "STAGE 1 INCOMPLETE: no pooled npz -- rerun to resume"; exit 1; }
say "stage 1 done. disk free: $(df -h ~ | tail -1 | awk '{print $4}')"

# ---- stage 2: pages -----------------------------------------------------
# Stop before the disk is actually full: pages average ~60 KB per residue, so
# the whole secretome is ~175 GB, and a run that fills the volume would leave
# truncated HTML behind rather than a clean partial gallery.
MIN_FREE_GB=60
say "stage 2: rendering (floor: stop if free disk drops under ${MIN_FREE_GB} GB)"
(
  while sleep 300; do
    free=$(df -g ~ | tail -1 | awk '{print $4}')
    if [[ $free -lt $MIN_FREE_GB ]]; then
      say "DISK FLOOR HIT (${free} GB free) -- stopping the render; rerun to resume"
      pkill -f "final_pages.R" ; pkill -f "allsec_pages.sh chunk"
      exit 0
    fi
  done
) & GUARD=$!
PROCS=${PROCS:-4} CHUNK=${CHUNK:-25} $S/allsec_pages.sh
kill $GUARD 2>/dev/null

say "disk free after: $(df -h ~ | tail -1 | awk '{print $4}')"
say "pages: $(ls $A2/pepend_pages_allsec/*.html 2>/dev/null | grep -vc index.html), $(du -sh $A2/pepend_pages_allsec 2>/dev/null | cut -f1)"
say "=== OVERNIGHT RUN END ==="
