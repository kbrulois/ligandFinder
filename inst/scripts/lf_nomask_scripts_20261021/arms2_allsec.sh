#!/bin/zsh
# Every secretome page in the two-arm layout (Kevin, 2026-11-03: "pages like
# these for all 5K+ genes"), the drafts in pepend_pages_arms2_draft at scale:
# final_pages.R ARMS_ONLY=1 with the AAswap and ESM-C only arms, candidate set
# arms3 capped at 20 deep. pepend_pages_allsec, the curated gallery, is left as
# it is.
#
#   caffeinate -is nohup ./arms2_allsec.sh > ~/AF2_analysis/arms2_allsec.log 2>&1 &
#
# RESUMABLE, stage by stage: after a crash, a reboot or the disk floor, the same
# command carries on and recomputes nothing finished.
#   1 gene lists
#   2 arm_gene_scans.sh    each fold's sd, for the arm tracks' bands
#   3 candidate_pages.sh   the pages, PROCS at a time
#   4 candidate_pages.sh   the 19 giants, one at a time, once 3 has finished: the
#                          levers below keep them openable and each needs the
#                          machine to itself (ALLSEC_README.md, "The 19 giants")
#   5 the index, and a check of the render logs
#
# Rendering is memory-bound: a process carrying all 200 members peaks near 10 GB
# on this 32 GB machine, so PROCS is low. Starting a second run with another
# PROCS is safe only after the first has stopped. To change it mid-run,
#   touch $OUT/.render/STOP        # chunks in flight finish, no more start
# then re-run with the new PROCS. MIN_FREE_GB does the same on its own.
#
# Run this from a copy of the checkout the run can keep (see lf_paths.sh): R
# reads scripts as it goes, and every chunk sources R/ afresh, so editing a
# checkout under a running render kills it or mixes two page versions.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
S=${0:A:h}
A2=~/AF2_analysis
OUT=${OUT:-$A2/pepend_pages_allsec_arms2}
OUT=${~OUT}
SFX=${SFX:-_allsec_arms2}
ARMS=${ARMS:-curated_AAswap,esmonly36}
PROCS=${PROCS:-3}
CHUNK=${CHUNK:-25}
SCAN_PROCS=${SCAN_PROCS:-4}
MIN_FREE_GB=${MIN_FREE_GB:-60}
GIANT_MIN=3000; GIANT_CAP=500
say() { echo "[$(date '+%F %H:%M:%S')] $*" }
free_gb() { df -g ~ | tail -1 | awk '{print $4}' }

say "=== ARMS2 RUN START: $OUT, code $W ==="
say "$PROCS page processes, chunks of $CHUNK, $(free_gb) GB free, floor $MIN_FREE_GB GB"

# ---- 1 gene lists --------------------------------------------------------
# the renderable list, plus the two names it left out for colliding on
# <gene>.html: the bundle holds one accession of each, so they do not
L=${LISTS:-$A2/lf_nomask_scripts_20261021}
EXTRA=(COL18A1@P39060 OTOA@Q7RTW8)
mkdir -p $OUT/.render
rm -f $OUT/.render/STOP
cut -f1 $L/gene_lengths.txt | grep -vxFf $L/giant_genes.txt > $OUT/genes.txt
print -l ${EXTRA%@*} >> $OUT/genes.txt
cp $L/giant_genes.txt $OUT/giant_genes.txt
{ cut -f1 $L/gene_lengths.txt; print -l $EXTRA } > $OUT/scan_genes.txt
say "stage 1: $(wc -l < $OUT/genes.txt | tr -d ' ') pages + $(wc -l < $OUT/giant_genes.txt | tr -d ' ') giants"

# ---- 2 arm scans ---------------------------------------------------------
if [[ -s $OUT/arm_scan_sfx.txt ]]; then
  say "stage 2: scans already complete ($(cat $OUT/arm_scan_sfx.txt))"
else
  say "stage 2: arm scans, $SCAN_PROCS at a time"
  ARMS="${ARMS//,/ }" GENES_FILE=$OUT/scan_genes.txt SFX=$SFX PROCS=$SCAN_PROCS LOGDIR=$OUT/.render/armscan \
    $S/arm_gene_scans.sh || { say "STAGE 2 FAILED -- re-run to resume"; exit 1; }
fi

# ---- 3, 4 pages ----------------------------------------------------------
# The floor stops the queue, it kills nothing: a chunk in flight finishes, and
# what it adds is small next to the slack the floor leaves
( while sleep 300; do
    if (( $(free_gb) < MIN_FREE_GB )); then
      say "DISK FLOOR: $(free_gb) GB free -- no new chunks; free space and re-run to resume"
      touch $OUT/.render/STOP; exit 0
    fi
  done ) & GUARD=$!
trap 'kill $GUARD 2>/dev/null' EXIT

common=(OUT=$OUT ARMS_ONLY=1 ARMS=$ARMS SKIP_INDEX=1)
say "stage 3: pages"
env $common PROCS=$PROCS CHUNK=$CHUNK $S/candidate_pages.sh
if [[ -e $OUT/.render/STOP ]]; then say "STOPPED before stage 4 -- re-run to resume"; exit 0; fi
say "stage 4: giants, one at a time"
env $common GENES=$OUT/giant_genes.txt PROCS=1 CHUNK=1 CAP=$GIANT_CAP \
  LF_ENS_ONLY_MIN=$GIANT_MIN LF_NO_ALIGN_MIN=$GIANT_MIN $S/candidate_pages.sh
if [[ -e $OUT/.render/STOP ]]; then say "STOPPED before stage 5 -- re-run to resume"; exit 0; fi

# ---- 5 index, check ------------------------------------------------------
say "stage 5: index"
LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/candidate_gallery_index.R --dir $OUT \
  --key xgb_pm_T0.5_tol4_arms3 --arms-only --title "LigandFinder — all secretome, AAswap and ESM-C only" \
  || say "INDEX FAILED (pages are fine; rerun candidate_gallery_index.R)"

# A surface that fails drops its tracks and the page still counts as written, so
# read the logs: every page must have said it carried popups
logs=($OUT/.render/logs/*.log(N))
all=$(cat $OUT/genes.txt $OUT/giant_genes.txt | sort -u)
popup=$(grep -h '\] popup: ' $logs | sed 's/^ *\[\([^ ]*\) .*/\1/' | sort -u)
nopopup=($(comm -23 <(print -r -- $all) <(print -r -- $popup)))
pages=($OUT/*.html(N)); pages=(${pages:#*/index.html})
fails=$(grep -hcE "extra series '.*(failed|no PLM|lacks)|per-index .* failed|popup SKIPPED|pepend tracks SKIPPED" $logs | paste -sd+ - | bc)
say "check: ${#pages} pages, $(du -sh $OUT | cut -f1), $(free_gb) GB free"
say "check: no popup logged for ${#nopopup} gene(s): ${nopopup[1,30]}"
say "check: $fails surface or track failure line(s) in the logs"
say "=== ARMS2 RUN END ==="
