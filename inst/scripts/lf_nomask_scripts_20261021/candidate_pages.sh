#!/bin/zsh
# Protein pages for a chosen set of precursors, drawn with a candidate set
# (Kevin, 2026-10-28), then the gallery index. Same final_pages.R and settings
# as the candidate gallery pepend_pages_candidates_arms3: the curated, ESM-C
# only and AAswap head tracks, the all-secretome per-gene scans (NPZ_SFX=_all),
# and the candidate set's end-segment and candidate tracks.
#
#   O=~/AF2_analysis/pepend_pages_cand_nocurated_mp0.6_nomp0.45
#   FILTER='curated_relation == "none" & mp_best_score > 0.6 & xgb_nomp_score > 0.45' \
#     OUT=$O Rscript candidate_select.R
#   GENES_FILE=$O/scan_genes.txt SFX=_cand_x ./arm_gene_scans.sh    # the arms' sd bands
#   OUT=$O ./candidate_pages.sh
#
# OUT holds what candidate_select.R wrote: genes.txt, selection.txt (the
# candidate set KEY, and the extra bundle LF_EXTRA_BUNDLE if the selection used
# one) and, for the index, focus.csv; arm_gene_scans.sh leaves arm_scan_sfx.txt
# there (ARM_SCAN_SFX). Each is taken from OUT unless given, and a given value
# that disagrees with OUT's is refused. RESUMABLE: a gene whose page is complete
# is skipped, so re-run the same command to continue after a crash.
#
# Each chunk renders into its own staging directory: a page is written twice
# (htmltools::save_html, then its JS inlined from dependency_files/), and
# save_html replaces dependency_files/ wholesale, so processes sharing one plot
# directory can pull it out from under each other -- the all-secretome run
# (shared directory, 5 processes) left IGSF10 and SPTLC1 still loading their
# JS from dependency_files/. A page reaches OUT only once it is complete and
# self-contained (page_ok), and only such a page counts as done.
#
# Rendering is MEMORY-bound -- every process carries plot_bundle.rds, the
# residue cache and TensorFlow (see allsec_giants.sh) -- so PROCS is low. The
# genes are dealt round-robin from the length-sorted list into chunks, so each
# chunk carries a share of the long precursors and no chunk trails the rest.
#   PROCS=3 CHUNK=10 ARMS=esmonly36,curated_AAswap
#   SKIP_INDEX=1   render only;  DRY=1  show the chunks only;  TITLE="..." the index title
#   FORCE=1        re-render every page (after a change to the page code). Until
#                  a forced pass has replaced them all, a page older than its
#                  start does not count as done, so a re-run resumes the pass
#   GENES=<file>   another gene list
# Chunk lists live in .render/lists; each chunk's log in .render/logs, appended
# run after run under a dated header.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
S=${0:A:h}
OUT=${OUT:?set OUT to the page directory (candidate_select.R writes its genes.txt)}
OUT=${~OUT}
GENES=${GENES:-$OUT/genes.txt}
# what OUT was built with: given values must agree with it
norm() { local x=${1/#\~/$HOME}; print -r -- ${x:A} }   # ~ and relative paths, resolved
from_out() {   # name, OUT's value, given value
  if [[ -n $3 && -n $2 && $3 != $2 && $(norm $3) != $(norm $2) ]]; then
    print -u2 -- "$1=$3 but $OUT was built with $2 -- refusing"; exit 1
  fi
  print -r -- ${3:-$2}
}
SEL=$OUT/selection.txt
KEY=$(from_out KEY "$([[ -s $SEL ]] && sed -n 's/^key  *//p' $SEL)" "${KEY:-}") || exit 1
KEY=${KEY:-xgb_pm_T0.5_tol4_arms3}
LF_EXTRA_BUNDLE=$(from_out LF_EXTRA_BUNDLE \
  "$([[ -s $SEL ]] && sed -n 's/^extra  *\(.*\) (LF_EXTRA_BUNDLE):.*/\1/p' $SEL)" "${LF_EXTRA_BUNDLE:-}") || exit 1
ARM_SCAN_SFX=$(from_out ARM_SCAN_SFX "$([[ -s $OUT/arm_scan_sfx.txt ]] && head -1 $OUT/arm_scan_sfx.txt)" \
  "${ARM_SCAN_SFX:-}") || exit 1
ARMS=${ARMS:-esmonly36,curated_AAswap}
PROCS=${PROCS:-3}
CHUNK=${CHUNK:-10}
RUN=$OUT/.render
export OUT GENES KEY ARMS RUN LF_EXTRA_BUNDLE ARM_SCAN_SFX
# the five precursors every other gallery renders at a reduced popup cap; each
# gets a chunk of its own so its cap does not spread to its chunk-mates
BIG_CAP=600
BIG=(ANO8 C3 C5 EGF COPA)
# genes the all-secretome scan (*_all.npz) skipped, with the suffix of their own
# scan files (own_scan_npz.py builds them); each gets a chunk of its own drawn
# with that NPZ_SFX. Without them a page has no curated xgb window head.
typeset -A OWN_NPZ
OWN_NPZ=(COL18A1 _col18a1)
# a forced pass in progress: pages older than this stamp are stale
STAMP=$RUN/.force_pass

# complete (ends in </html>) and self-contained (no JS still loaded from
# dependency_files/); a kill mid-write or a failed inlining fails it
page_ok() {
  [[ -s $1 && $(tail -c 16 $1) == *'</html>'* ]] || return 1
  ! head -c 2000000 $1 | grep -q 'src="dependency_files/'
}
# done = complete, and newer than an unfinished forced pass
done_ok() { page_ok $1 && { [[ ! -e $STAMP ]] || [[ $1 -nt $STAMP ]]; } }

if [[ ${1:-} == chunk ]]; then
  f=$2
  genes=($(cat $f))
  (( ${#genes} )) || exit 0
  cap=2000
  for g in $genes; do (( ${BIG[(Ie)$g]} )) && cap=$BIG_CAP; done
  sfx=_all
  [[ ${f:t} == own_* ]] && sfx=${OWN_NPZ[${${f:t}#own_}]}
  stage=$RUN/stage_${f:t}
  # a stage left by an interrupted run could hold a finished OLD page that the
  # check below would promote: start every chunk from an empty one
  [[ $stage == */.render/stage_* ]] && rm -rf -- $stage
  mkdir -p $stage $RUN/logs
  log=$RUN/logs/${f:t}.log
  print -r -- "### $(date '+%F %H:%M:%S') ${genes[*]}" >> $log
  MODE=cv5 NPZ_SFX=$sfx ARMS=$ARMS CAND_SETS=$KEY LF_PLOT_DIR=$stage LF_POPUP_MAX=$cap PROV_N=0 \
    LF_SRC=$W LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/final_pages.R $genes >> $log 2>&1
  rc=$?
  # a page script exits 0 having written nothing when make_protein_plot_win
  # fails inside its own tryCatch: judge the files, not the exit code
  ok=0; bad=()
  for g in $genes; do
    if page_ok $stage/$g.html; then mv -f $stage/$g.html $OUT/$g.html && (( ok++ ))
    else bad+=$g; fi
  done
  msg="  $(date '+%H:%M:%S') ${f:t}: $ok/${#genes} pages written, rc=$rc"
  if (( ${#bad} )); then msg+="; not written: ${bad[*]} (see $log)"; else rm -rf -- $stage; fi
  echo $msg
  exit 0
fi

[[ -s $GENES ]] || { echo "no gene list at $GENES -- run candidate_select.R first"; exit 1; }
mkdir -p $RUN/lists
# (N): a glob that matches nothing is empty, not zsh's "no matches found" error
rm -f $RUN/lists/*(N)
todo=()
while read g; do
  [[ -n $g ]] || continue
  [[ ${FORCE:-0} == 1 ]] || ! done_ok $OUT/$g.html && todo+=$g
done < $GENES
n=${#todo}
echo "CANDIDATE PAGES $(date '+%F %H:%M:%S')  $OUT"
echo "  set $KEY, arms $ARMS${ARM_SCAN_SFX:+ (arm scans $ARM_SCAN_SFX)}${LF_EXTRA_BUNDLE:+, extra bundle ${LF_EXTRA_BUNDLE:t}}"
echo "  $(grep -c . $GENES) genes, $n to render, $PROCS processes, chunks of ~$CHUNK"
[[ -z $ARM_SCAN_SFX ]] && echo "  (no ARM_SCAN_SFX: the ESM-C only and AAswap tracks get no sd bands -- arm_gene_scans.sh)"
if (( n )); then
  nch=$(( (n + CHUNK - 1) / CHUNK ))
  i=0
  for g in $todo; do
    if (( ${BIG[(Ie)$g]} )); then echo $g > $RUN/lists/big_$g; continue; fi
    if [[ -n ${OWN_NPZ[$g]:-} ]]; then
      if [[ -f $HOME/AF2_analysis/lf_winxgb_scan_N_cv5all50_nomask_flatten${OWN_NPZ[$g]}.npz ]]; then
        echo $g > $RUN/lists/own_$g; continue
      fi
      echo "  WARNING: no ${OWN_NPZ[$g]} scan files for $g -- its page will lack the curated xgb window head (own_scan_npz.py)"
    fi
    echo $g >> $RUN/lists/c_$(printf %03d $(( i % nch )))
    (( i++ ))
  done
  if [[ ${DRY:-0} == 1 ]]; then        # DRY=1: show the chunks, render nothing
    for f in $RUN/lists/*(N); do echo "  ${f:t}: $(tr '\n' ' ' < $f)"; done
    exit 0
  fi
  [[ ${FORCE:-0} == 1 ]] && : > $STAMP
  # NUL-separated: xargs would split an OUT holding a blank or a quote
  print -rN -- $RUN/lists/big_*(N) $RUN/lists/own_*(N) $RUN/lists/c_*(N) | xargs -0 -P $PROCS -n 1 ${0:A} chunk
fi

missing=()
while read g; do [[ -n $g ]] && ! done_ok $OUT/$g.html && missing+=$g; done < $GENES
echo "CANDIDATE PAGES DONE $(date '+%F %H:%M:%S')"
pages=($OUT/*.html(N)); pages=(${pages:#*/index.html})
echo "  ${#pages} pages, $(du -sh $OUT | cut -f1); missing or stale: ${#missing}"
if (( ${#missing} )); then
  echo "  not written: ${missing[*]} -- see $RUN/logs, then re-run (without FORCE) to retry just these"
else
  rm -f $STAMP
fi

if [[ ${SKIP_INDEX:-0} != 1 ]]; then
  sel=""
  [[ -s $SEL ]] && sel=$(sed -n 's/^filter  *//p' $SEL)
  fo=(); [[ -s $OUT/focus.csv ]] && fo=(--focus $OUT/focus.csv)
  LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/candidate_gallery_index.R --dir $OUT --key $KEY \
    $fo --selection "$sel" \
    --title "${TITLE:-LigandFinder — selected candidate peptides}" \
    || echo "INDEX FAILED (pages are fine; rerun candidate_gallery_index.R)"
fi
