#!/bin/zsh
# All 104 known-peptide pages with the ESM-C-only track, plus the gallery index
# (Kevin, 2026-10-24). The counterpart of allsec_pages.sh for the ESM arm.
#
#   ./esm_pages.sh              -> render everything still missing, then index
#   ./esm_pages.sh chunk <file> -> one chunk (what xargs calls)
#   PROCS=4 CHUNK=8 ./esm_pages.sh
#
# RESUMABLE: a gene whose .html already exists is skipped, so this survives an
# interruption or an OOM kill and picks up where it stopped. That matters here --
# this machine is memory-bound and the giant-precursor render may still be going.
#
# PROCS defaults to 4, as allsec_pages.sh uses. Drop it if something large is
# rendering alongside -- this machine is memory-bound under concurrent load.
#
# It renders through final_pages_esm.R, NOT final_pages.R: the ESM specs and the
# PLM-aware R/pepend_tracks.R live there, and final_pages.R is the record of
# what produced the curated galleries.
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
ARM=${ARM:-esmonly36}
OUT=${OUT:-$A2/pepend_pages_${ARM}}
## which page script, and which extra arms it draws. final_pages_arms.R
## takes any number of arms; final_pages_esm.R is the ESM-C-only record.
export PAGES_R=${PAGES_R:-final_pages_arms.R}
export ARMS=${ARMS:-esmonly36,curated_wAAid}
## the gallery title names the arms it actually drew, so a gallery cannot
## claim an arm set it does not have.
TITLE=${TITLE:-"LigandFinder — peptide-end pages, linear window head (5-fold CV): curated, ESM-C only and curated + AA ID"}
PROCS=${PROCS:-4}
CHUNK=${CHUNK:-6}
# the precursors the curated gallery rendered at a reduced popup cap
BIG_CAP=600
BIG=(ANO8 C3 C5 EGF COPA)

if [[ ${1:-} == chunk ]]; then
  f=$2
  genes=($(cat $f))
  (( ${#genes} )) || exit 0
  cap=2000
  for g in $genes; do (( ${BIG[(Ie)$g]} )) && cap=$BIG_CAP; done
  MODE=cv5 LF_PLOT_DIR=$OUT LF_POPUP_MAX=$cap ARMS=$ARMS \
    LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/$PAGES_R $genes \
    >> $S/esmpages_$(basename $f).log 2>&1
  echo "chunk $(basename $f): ${#genes} genes rc=$?"
  exit 0
fi

mkdir -p $OUT $S/esm_chunks
rm -f $S/esm_chunks/*
todo=$S/esm_chunks/todo.txt
: > $todo
while read -r g; do
  [[ -z $g ]] && continue
  [[ -s $OUT/$g.html ]] || echo $g >> $todo
done < $S/page_genes.txt
n=$(wc -l < $todo | tr -d ' ')
echo "ESM-C PAGES start $(date '+%F %H:%M:%S')"
echo "  $(wc -l < $S/page_genes.txt | tr -d ' ') known-peptide genes, $n to render, $PROCS at a time"
if [[ $n -eq 0 ]]; then
  echo "  nothing to render"
else
  split -l $CHUNK $todo $S/esm_chunks/c_
  ls $S/esm_chunks/c_* 2>/dev/null | xargs -P $PROCS -n 1 ${0:A} chunk
fi

echo "--- gallery index over whatever is on disk"
LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $W/inst/scripts/make_gallery_index.R \
  --dir $OUT \
  --title "$TITLE" \
  --blurb $S/final_blurb_tpfn.html --footer $S/final_footer.html --no-stock-desc \
  || echo "INDEX FAILED (pages are fine; rerun make_gallery_index.R)"

echo "ESM-C PAGES DONE $(date '+%F %H:%M:%S')"
echo "  rendered: $(ls $OUT/*.html 2>/dev/null | grep -vc index.html) pages, $(du -sh $OUT | cut -f1)"
missing=$(comm -23 <(sort $S/page_genes.txt) <(ls $OUT | sed 's/\.html$//' | sort) | grep -v '^index$' | wc -l | tr -d ' ')
echo "  still missing: $missing  (re-run esm_pages.sh to resume)"
