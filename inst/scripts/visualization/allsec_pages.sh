#!/bin/zsh
# Protein pages for the WHOLE secretome (Kevin, 2026-10-22), in exactly the
# state the 104-gene gallery is in: same final_pages.R, same options, same
# LF_POPUP_MAX default, the only change being which per-gene scan npz it reads
# (NPZ_SFX=_all, from allsec_scan.sh) and where pages land.
#
# Output goes to a NEW directory, so pepend_pages_flatten_cv5/ stays as the
# verified reference -- a gene rendered both ways must come out byte-identical
# apart from the build timestamp (see allsec_verify.sh).
#
# RESUMABLE, and that is the point: a gene whose .html already exists is
# skipped, so this can be interrupted, re-run, or restarted after a crash and
# it picks up where it stopped. Nothing is recomputed.
#
#   allsec_pages.sh              -> render everything still missing, then index
#   allsec_pages.sh chunk <file> -> one chunk (what xargs calls)
#   PROCS=4 CHUNK=25 allsec_pages.sh
set -u
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
S=${0:A:h}
OUT=$A2/pepend_pages_allsec
PROCS=${PROCS:-4}
CHUNK=${CHUNK:-25}
# the five precursors the 104-gene gallery rendered at a reduced popup cap
BIG_CAP=600
BIG=(ANO8 C3 C5 EGF COPA)

if [[ ${1:-} == chunk ]]; then
  f=$2
  genes=($(cat $f))
  (( ${#genes} )) || exit 0
  cap=2000
  for g in $genes; do (( ${BIG[(Ie)$g]} )) && cap=$BIG_CAP; done
  MODE=cv5 NPZ_SFX=_all LF_PLOT_DIR=$OUT LF_POPUP_MAX=$cap \
    LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/final_pages.R $genes \
    >> $S/allpages_$(basename $f).log 2>&1
  echo "chunk $(basename $f): ${#genes} genes rc=$?"
  exit 0
fi

mkdir -p $OUT $S/allsec_chunks
rm -f $S/allsec_chunks/*

# only what is still missing -- this is what makes a re-run cheap
# shortest precursors first (gene_lengths.txt is sorted by length): the useful
# bulk lands early, and anything that does not finish is the giant pages, whose
# HTML is too big to open anyway
todo=$S/allsec_chunks/todo.txt
: > $todo
while IFS=$'\t' read -r g L; do
  [[ -z $g ]] && continue
  [[ -s $OUT/$g.html ]] || echo $g >> $todo
done < $S/gene_lengths.txt
n=$(wc -l < $todo | tr -d ' ')
echo "ALL-SECRETOME PAGES start $(date '+%F %H:%M:%S')"
echo "  $(wc -l < $S/renderable_genes.txt | tr -d ' ') genes total, $n still to render, $PROCS processes, chunks of $CHUNK"
if [[ $n -eq 0 ]]; then
  echo "  nothing to render"
else
  # the five reduced-cap precursors each get a chunk of their own, so one
  # oversized page cannot drag a whole chunk down to its cap
  for g in $BIG; do
    if grep -qx $g $todo; then echo $g > $S/allsec_chunks/big_$g; sed -i '' "/^$g\$/d" $todo; fi
  done
  split -l $CHUNK $todo $S/allsec_chunks/c_
  ls $S/allsec_chunks/big_* $S/allsec_chunks/c_* 2>/dev/null | xargs -P $PROCS -n 1 ${0:A} chunk
fi

# the index, over whatever is actually on disk
LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $W/inst/scripts/make_gallery_index.R \
  --dir $OUT \
  --title "LigandFinder — peptide-end pages, linear window head (5-fold CV), all secretome" \
  --blurb $S/final_blurb_tpfn.html --footer $S/final_footer.html --no-stock-desc \
  || echo "INDEX FAILED (pages are fine; rerun make_gallery_index.R)"

echo "ALL-SECRETOME PAGES DONE $(date '+%F %H:%M:%S')"
echo "  rendered: $(ls $OUT/*.html 2>/dev/null | grep -vc index.html) pages, $(du -sh $OUT | cut -f1)"
missing=$(comm -23 <(sort $S/renderable_genes.txt) <(ls $OUT | sed 's/\.html$//' | sort) | grep -v '^index$' | wc -l | tr -d ' ')
echo "  still missing: $missing  (re-run allsec_pages.sh to resume)"
