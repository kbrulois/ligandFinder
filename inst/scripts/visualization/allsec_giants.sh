#!/bin/zsh
# The 19 longest secretome precursors (Kevin, 2026-10-24: "re-render the
# mucin-sized genes ... so the pages are openable").
#
# These were chunk c_hy of the overnight run. That chunk was EMPTIED while it
# was still pending and its gene list saved to giant_genes.txt, for two reasons:
# allsec_pages.sh would have rendered them at the default settings into pages
# of 0.2-1.0 GB, and being one chunk they would have run single-threaded on one
# worker for ~2 h while the other four idled.
#
# WHY NOT JUST A LOWER POPUP CAP, which is the obvious lever:
#   ADAMTS10 (1,103 aa), same code path, only LF_POPUP_MAX changed:
#       cap 2000 (not binding)  85.2 MB
#       cap   50                55.0 MB
#   i.e. 28.7 KB per suppressed anchor, but a 55 MB floor no cap can reach --
#   49.9 MB per 1,000 residues, which is ~720 MB for MUC16 on its own.
#
# WHAT THE PAGE ACTUALLY IS: counted in a loaded browser, CP (1,065 aa) has
# 315,878 SVG elements and 273,373 of them -- 87% -- are the species-alignment
# panel, which draws a <rect> and a <text> per species per residue over the 63
# species in species_dat. By bytes it is ~30% of a page, steady in every length
# band. At 14,507 residues a residue is well under a pixel wide, so none of it
# can be read anyway.
#
# THE RECIPE, three independent levers, all in final_pages_giant.R:
#   LF_NO_ALIGN_MIN   drop the species-alignment panel        (-87% of elements)
#   LF_ENS_ONLY_MIN   pooled all-50 ensembles only, no        (-43% measured on
#                     per-fold series: 2 series, not 12        ADAMTS10)
#   LF_POPUP_MAX      bound the popup payload                 (-28.7 KB/anchor)
#
# final_pages_giant.R is a COPY of final_pages.R. It has to be: R reads a
# script incrementally as it runs, so editing final_pages.R while the overnight
# render is using it kills that run at the end.
#
#   ./allsec_giants.sh              # render whatever is still missing
#   PROCS=1 ./allsec_giants.sh      # one at a time; safest, and MUC16 needs it
#   DRY=1 ./allsec_giants.sh        # just say what it would do
#
# RUN THIS ONLY AFTER THE MAIN RENDER HAS FINISHED. Page rendering is
# MEMORY-bound, not CPU-bound: every R process carries plot_bundle.rds, the
# residue cache and TensorFlow, so it costs GBs resident, and on this 32 GB
# machine the main run's five workers are already at the limit. Two giants
# beside them took swap to 76,810 MB of 78,848 (97%), left 124 MB of RAM free,
# put the machine at 3.45% user / 26.88% sys with 70% of cores IDLE -- every
# process blocked on paging, the OOM killer within 2 GB of the main run -- and
# halved its throughput from ~111,000 to ~57,000 residues/h. Killing the two
# freed 15 GB of swap and took user CPU straight back to 17.5%.
set -u
S=${0:A:h}
A2=~/AF2_analysis
OUT=$A2/pepend_pages_allsec
PROCS=${PROCS:-1}   # memory, not cores, is the limit -- see the note above
MIN=${MIN:-3000}                 # the length above which the levers kick in
CAP=${CAP:-500}

# GENES defaults to the whole list; point it elsewhere to render a subset
# (e.g. to leave out one gene that is already rendering in another process --
# a page only counts as done once its .html exists, so an in-flight gene would
# otherwise be started a second time and the two would fight over the file).
GENES=${GENES:-$S/giant_genes.txt}
[[ -s $GENES ]] || { echo "no gene list at $GENES"; exit 1; }

if [[ ${1:-} == one ]]; then
  g=$2
  MODE=cv5 NPZ_SFX=_all LF_PLOT_DIR=$OUT LF_POPUP_MAX=$CAP \
    LF_ENS_ONLY_MIN=$MIN LF_NO_ALIGN_MIN=$MIN PROV_N=0 \
    LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/final_pages_giant.R $g \
    >> $S/giant_$g.log 2>&1
  rc=$?
  sz=$([[ -s $OUT/$g.html ]] && du -h $OUT/$g.html | cut -f1 || echo "-")
  echo "$g rc=$rc $sz"
  exit $rc
fi

todo=()
while read g; do [[ -n $g && ! -s $OUT/$g.html ]] && todo+=$g; done < $GENES
echo "GIANTS  $(date '+%F %H:%M:%S')  ${#todo} of $(wc -l < $GENES | tr -d ' ') still to render"
echo "  recipe: no alignment + ensembles only above $MIN aa, popup cap $CAP, $PROCS processes"
(( ${#todo} )) || { echo "  nothing to do"; exit 0; }
# shortest first, so the cheap ones land while the monsters are still going
# DRY runs the SAME ordering pipeline and only swaps what it feeds, so a dry
# run actually tests the ordering rather than exiting before it
order() {
  for g in $todo; do
    printf "%s\t%s\n" "$(awk -F'\t' -v g=$g '$1==g{print $2}' $S/gene_lengths.txt)" "$g"
  done | sort -n
}
if [[ ${DRY:-0} == 1 ]]; then
  echo "  would render, shortest first:"
  order | awk -F'\t' '{printf "    %-10s %6s aa\n", $2, $1}'
  exit 0
fi
order | cut -f2 | xargs -P $PROCS -n 1 ${0:A} one

echo "GIANTS DONE $(date '+%F %H:%M:%S')"
printf "%-10s %8s %s\n" GENE AA PAGE
while read g; do
  printf "%-10s %8s %s\n" $g "$(awk -F'\t' -v g=$g '$1==g{print $2}' $S/gene_lengths.txt)" \
    "$([[ -s $OUT/$g.html ]] && du -h $OUT/$g.html | cut -f1 || echo MISSING)"
done < $S/giant_genes.txt
# The main run builds the index when IT finishes, over whatever is on disk at
# that moment. These pages may land afterwards, so rebuild it here -- same
# command and same arguments as allsec_pages.sh, or the gallery would differ.
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
if [[ ${SKIP_INDEX:-0} != 1 ]]; then
  echo
  echo "rebuilding the gallery index over $(ls $OUT/*.html | grep -vc index.html) pages"
  LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $W/inst/scripts/make_gallery_index.R \
    --dir $OUT \
    --title "LigandFinder — peptide-end pages, linear window head (5-fold CV), all secretome" \
    --blurb $S/final_blurb_tpfn.html --footer $S/final_footer.html --no-stock-desc \
    || echo "INDEX FAILED (pages are fine; rerun make_gallery_index.R)"
fi
