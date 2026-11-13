#!/bin/zsh
## ---- the cross-validated gallery: held-out pages + fixed-fold pages, one index ----
## Gathers ~/AF2_analysis/cv_heldout/<GENE>.html and every
## ~/AF2_analysis/cv_fold_pages/fold<F>/<GENE>.html (as <GENE>_fold<F>.html) into
## ~/AF2_analysis/cv_gallery, then writes its index.html. Both page sets come from
## 10_9w_cv_fold_pages.R.
##
## Pages are HARD LINKS, not copies: ~500 MB of pages costs no extra disk, and
## deleting the gallery leaves the source directories intact. Re-run after adding
## pages; it relinks everything and rebuilds the index.
##
##   zsh inst/scripts/build_cv_gallery.sh
## ------------------------------------------------------------------------------
set -e
set -o pipefail
## R in a bare "C" locale drops every non-ASCII character from the index (the
## title's em dash, the notes' middle dot) without a warning
export LC_ALL=${LC_ALL:-en_US.UTF-8}
A2=~/AF2_analysis
G=$A2/cv_gallery
ROOT=${0:A:h:h:h}
mkdir -p $G

for f in $A2/cv_heldout/*.html; do ln -f "$f" $G/; done
for d in $A2/cv_fold_pages/fold*(N/); do
  F=${d:t}
  for f in $d/*.html; do ln -f "$f" $G/${${f:t}%.html}_$F.html; done
done
## every page set ships identical htmlwidgets dependencies; link one copy
(cd $A2/cv_heldout/dependency_files &&
   find . -type d -exec mkdir -p $G/dependency_files/{} \; &&
   find . -type f -exec ln -f {} $G/dependency_files/{} \;)

## card notes: the held-out set's own (how many of each gene's known ends the page
## really holds out), plus one per fixed-fold page
T=$(mktemp -d)
## (a file, not Rscript -e: -e re-reads backslashes, so "\\." arrives as "\.")
cat > $T/notes.R <<'R'
a <- commandArgs(TRUE)
n <- as.data.frame(readRDS(file.path(a[1], "cv_meta", "cv_notes.rds")))[, c("gene", "note")]
fp <- sub("\\.html$", "", list.files(a[2], pattern = "_fold[AB]\\.html$"))
if (length(fp)) n <- rbind(n, data.frame(gene = fp,
  note = sprintf("fold %s trunks \u00b7 all held out (not a known)", sub("^.*_fold", "", fp))))
saveRDS(n, file.path(a[3], "notes.rds"))
R
Rscript $T/notes.R $A2 $G $T

cat > $T/blurb.html <<'HTML'
Every page here is drawn by models trained on only half the knowns (two
complementary folds), with each terminus taken from the fold that held that gene
out &mdash; so the scores are what a model that never saw the peptide says. The
card note says how many of the gene&rsquo;s known ends are actually held out;
<span style="color:#C8690B">orange</span> means some were trained on. Pages named
<code>_foldA</code> / <code>_foldB</code> are genes that are not knowns in either
fold, drawn once per fold, both fully held out.<br><br>
HTML
cat > $T/footer.html <<'HTML'
Trunks: <code>lf_pepend_run_{N,C}_fold{A,B}</code>, 20 members each; xgb heads
<code>lf_winxgb_cv_*</code>; residue models <code>lf_resid_preds_fold{A,B}</code>
(their own split, chosen per gene so they too are held out).
HTML

Rscript $ROOT/inst/scripts/make_gallery_index.R --dir $G \
  --title "LigandFinder — cross-validated, held-out predictions" \
  --notes $T/notes.rds --blurb $T/blurb.html --footer $T/footer.html
rm -rf $T
