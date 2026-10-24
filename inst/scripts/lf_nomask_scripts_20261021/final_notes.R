## Per-gene card notes for the final (flatten, 5-fold) gallery: which of the
## gene's known ends trained the models, and which fold held each out. A page
## where every end is out-of-sample is not interchangeable with one whose ends
## were trained on (4 of 5 folds saw each). The index greens notes containing
## "all", so the out-of-sample wording carries it.
suppressMessages({library(data.table); library(tibble)})
A2 <- path.expand("~/AF2_analysis")
genes <- readLines(commandArgs(TRUE)[1])
kn <- rbindlist(lapply(c("N", "C"), function(t) {
  d <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", t)))
  ## a window can start before residue 1 (ANXA1_w-5-30), hence the -?
  d[, .(gene = sub("_w-?[0-9]+--?[0-9]+$", "", peps), term = t, anchor, fold)]
}))
## the never-seen peptides drawn on each page (the magenta line), one span each
u <- fread(file.path(dirname(commandArgs(TRUE)[1]), "cv5_unseen_ends.csv"))
u[, `:=`(start = ifelse(term == "N", anchor, anchor - len + 1L),
         end   = ifelse(term == "N", anchor + len - 1L, anchor))]
u <- unique(u[, .(gene, start, end)])
note <- vapply(genes, function(g) {
  k <- kn[gene == g][order(term, anchor)]
  ns <- u[gene == g][order(start)]
  ns_txt <- if (nrow(ns)) paste0("never-seen peptide(s): ", paste(sprintf("%d-%d", ns$start, ns$end), collapse = ", ")) else NULL
  tr_txt <- if (!nrow(k)) "no known end in any fold: all tracks out-of-sample" else
    paste0("training ends (4 of 5 folds saw each; held out by): ",
           paste(sprintf("%s %d (fold %d)", k$term, k$anchor, k$fold), collapse = ", "))
  paste(c(tr_txt, ns_txt), collapse = "; ")
}, character(1))
saveRDS(tibble(gene = genes, note = unname(note)), file.path(dirname(commandArgs(TRUE)[1]), "final_notes.rds"))
cat(sprintf("%d genes: %d with training ends, %d all out-of-sample\n",
            length(genes), sum(!grepl("^no known", note)), sum(grepl("^no known", note))))
