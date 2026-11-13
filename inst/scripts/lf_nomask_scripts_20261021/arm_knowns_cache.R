## Cache an INPUT_SFX window set's known ends with the fold that holds each one
## out, so a script that needs an arm's OWN known set (end_score_cards.R) reads
## a small csv instead of five 92 MB fold inputs. One row per known end: every
## known is `val` in exactly one of the five inputs.
## Writes <OUT_DIR>/lf_pepend_cv5_knowns_<T><SFX>.csv (OUT_DIR: ~/AF2_analysis).
##   INPUT_SFX=_msa LC_ALL=en_US.UTF-8 Rscript arm_knowns_cache.R
suppressMessages(library(data.table))
## Kevin, 2026-11-12: required, no default -- unset or "" would write over
## lf_pepend_cv5_knowns_<T>.csv, the shared end table every evaluator reads
SFX <- Sys.getenv("INPUT_SFX", "")
if (!nzchar(SFX)) stop("set INPUT_SFX, the window set to cache (e.g. INPUT_SFX=_msa)", call. = FALSE)
OUT_DIR <- path.expand(Sys.getenv("OUT_DIR", "~/AF2_analysis"))
for (tm in c("C", "N")) {
  d <- rbindlist(lapply(1:5, function(k) {
    x <- readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_nn_input_cv5_f%d%s.rds", k, SFX)))$knowns
    x <- x[x$term == tm & x$split == "val", ]
    data.table(accession = x$accession, anchor = as.integer(x$anchor), fold = k,
               peps = x$peps, gene = x$gene)
  }))
  stopifnot(!anyDuplicated(d[, .(accession, anchor)]))
  f <- file.path(OUT_DIR, sprintf("lf_pepend_cv5_knowns_%s%s.csv", tm, SFX))
  fwrite(d, f); cat(tm, nrow(d), "ends ->", f, "\n")
}
