## Cache each INPUT_SFX arm's (accession, anchor, fold) knowns, so arm_fix_knowns
## does not read five 92 MB fold inputs in every R process (a render spawns ~18).
## Writes ~/AF2_analysis/lf_pepend_cv5_knowns_<T><SFX>.csv.
suppressMessages(library(data.table))
SFX <- Sys.getenv("INPUT_SFX", "_padg")
for (tm in c("C", "N")) {
  d <- rbindlist(lapply(1:5, function(k) {
    x <- readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_nn_input_cv5_f%d%s.rds", k, SFX)))$knowns
    x <- x[x$term == tm & x$split == "val", ]
    data.table(accession = x$accession, anchor = as.integer(x$anchor), fold = k,
               peps = x$peps, gene = x$gene, pep_end_mature = x$pep_end_mature, g_amid = x$g_amid)
  }))
  stopifnot(!anyDuplicated(d[, .(accession, anchor)]))
  f <- path.expand(sprintf("~/AF2_analysis/lf_pepend_cv5_knowns_%s%s.csv", tm, SFX))
  fwrite(d, f); cat(tm, nrow(d), "ends ->", f, "\n")
}
