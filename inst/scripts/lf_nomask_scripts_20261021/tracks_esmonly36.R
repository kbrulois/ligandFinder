## Assemble the per-anchor tracks rds for an arm, from its ten full scans
## (Kevin, 2026-10-24). This is metrics_final.R's own `flat_tracks()` lifted into
## a standalone, TAG-parameterised script so the ESM-C-only arm can be built
## WITHOUT editing metrics_final.R -- which PR #14 is also changing, and which
## hardcodes the curated arm's filename.
##
##   TAG=_nomask_flatten_esmonly36 Rscript tracks_esmonly36.R
##   TAG=_nomask_flatten           Rscript tracks_esmonly36.R   # reproduces the
##                                                              # curated arm
##
## Writes ~/AF2_analysis/lf_pepend_cv5_tracks_<T>_flatten<SUFFIX>.rds with
## columns accession, anchor, {cnn,xgb_pm}__f1..5 and {cnn,xgb_pm}__all50 --
## exactly the shape segments_roc.R reads, so the segmentation is identical
## apart from which scores it is given.
##
## The two heads come out of ONE scan: full_scan writes `cnn` (the flatten
## window head) and `xgb_pm` (the per-member xgboost head) into the same npz.
suppressMessages({library(data.table); library(reticulate)})
A2  <- path.expand("~/AF2_analysis")
TAG <- Sys.getenv("TAG", "_nomask_flatten_esmonly36")
## the tracks file keeps the established name shape: _flatten plus whatever the
## arm adds after it, so the curated arm round-trips to its own existing file
SUF <- sub("^_nomask_flatten", "", TAG)
np   <- import("numpy", convert = TRUE)
prec <- fread(file.path(A2, "pepend_locus", "prec_order.csv"))

for (tm in c("C", "N")) {
  p  <- file.path(A2, sprintf("lf_pepend_cv5_tracks_%s_flatten%s.rds", tm, SUF))
  if (file.exists(p)) { cat(sprintf("%s: %s already built\n", tm, basename(p))); next }
  fs <- sprintf("%s/lf_winxgb_fullscan_%s_cv5_f%d%s.npz", A2, tm, 1:5, TAG)
  if (!all(file.exists(fs))) {
    cat(sprintf("%s: SKIPPED, %d of 5 scans present\n", tm, sum(file.exists(fs)))); next
  }
  z1 <- np$load(fs[1], allow_pickle = TRUE)
  b  <- data.table(accession = prec$accession[z1[["prot_idx"]] + 1L],
                   anchor    = as.integer(z1[["anchor"]]))
  for (k in 1:5) {
    z <- np$load(fs[k], allow_pickle = TRUE)
    ## the folds must agree on the anchor order or the row means mix proteins
    stopifnot(identical(as.integer(z[["anchor"]]), b$anchor))
    b[[sprintf("cnn__f%d", k)]]    <- as.numeric(z[["cnn"]])
    b[[sprintf("xgb_pm__f%d", k)]] <- as.numeric(z[["xgb_pm"]])
  }
  b[, cnn__all50    := rowMeans(.SD), .SDcols = sprintf("cnn__f%d", 1:5)]
  b[, xgb_pm__all50 := rowMeans(.SD), .SDcols = sprintf("xgb_pm__f%d", 1:5)]
  saveRDS(b, p)
  cat(sprintf("%s: %s anchors -> %s\n", tm, format(nrow(b), big.mark = ","), basename(p)))
}
