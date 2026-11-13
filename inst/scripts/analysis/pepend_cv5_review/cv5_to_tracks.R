## Put the CV tracks into the two shapes the protein page reads:
##   scan rds  $scan(accession, gene, anchor, score, sd)   -> CNN attention track
##   winxgb npz <gene>__anchor/__xgb_pm/__xgb_pm__sd       -> xgb window track
## all50 for the main tracks, and one npz per fold for the per-fold series.
suppressMessages({library(dplyr); library(data.table); library(reticulate)})
A2 <- path.expand("~/AF2_analysis"); np <- import("numpy", convert = FALSE)
GENES <- readLines(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "pages_all.txt"))
for (tm in c("C","N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  setkey(b, gene)
  ## the CNN track: all50 mean, sd = fold-to-fold spread
  scan_t <- b[, .(accession, gene, anchor, score = cnn__all50, sd = cnn__sd)]
  scan_t[, rank_all := frank(-score, ties.method = "first")]
  saveRDS(list(term = tm, n_members = 50L, scan = as_tibble(scan_t)),
          file.path(A2, sprintf("lf_pepend_scan_%s_cv5all50.rds", tm)))
  ## the xgb track, per gene, all50 and each fold
  ## Every series carries a real band. The per-fold ones take that fold's own
  ## spread over its 10 seeds (`xgb_pm__f<k>_sd`, from full_scan's xgb_pm_sd);
  ## writing zeros there, as a first pass did, draws a flat band and makes a
  ## single fold look far more certain than it is.
  mk <- function(col_mean, col_sd, path) {
    stopifnot(!is.null(col_sd), col_sd %in% names(b))
    args <- list()
    for (g in GENES) {
      d <- b[gene == g]
      if (!nrow(d)) { message("  ", g, ": not in the ", tm, " scan"); next }
      args[[paste0(g, "__anchor")]]     <- np$asarray(as.integer(d$anchor))
      args[[paste0(g, "__xgb_pm")]]     <- np$asarray(as.numeric(d[[col_mean]]))
      args[[paste0(g, "__xgb_pm__sd")]] <- np$asarray(as.numeric(d[[col_sd]]))
    }
    do.call(np$savez_compressed, c(list(path), args))
  }
  mk("xgb_pm__all50", "xgb_pm__sd", file.path(A2, sprintf("lf_winxgb_scan_%s_cv5all50.npz", tm)))
  for (k in 1:5)
    mk(sprintf("xgb_pm__f%d", k), sprintf("xgb_pm__f%d_sd", k),
       file.path(A2, sprintf("lf_winxgb_scan_%s_cv5_f%d.npz", tm, k)))
  message(sprintf("%s: scan rds + 6 winxgb npz written", tm))
  for (g in head(GENES, 0)) {
    d <- b[gene == g]
    if (!nrow(d)) next
    message(sprintf("   %-7s %4d anchors | all50 cnn peak %.3f @%d | xgb peak %.3f @%d",
      g, nrow(d), max(d$cnn__all50), d$anchor[which.max(d$cnn__all50)],
      max(d$xgb_pm__all50), d$anchor[which.max(d$xgb_pm__all50)]))
  }
}
