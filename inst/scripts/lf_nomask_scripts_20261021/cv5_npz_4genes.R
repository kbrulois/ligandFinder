## The current model's xgb tracks (cv5: all-50 + one per fold) for the four page
## genes, in the npz shape the page reads. Same columns cv5_to_tracks.R wrote;
## CXCL17 was never in its gene list, so its pages had no xgb track. New file
## names -- the cv5 npz files the existing pages use are left alone.
suppressMessages({library(data.table); library(reticulate)})
A2 <- path.expand("~/AF2_analysis"); np <- import("numpy", convert = FALSE)
GENES <- c("CXCL14", "CXCL17", "ANO8", "NPPA")
for (tm in c("C", "N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  setkey(b, gene)
  mk <- function(col_mean, col_sd, path) {
    args <- list()
    for (g in GENES) {
      d <- b[gene == g]
      stopifnot(nrow(d) > 0)
      args[[paste0(g, "__anchor")]]     <- np$asarray(as.integer(d$anchor))
      args[[paste0(g, "__xgb_pm")]]     <- np$asarray(as.numeric(d[[col_mean]]))
      args[[paste0(g, "__xgb_pm__sd")]] <- np$asarray(as.numeric(d[[col_sd]]))
    }
    do.call(np$savez_compressed, c(list(path), args))
  }
  mk("xgb_pm__all50", "xgb_pm__sd",
     file.path(A2, sprintf("lf_winxgb_scan_%s_cv5all50_nomaskpages.npz", tm)))
  for (k in 1:5)
    mk(sprintf("xgb_pm__f%d", k), sprintf("xgb_pm__f%d_sd", k),
       file.path(A2, sprintf("lf_winxgb_scan_%s_cv5_f%d_nomaskpages.npz", tm, k)))
  ## the three genes the old files carry must come out identical
  old <- np$load(file.path(A2, sprintf("lf_winxgb_scan_%s_cv5all50.npz", tm)))
  new <- np$load(file.path(A2, sprintf("lf_winxgb_scan_%s_cv5all50_nomaskpages.npz", tm)))
  for (g in setdiff(GENES, "CXCL17")) {
    k <- paste0(g, "__xgb_pm")
    dd <- max(abs(py_to_r(old[[k]]) - py_to_r(new[[k]])))
    message(sprintf("%s %-6s max |old - new| = %g", tm, g, dd))
  }
}
