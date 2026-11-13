## Distance from each true end to the nearest PEAK in the xgb track.
##
## A residue is a peak when its score is at least the max over +/-5 (PW), the
## same rule cv5_scorers.R used. Two scoring rules are compared:
##   all50    the mean of all 50 models at that residue  (the ensemble average)
##   maxfold  the MAX of the 5 per-fold means            (top fold per position)
suppressMessages({library(data.table); library(dplyr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis"); PW <- 5L
u  <- fread(file.path(S2, "cv5_unseen_ends.csv"))
out <- list()
for (tm in c("N","C")) {
  b  <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  setorder(b, accession, anchor)
  b[, xgb_pm__maxfold := pmax(xgb_pm__f1, xgb_pm__f2, xgb_pm__f3, xgb_pm__f4, xgb_pm__f5)]
  SETS <- list(
    "cv"         = kn[, .(accession, anchor)],
    "unseen_pep" = u[set == "never seen (peptide, < 60 aa)" & term == tm, .(accession, anchor)],
    "unseen_big" = u[set == "large cleavage product ends"   & term == tm, .(accession, anchor)])
  for (col in c("xgb_pm__all50", "xgb_pm__maxfold")) {
    sh <- c(lapply(1:PW, \(j) list(j, "lag")), lapply(1:PW, \(j) list(j, "lead")))
    b[, rmax := do.call(pmax, c(list(get(col)),
         lapply(sh, \(z) shift(get(col), z[[1]], fill = -Inf, type = z[[2]])))), by = accession]
    pk <- split(b[get(col) >= rmax]$anchor, b[get(col) >= rmax]$accession)
    npk <- vapply(pk, length, 1L)
    for (sk in names(SETS)) {
      pos <- SETS[[sk]]
      off <- mapply(function(a, an) { p <- pk[[a]]
        if (is.null(p)) NA_integer_ else p[which.min(abs(p - an))] - an },
        pos$accession, pos$anchor)
      out[[length(out)+1L]] <- data.table(term = tm, set = sk,
        scorer = sub("^xgb_pm__", "", col), n = sum(!is.na(off)),
        exact = mean(off == 0, na.rm = TRUE), w1 = mean(abs(off) <= 1, na.rm = TRUE),
        w2 = mean(abs(off) <= 2, na.rm = TRUE), w5 = mean(abs(off) <= 5, na.rm = TRUE),
        med = median(abs(off), na.rm = TRUE),
        peaks_per_prot = median(npk[unique(pos$accession)], na.rm = TRUE))
    }
    b[, rmax := NULL]
  }
  rm(b); invisible(gc()); message(tm, " done")
}
R <- rbindlist(out); fwrite(R, file.path(S2, "peak_dist.csv"))
print(as.data.frame(R %>% mutate(across(c(exact,w1,w2,w5), ~round(.x,4)))), row.names = FALSE)
