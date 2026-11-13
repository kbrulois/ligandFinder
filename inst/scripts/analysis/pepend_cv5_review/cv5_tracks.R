## Assemble the two tracks from the 5 folds x 2 termini x 10 seeds.
##
##   all50   mean over every fold's ensemble -- the track to SCAN a novel
##           protein with. For a KNOWN end it is partly in-sample: 4 of the 5
##           folds trained on that end, so never quote it as held out.
##   perfold one column per fold, and `cv` -- each known end taken from the ONE
##           fold that held it out, which is the honest cross-validated score.
suppressMessages({library(dplyr); library(data.table); library(reticulate)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
np <- import("numpy", convert = TRUE)
prec <- fread(file.path(A2, "pepend_locus", "prec_order.csv"))
for (tm in c("C","N")) {
  fs <- sprintf("%s/lf_winxgb_fullscan_%s_cv5_f%d.npz", A2, tm, 1:5)
  if (!all(file.exists(fs))) { message("missing scans for ", tm, " -- skipped"); next }
  z1 <- np$load(fs[1], allow_pickle = TRUE)
  base <- data.table(accession = prec$accession[z1[["prot_idx"]] + 1L],
                     gene = prec$gene[z1[["prot_idx"]] + 1L],
                     anchor = as.integer(z1[["anchor"]]))
  for (hd in c("cnn","xgb_pm")) {
    M <- sapply(fs, function(f) as.numeric(np$load(f, allow_pickle=TRUE)[[hd]]))
    stopifnot(nrow(M) == nrow(base))
    base[[paste0(hd, "__all50")]] <- rowMeans(M)
    base[[paste0(hd, "__sd_fold")]] <- apply(M, 1, sd)   # spread BETWEEN folds
    for (k in 1:5) base[[sprintf("%s__f%d", hd, k)]] <- M[, k]
    ## within-fold spread over that fold's 10 seeds. full_scan emits it for the
    ## xgb head only, so the CNN has no seed sd here and `__sd` stays the
    ## between-fold number for it -- do not read the two as the same quantity.
    if (hd == "xgb_pm") {
      S <- sapply(fs, function(f) as.numeric(np$load(f, allow_pickle=TRUE)[["xgb_pm_sd"]]))
      for (k in 1:5) base[[sprintf("xgb_pm__f%d_sd", k)]] <- S[, k]
      ## sd over all 50 models, by the law of total variance:
      ## E[within-fold variance] + variance of the fold means
      base[["xgb_pm__sd"]] <- sqrt(rowMeans(S^2) + apply(M, 1, var) * 4/5)
    } else {
      base[["cnn__sd"]] <- base[[paste0(hd, "__sd_fold")]]
    }
  }
  fwrite(base, file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.csv.gz", tm)))
  saveRDS(base, file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm)))
  message(sprintf("%s: %s anchors x (2 heads x [all50, sd, f1..f5]) -> lf_pepend_cv5_tracks_%s.rds",
                  tm, format(nrow(base), big.mark=","), tm))

  ## the honest CV score for each known end: the fold where it was validation
  kn <- rbindlist(lapply(1:5, function(k) {
    x <- readRDS(sprintf("%s/lf_pepend_nn_input_cv5_f%d.rds", A2, k))$knowns
    as.data.table(x)[term == tm & split == "val", .(accession, anchor, peps, gene, motif, fold = k)]
  }))
  kn <- unique(kn, by = c("accession","anchor"))
  kn <- merge(kn, base, by = c("accession","anchor"), all.x = TRUE)
  for (hd in c("cnn","xgb_pm")) {
    kn[[paste0(hd, "__cv")]] <- mapply(function(f, i) kn[[sprintf("%s__f%d", hd, f)]][i],
                                       kn$fold, seq_len(nrow(kn)))
  }
  negs <- base[!paste(accession, anchor) %in% paste(kn$accession, kn$anchor)]
  for (hd in c("cnn","xgb_pm")) {
    nv <- negs[[paste0(hd, "__all50")]]
    kn[[paste0(hd, "__cv_rank")]] <- sapply(kn[[paste0(hd, "__cv")]],
      function(v) if (is.na(v)) NA_integer_ else sum(nv > v) + 1L)
  }
  fwrite(kn, file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  ## an end whose anchor is outside the scanned mature range has no score --
  ## SAA1's N anchor is residue 18 and its scan starts at 19 (signal peptide)
  .ok <- !is.na(kn$cnn__cv_rank)
  message(sprintf("  %d held-out ends (%d scorable); median CV rank  cnn %s | xgb_pm %s",
                  nrow(kn), sum(.ok), format(median(kn$cnn__cv_rank, na.rm=TRUE)),
                  format(median(kn$xgb_pm__cv_rank, na.rm=TRUE))))
  if (any(!.ok)) message("    not scanned: ", paste(kn$peps[!.ok], collapse=", "))
}
