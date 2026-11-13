suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")
ap <- fc$all_params3; scale_cols <- setdiff(ap, c("AA_hydro","AA_charge","AA_mw","AA_pI","padding"))
nd <- readRDS("~/AF2_analysis/nn_dat_cache.rds"); cr <- nd$chan_range; rm(nd); invisible(gc())
x <- readRDS(file.path(S2, "sa1_slim.rds"))
for (cn in scale_cols) x[[cn]] <- pmin(pmax((x[[cn]] - cr[1,cn])/(cr[2,cn]-cr[1,cn]), 0), 1)
rn <- split(seq_len(nrow(x)), x$accession)
old <- readRDS(file.path(S2, "sa_slim.rds")); ro <- split(seq_len(nrow(old)), old$accession)
probe <- c("relASA","mean_afm","Phi_cos")     # non-energy, non-rescaled
set.seed(1); samp <- sample(intersect(names(fc$feats), names(rn)), 300)
hits <- list()
for (a in samp) {
  o <- fc$feats[[a]]; n <- x[rn[[a]], ]
  if (nrow(o) != nrow(n)) next
  had <- rowSums(abs(o[, scale_cols, drop=FALSE])) > 0
  bad <- which(had & apply(abs(o[, probe, drop=FALSE]) -
          abs(as.matrix(replace(n[, probe], is.na(n[, probe]), 0))), 1, function(v) max(abs(v)) > 1e-6))
  if (length(bad)) hits[[a]] <- bad
}
cat(sprintf("precursors with any such residue: %d of %d; residues: %d\n",
            length(hits), length(samp), sum(lengths(hits))))
cat("\nper precursor (top 12):\n")
h <- sort(lengths(hits), decreasing = TRUE)
for (a in names(h)[1:min(12,length(h))]) {
  L <- nrow(fc$feats[[a]]); b <- hits[[a]]
  oi <- if (!is.null(ro[[a]])) old$index[ro[[a]]] else integer(0)
  cat(sprintf("  %-8s %-9s len %5d  n_bad %4d  bad range %5d-%5d  old table max index %s  dup old idx %d\n",
      a, fc$prec$gene[fc$prec$accession == a], L, length(b), min(b), max(b),
      if (length(oi)) max(oi, na.rm=TRUE) else "none", sum(duplicated(oi[!is.na(oi)]))))
}
