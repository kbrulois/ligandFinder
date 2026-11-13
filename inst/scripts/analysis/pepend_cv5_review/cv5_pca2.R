## Redo, with the CV structure stated correctly.
##
## Under 5-fold CV there is no "train set" and "val set": every known end is
## TRAINING in 4 folds and VALIDATION in exactly 1. So a known gets two scores,
## and they are not interchangeable:
##   cv     the fold that held it out -- the honest one
##   all50  the mean of all 50 models, 40 of which trained on it -- in-sample
## The 37 never-seen peptides have only one kind of score, since no fold ever
## trained on them: for THEM all50 is already honest.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
source("R/pepend_windows.R")
fc <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); prec <- fc$prec
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
ap <- fc$all_params3
rows <- list()
for (tm in c("C","N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  b[, key := paste(accession, anchor)]; setkey(b, key)
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))   # carries fold + __cv
  pa <- if (tm=="N") tst$start else tst$end
  grab <- function(a, anc) { s <- b[.(paste(a, anc))]; if (!nrow(s)) NULL else s }
  for (i in seq_len(nrow(tst))) {
    a <- tst$accession[i]; j <- match(a, prec$accession); if (is.na(j)) next
    x <- lf_pepend_slice(fc$feats[[a]], lf_pepend_start(tm, pa[i]), prec$n_prot[j], prec$c_prot[j])
    s <- grab(a, pa[i])
    rows[[length(rows)+1]] <- list(set="never seen (no fold trained on it)", term=tm,
      label=tst$name[i], fold=NA_integer_,
      motif=lf_pepend_motif(list(prec$seq[j]), tm, list(pa[i]), list(prec$n_prot[j])),
      score_cv=if (is.null(s)) NA_real_ else s$xgb_pm__all50[1],   # honest for these
      score_all50=if (is.null(s)) NA_real_ else s$xgb_pm__all50[1],
      v=as.numeric(as.matrix(x[, ap])))
  }
  for (i in seq_len(nrow(kn))) {
    a <- kn$accession[i]; j <- match(a, prec$accession); if (is.na(j)) next
    x <- lf_pepend_slice(fc$feats[[a]], lf_pepend_start(tm, kn$anchor[i]), prec$n_prot[j], prec$c_prot[j])
    s <- grab(a, kn$anchor[i])
    rows[[length(rows)+1]] <- list(set="known (train in 4 folds, held out in 1)", term=tm,
      label=kn$peps[i], fold=kn$fold[i], motif=kn$motif[i],
      score_cv=kn$xgb_pm__cv[i],
      score_all50=if (is.null(s)) NA_real_ else s$xgb_pm__all50[1],
      v=as.numeric(as.matrix(x[, ap])))
  }
}
M <- do.call(rbind, lapply(rows,`[[`,"v"))
meta <- tibble(set=sapply(rows,`[[`,"set"), term=sapply(rows,`[[`,"term"),
               label=sapply(rows,`[[`,"label"), fold=sapply(rows,`[[`,"fold"),
               motif=sapply(rows,`[[`,"motif"),
               score_cv=sapply(rows,`[[`,"score_cv"), score_all50=sapply(rows,`[[`,"score_all50"))
keep <- apply(M,2,sd) > 0
pc <- prcomp(M[,keep], center=TRUE, scale.=TRUE)
ve <- round(100*summary(pc)$importance[2,1:4],1)
meta <- bind_cols(meta, as_tibble(pc$x[,1:4]))
cen <- colMeans(pc$x[, 1:10, drop=FALSE])
meta$d_cen <- sqrt(rowSums(sweep(pc$x[,1:10,drop=FALSE], 2, cen)^2))
saveRDS(list(meta=meta, ve=ve), file.path(S2,"cv5_pca2.rds"))
cat("\nend types present:\n"); print(table(meta$set, meta$motif))
cat(sprintf("windows %d | PC1-4 %s\n", nrow(M), paste0(ve,"%",collapse=" ")))
k <- meta %>% filter(grepl("^known", set), !is.na(score_cv), !is.na(score_all50))
cat("\n=== for the KNOWNS, what does training on them buy? ===\n")
print(k %>% group_by(term) %>% summarise(n=n(),
  `median held-out (cv)`=round(median(score_cv),3),
  `median all50 (in-sample)`=round(median(score_all50),3),
  `% where all50 > cv`=sprintf("%.0f%%", 100*mean(score_all50 > score_cv)), .groups="drop") %>%
  as.data.frame(), row.names=FALSE)
cat("\n=== score vs distance to the centroid of ALL peptide windows ===\n")
for (s in unique(meta$set)) {
  d <- meta %>% filter(set==s, !is.na(score_cv))
  cat(sprintf("  %-42s Spearman %.3f (n=%d)\n", s, cor(d$score_cv, d$d_cen, method="spearman"), nrow(d)))
}
