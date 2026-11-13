## Where do the 37 never-seen peptide windows sit in feature space, and are the
## ones that score badly unusual -- or ordinary windows the model simply misses?
## Each window is its 36 x 26 channel block flattened to 936 numbers. Training
## windows go in too, so "unusual" means unusual RELATIVE TO WHAT WAS TRAINED ON.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
source("R/pepend_windows.R")
fc  <- readRDS(file.path(A2,"lf_pepend_residue_cache.rds")); prec <- fc$prec
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
nn  <- readRDS(file.path(A2,"lf_pepend_nn_input_near.rds"))
ap  <- fc$all_params3

rows <- list()
for (tm in c("C","N")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  b[, key := paste(accession, anchor)]
  setkey(b, key)
  ## --- the 37 never-seen ends
  pa <- if (tm=="N") tst$start else tst$end
  for (i in seq_len(nrow(tst))) {
    a <- tst$accession[i]; j <- match(a, prec$accession); if (is.na(j)) next
    ws <- lf_pepend_start(tm, pa[i])
    x  <- lf_pepend_slice(fc$feats[[a]], ws, prec$n_prot[j], prec$c_prot[j])
    sc <- b[.(paste(a, pa[i]))]
    rows[[length(rows)+1]] <- list(set="never seen", term=tm, gene=tst$gene[i],
      label=tst$name[i], anchor=pa[i],
      score=if (nrow(sc)) sc$xgb_pm__all50[1] else NA_real_,
      folds=if (nrow(sc)) unlist(sc[, paste0("xgb_pm__f",1:5), with=FALSE]) else rep(NA_real_,5),
      v=as.numeric(as.matrix(x[, ap])))
  }
  ## --- the training/validation knowns, same construction
  k <- nn$knowns %>% filter(term == tm) %>% distinct(peps, .keep_all = TRUE)
  for (i in seq_len(nrow(k))) {
    a <- k$accession[i]; j <- match(a, prec$accession); if (is.na(j)) next
    ws <- lf_pepend_start(tm, k$anchor[i])
    x  <- lf_pepend_slice(fc$feats[[a]], ws, prec$n_prot[j], prec$c_prot[j])
    sc <- b[.(paste(a, k$anchor[i]))]
    rows[[length(rows)+1]] <- list(set=paste0("known (", k$split[i], ")"), term=tm,
      gene=k$gene[i], label=k$pep_id[i], anchor=k$anchor[i],
      score=if (nrow(sc)) sc$xgb_pm__all50[1] else NA_real_,
      folds=if (nrow(sc)) unlist(sc[, paste0("xgb_pm__f",1:5), with=FALSE]) else rep(NA_real_,5),
      v=as.numeric(as.matrix(x[, ap])))
  }
}
M <- do.call(rbind, lapply(rows, `[[`, "v"))
meta <- tibble(set=sapply(rows,`[[`,"set"), term=sapply(rows,`[[`,"term"),
               gene=sapply(rows,`[[`,"gene"), label=sapply(rows,`[[`,"label"),
               anchor=sapply(rows,`[[`,"anchor"), score=sapply(rows,`[[`,"score"))
fo <- do.call(rbind, lapply(rows, `[[`, "folds")); colnames(fo) <- paste0("f",1:5)
meta <- bind_cols(meta, as_tibble(fo))
cat(sprintf("windows: %d  (%s)\n", nrow(M), paste(names(table(meta$set)), table(meta$set), collapse=", ")))
keep <- apply(M, 2, sd) > 0
pc <- prcomp(M[, keep], center = TRUE, scale. = TRUE)
ve <- round(100*summary(pc)$importance[2,1:4], 1)
cat("variance explained PC1-4:", paste0(ve, "%", collapse=" "), "\n")
meta <- bind_cols(meta, as_tibble(pc$x[,1:4]))
## distance from each window to the centre of the TRAINING windows, in PC space
trc <- colMeans(pc$x[meta$set == "known (train)", 1:10, drop=FALSE])
meta$d_train <- sqrt(rowSums((sweep(pc$x[,1:10,drop=FALSE], 2, trc))^2))
saveRDS(list(meta=meta, ve=ve), file.path(S2,"cv5_pca.rds"))
ts <- meta %>% filter(set=="never seen", !is.na(score))
cat("\n=== never-seen ends that score low in EVERY fold (bottom third of the 37) ===\n")
ts <- ts %>% mutate(min_fold=pmin(f1,f2,f3,f4,f5), max_fold=pmax(f1,f2,f3,f4,f5))
print(ts %>% arrange(max_fold) %>% head(10) %>%
  transmute(term, gene, label, score=round(score,3), best_fold=round(max_fold,3),
            worst_fold=round(min_fold,3), PC1=round(PC1,1), dist_to_train=round(d_train,1)) %>%
  as.data.frame(), row.names=FALSE)
cat("\n=== is a low score explained by being far from the training windows? ===\n")
cat(sprintf("  Spearman(score, distance to training centre), never-seen: %.3f  (n=%d)\n",
            cor(ts$score, ts$d_train, method="spearman"), nrow(ts)))
kn <- meta %>% filter(grepl("^known", set), !is.na(score))
cat(sprintf("  same for the known ends: %.3f  (n=%d)\n", cor(kn$score, kn$d_train, method="spearman"), nrow(kn)))
cat(sprintf("  median distance to training centre: never-seen %.1f vs known %.1f\n",
            median(ts$d_train), median(kn$d_train)))
