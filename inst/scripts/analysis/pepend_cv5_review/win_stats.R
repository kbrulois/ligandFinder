## Held-out stats for the WINDOW models, fold by fold.
##
## A known end is honest only on the fold that held it out, so each fold k is
## scored on its OWN positives with its OWN column f{k}. Other annotated ends
## are removed from the pool entirely: they are neither clean positives (this
## fold trained on them) nor honest negatives (they are real ends), and counting
## them as false positives understates the model.
suppressMessages({library(data.table); library(dplyr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
ll <- as.data.frame(readRDS("inst/extdata/ligand_list.rds"))
im <- readRDS("data/id_mapping.rds")
ll$len <- ll$end - ll$start + 1
fast_auc <- function(score, y) {           # rank-based, so 3M rows is instant
  r <- frank(score, ties.method = "average")
  n1 <- sum(y); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fast_ap <- function(s, y) { o <- order(-s); yy <- y[o]
  sum((cumsum(yy)/seq_along(yy))[yy == 1]) / sum(y) }
## AP re-derived at a COMMON prevalence. Recall is prevalence-free; precision is
## not, so the negatives are reweighted until the positive rate equals `target`
## and precision is recomputed along the same ranking. This makes AP comparable
## across sets whose prevalence differs 17-fold, which raw AP is not.
TARGET <- 0.01
fast_ap_at <- function(s, y, target = TARGET) {
  o <- order(-s); yy <- y[o]
  n1 <- sum(y); n0 <- length(y) - n1
  w  <- n1 * (1 - target) / (target * n0)      # weight on each negative
  tp <- cumsum(yy); fp <- cumsum(1 - yy)
  sum((tp / (tp + w * fp))[yy == 1]) / n1
}
res <- list()
for (tm in c("N", "C")) {
  b  <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  b[, key := paste(accession, anchor)]
  kn[, key := paste(accession, anchor)]
  ## every OTHER annotated end on the same terminus, from ligand_list
  ann <- ll %>% filter(!is.na(start), !is.na(end)) %>%
    transmute(accession, anchor = if (tm == "N") start else end) %>%
    distinct() %>% mutate(key = paste(accession, anchor))
  drop_other <- union(ann$key, kn$key)
  for (k in 1:5) {
    pos_k <- kn[fold == k]
    keep  <- !(b$key %in% setdiff(drop_other, pos_k$key))
    d <- b[keep]
    y <- as.integer(d$key %in% pos_k$key)
    for (hd in c("cnn", "xgb_pm")) {
      s <- d[[sprintf("%s__f%d", hd, k)]]
      ok <- !is.na(s)
      rk <- d[y == 1, .(key, s = get(sprintf("%s__f%d", hd, k)), accession)]
      ranks <- vapply(seq_len(nrow(rk)), function(i) {
        g <- d[accession == rk$accession[i]]
        sum(g[[sprintf("%s__f%d", hd, k)]] > rk$s[i], na.rm = TRUE) + 1L }, 1L)
      res[[length(res) + 1L]] <- data.table(
        term = tm, head = hd, fold = k, npos = sum(y[ok]), pool = sum(ok),
        auc = fast_auc(s[ok], y[ok]), ap = fast_ap(s[ok], y[ok]), ap1 = fast_ap_at(s[ok], y[ok]),
        med_rank = median(ranks), top5 = mean(ranks <= 5))
    }
  }
  rm(b); invisible(gc())
  message(tm, " done")
}
R <- rbindlist(res)
R[, prev := npos / pool]
fwrite(R, file.path(S2, "win_stats_folds.csv"))
agg <- R[, .(folds = .N, npos = sum(npos), pool = round(mean(pool)),
             auc = mean(auc), auc_sd = sd(auc), ap = mean(ap), ap_sd = sd(ap), ap1 = mean(ap1),
             prev = mean(prev), med_rank = median(med_rank), top5 = mean(top5)),
         by = .(term, head)]
fwrite(agg, file.path(S2, "win_stats.csv"))
print(as.data.frame(agg %>% mutate(across(c(auc, auc_sd, ap, ap_sd, top5), ~round(.x, 4)),
                                   prev = signif(prev, 3))), row.names = FALSE)
