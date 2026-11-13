## Final model vs the one it replaces, 5-fold CV x 10 seeds, everything HELD OUT
## (each known end scored only by the fold whose val split holds it):
##   attention  lf_pepend_run_<T>_cv5_f<k>_slack8         (the current pages)
##   flatten    lf_pepend_run_<T>_cv5_f<k>_nomask_flatten (the new model)
## Ranking: 10_8's candidate rank. Shape: every anchor of the end's own precursor.
suppressMessages({library(data.table); library(yardstick); library(reticulate)})
A2 <- path.expand("~/AF2_analysis"); np <- import("numpy", convert = TRUE)
arms <- list(attention = list(sum = "lf_pepend_run_%s_cv5_f%d_slack8_summary.rds", npz = "attn"),
             flatten   = list(sum = "lf_pepend_run_%s_cv5_f%d_nomask_flatten_summary.rds", npz = "flatten"))
rows <- list(); prof <- list()
for (t in c("C", "N")) for (arm in names(arms)) {
  ends <- list(); pr <- c(); roc <- c()
  for (k in 1:5) {
    p <- file.path(A2, sprintf(arms[[arm]]$sum, t, k)); if (!file.exists(p)) next
    s <- readRDS(p)
    tr <- factor(s$val$truth, levels = c(0, 1))
    pr  <- c(pr,  pr_auc_vec(tr, s$val$score, event_level = "second"))
    roc <- c(roc, roc_auc_vec(tr, s$val$score, event_level = "second"))
    kk <- as.data.table(s$knowns_ranked)[split == "val"]
    kk[, anchor := if (t == "C") pep_end else pep_start]
    kk <- unique(kk, by = c("accession", "anchor"))
    zf <- file.path(A2, sprintf("lf_headcmp/cv5_anchors_%s_f%d_%s.npz", t, k, arms[[arm]]$npz))
    if (!file.exists(zf)) next
    z <- np$load(zf); fz <- z$files
    for (i in seq_len(nrow(kk))) {
      ka <- paste0(kk$accession[i], "__anchor"); if (!ka %in% fz) next
      an <- z[[ka]]; v <- z[[paste0(kk$accession[i], "__cnn")]]
      e <- match(kk$anchor[i], an); if (is.na(e)) next
      d <- an - kk$anchor[i]; loc <- abs(d) <= 10; ring <- abs(d) >= 3 & abs(d) <= 8
      ends[[length(ends) + 1]] <- data.table(
        rank_cand = kk$rank_cand[i], wrank = sum(v > v[e]) + 1L,
        near2 = abs(an[which.max(v)] - kk$anchor[i]) <= 2,
        loc_on_end = abs(d[loc][which.max(v[loc])]) <= 1,
        contrast = v[e] / mean(v[ring]), width = sum(loc & v >= v[e] / 2), end_score = v[e])
      prof[[length(prof) + 1]] <- data.table(term = t, arm = arm, d = d[loc], rel = v[loc] / max(v[loc]))
    }
  }
  if (!length(ends)) { message(t, " ", arm, ": nothing to score yet"); next }
  E <- rbindlist(ends)
  rows[[length(rows) + 1]] <- data.table(
    term = t, arm = arm, folds = length(pr), n = nrow(E),
    val_PR = sprintf("%.3f", mean(pr)), val_ROC = sprintf("%.3f", mean(roc)),
    cand_rank_med = median(E$rank_cand, na.rm = TRUE), top100 = sum(E$rank_cand <= 100, na.rm = TRUE),
    top500 = sum(E$rank_cand <= 500, na.rm = TRUE),
    gene_rank_med = median(E$wrank), gene_rank1 = sum(E$wrank == 1), gene_peak_2 = sum(E$near2),
    local_peak_on_end = sum(E$loc_on_end), contrast_med = round(median(E$contrast), 2),
    width_med = median(E$width), end_med = round(median(E$end_score), 3))
}
R <- rbindlist(rows)
options(width = 250); print(R, row.names = FALSE)
P <- rbindlist(prof)[, .(rel = mean(rel)), by = .(term, arm, d)]
cat("\nmean local profile around the held-out end (score / local max, d = anchor - end):\n")
print(dcast(P[d %in% c(-8, -4, -2, -1, 0, 1, 2, 4, 8)], term + arm ~ d, value.var = "rel")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE)
fwrite(R, file.path(A2, "lf_headcmp", "final_vs_attention_cv5.csv"))
