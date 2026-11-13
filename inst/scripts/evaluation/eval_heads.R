## Window head comparison, 2-fold CV x 5 seeds, no mask: attn vs flatten (x1, x10 lr).
## Every number is HELD OUT: each known end is scored by the one fold whose val
## split it sits in. Ranking: 10_8's candidate rank (against the motif-anchored
## candidate windows). Shape: every anchor of the end's own precursor
## (score_anchors.py), so peak sharpness can be read off directly.
suppressMessages({library(data.table); library(yardstick); library(reticulate)})
A2 <- path.expand("~/AF2_analysis"); np <- import("numpy", convert = TRUE)
arms <- c(attn = "", flatten = "_flatten", flatten_x10 = "_flatten_hlr10")
rows <- list(); prof <- list()
for (t in c("C", "N")) for (arm in names(arms)) {
  ends <- list(); pr <- c(); roc <- c()
  for (f in c("A", "B")) {
    p <- sprintf("%s/lf_pepend_run_%s_fold%s_nomask%s_summary.rds", A2, t, f, arms[[arm]])
    if (!file.exists(p)) next
    s <- readRDS(p)
    tr <- factor(s$val$truth, levels = c(0, 1))
    pr  <- c(pr,  pr_auc_vec(tr, s$val$score, event_level = "second"))
    roc <- c(roc, roc_auc_vec(tr, s$val$score, event_level = "second"))
    k <- as.data.table(s$knowns_ranked)[split == "val"]
    k[, anchor := if (t == "C") pep_end else pep_start]
    k <- unique(k, by = c("accession", "anchor"))
    z <- np$load(sprintf("%s/lf_headcmp/anchors_%s_fold%s_%s.npz", A2, t, f, arm))
    fz <- z$files
    for (i in seq_len(nrow(k))) {
      ka <- paste0(k$accession[i], "__anchor"); if (!ka %in% fz) next
      an <- z[[ka]]; v <- z[[paste0(k$accession[i], "__cnn")]]
      e <- match(k$anchor[i], an); if (is.na(e)) next          # end outside the mature range
      d <- an - k$anchor[i]; loc <- abs(d) <= 10; ring <- abs(d) >= 3 & abs(d) <= 8
      ends[[length(ends) + 1]] <- data.table(
        rank_cand = k$rank_cand[i],
        wrank = sum(v > v[e]) + 1L,
        near2 = abs(an[which.max(v)] - k$anchor[i]) <= 2,
        loc_on_end = abs(d[loc][which.max(v[loc])]) <= 1,
        contrast = v[e] / mean(v[ring]),
        width = sum(loc & v >= v[e] / 2),
        end_score = v[e], gene_max = max(v))
      pp <- v[loc] / max(v[loc]); prof[[length(prof) + 1]] <- data.table(term = t, arm = arm, d = d[loc], rel = pp)
    }
  }
  if (!length(ends)) next
  E <- rbindlist(ends)
  rows[[length(rows) + 1]] <- data.table(
    term = t, arm = arm, n = nrow(E),
    val_PR = sprintf("%.3f", mean(pr)), val_ROC = sprintf("%.3f", mean(roc)),
    cand_rank_med = median(E$rank_cand, na.rm = TRUE), top100 = sum(E$rank_cand <= 100, na.rm = TRUE),
    top500 = sum(E$rank_cand <= 500, na.rm = TRUE),
    gene_rank_med = median(E$wrank), gene_rank1 = sum(E$wrank == 1), gene_peak_2 = sum(E$near2),
    local_peak_on_end = sum(E$loc_on_end), contrast_med = round(median(E$contrast), 2),
    width_med = median(E$width), end_med = round(median(E$end_score), 3),
    gene_max_med = round(median(E$gene_max), 3))
}
R <- rbindlist(rows)
options(width = 250); print(R, row.names = FALSE)
P <- rbindlist(prof)[, .(rel = mean(rel)), by = .(term, arm, d)]
cat("\nmean local profile around the held-out end (score / local max, d = anchor - end):\n")
print(dcast(P[d %in% c(-8, -4, -2, -1, 0, 1, 2, 4, 8)], term + arm ~ d, value.var = "rel")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE)
fwrite(R, file.path(A2, "lf_headcmp", "head_comparison.csv"))
