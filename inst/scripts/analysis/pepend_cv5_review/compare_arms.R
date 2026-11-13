## old vs rebuilt trunks on the SAME held-out ends: val PR/ROC and retrieval,
## plus a paired test on log rank (the project's usual currency)
suppressMessages({library(dplyr); library(yardstick)})
A2 <- path.expand("~/AF2_analysis")
arms <- tibble::tribble(~label, ~file,
  "C old (near_slack8)",  "lf_pepend_run_C_near_slack8_summary.rds",
  "C rebuilt (allres)",   "lf_pepend_run_C_allres_slack8_summary.rds",
  "N old (near, NO slack)", "lf_pepend_run_N_near_summary.rds",
  "N rebuilt (allres)",   "lf_pepend_run_N_allres_slack8_summary.rds")
arms$path <- file.path(A2, arms$file)
arms <- arms %>% filter(file.exists(path))
res <- lapply(seq_len(nrow(arms)), function(i) {
  s <- readRDS(arms$path[i]); h <- s$knowns_ranked %>% filter(split == "val")
  tibble(label = arms$label[i], term = s$term, n_seeds = s$n_seeds,
         n_val = nrow(distinct(h, peps)),
         pr = pr_auc_vec(s$val$truth, s$val$score, event_level = "second"),
         roc = roc_auc_vec(s$val$truth, s$val$score, event_level = "second"),
         med = median(h$rank_cand), t100 = sum(h$rank_cand <= 100),
         t500 = sum(h$rank_cand <= 500), t1000 = sum(h$rank_cand <= 1000))
}) %>% bind_rows()
cat("\n=== validation and retrieval ===\n")
res %>% mutate(across(c(pr, roc), ~ round(.x, 3))) %>% as.data.frame() %>% print(row.names = FALSE)

for (tm in c("C","N")) {
  o <- arms$path[grepl(sprintf("^%s old", tm), arms$label)]
  n <- arms$path[arms$label == sprintf("%s rebuilt (allres)", tm)]
  if (!length(o) || !length(n)) next
  ## distinct(peps): summaries written before 28f9e57 double-count an end whose
  ## N and C windows share a `peps` name (EDN1/EDN2) -- 31 rows for 29 ends
  ho <- readRDS(o)$knowns_ranked %>% filter(split == "val") %>%
    distinct(peps, .keep_all = TRUE) %>% select(peps, motif, rank_old = rank_cand)
  hn <- readRDS(n)$knowns_ranked %>% filter(split == "val") %>%
    distinct(peps, .keep_all = TRUE) %>% select(peps, rank_new = rank_cand)
  j <- inner_join(ho, hn, by = "peps")
  if (!nrow(j)) next
  w <- suppressWarnings(wilcox.test(log(j$rank_old), log(j$rank_new), paired = TRUE))
  cat(sprintf("\n=== %s: %d paired held-out ends ===\n", tm, nrow(j)))
  cat(sprintf("  median rank %s -> %s;  better %d / worse %d / tied %d;  paired Wilcoxon on log rank p = %.3f\n",
      format(median(j$rank_old)), format(median(j$rank_new)),
      sum(j$rank_new < j$rank_old), sum(j$rank_new > j$rank_old),
      sum(j$rank_new == j$rank_old), w$p.value))
  cat("  by motif:\n")
  j %>% group_by(motif) %>% summarise(n = n(), old = median(rank_old), new = median(rank_new),
        better = sum(rank_new < rank_old), .groups = "drop") %>%
    arrange(new) %>% as.data.frame() %>% print(row.names = FALSE)
  cat("  biggest moves:\n")
  j %>% mutate(lr = log2(rank_new/rank_old)) %>% arrange(lr) %>%
    { bind_rows(head(., 4), tail(., 4)) } %>%
    transmute(peps, motif, rank_old, rank_new, log2_ratio = round(lr, 2)) %>%
    as.data.frame() %>% print(row.names = FALSE)
}
