suppressMessages({library(dplyr); library(purrr)})
A2 <- path.expand("~/AF2_analysis")
arms <- c(main_near = "lf_pepend_nn_input_near.rds",
          rebuilt   = "lf_pepend_nn_input_allres.rds",
          foldA     = "lf_pepend_nn_input_foldA.rds",
          foldB     = "lf_pepend_nn_input_foldB.rds")
arms <- arms[file.exists(file.path(A2, arms))]
K <- imap_dfr(arms, function(f, nm) {
  x <- readRDS(file.path(A2, f))
  x$knowns %>% distinct(peps, term, .keep_all = TRUE) %>%
    transmute(arm = nm, gene, pep_id, term, anchor, split, motif)
})
cat("=== windows per arm / split / terminus ===\n")
print(K %>% count(arm, term, split) %>% tidyr::pivot_wider(names_from = split, values_from = n) %>% as.data.frame(), row.names = FALSE)
cat("\n=== val_frac actually realised ===\n")
print(K %>% group_by(arm, term) %>% summarise(n = n(), val = sum(split=="val"),
      pct_val = round(100*mean(split=="val"),1), .groups="drop") %>% as.data.frame(), row.names = FALSE)

## has any end ever been VAL-only across every arm?
cat("\n=== across main + foldA + foldB: has each end ever been in a TRAIN fold? ===\n")
ev <- K %>% group_by(gene, pep_id, term) %>%
  summarise(arms = n(), ever_train = any(split=="train"), ever_val = any(split=="val"), .groups="drop")
print(ev %>% count(ever_train, ever_val) %>% as.data.frame(), row.names = FALSE)
cat(sprintf("\nends NEVER in any training fold: %d of %d\n", sum(!ev$ever_train), nrow(ev)))
if (any(!ev$ever_train)) print(as.data.frame(ev[!ev$ever_train, c("gene","pep_id","term","arms")]), row.names = FALSE)
saveRDS(K, file.path(commandArgs(TRUE)[1], "splits_K.rds"))
