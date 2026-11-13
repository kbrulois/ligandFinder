## ROC / PR / offset for the 37 never-seen peptides, old arm vs rebuilt arm.
## Positives: the exact N (start) or C (end) anchor of each test peptide.
## EXCLUDED (neither class): anchors within +/-tol of a test end but not the end
## itself, and anchors within +/-tol of ANY other known end (docked knowns,
## ligand_list, UniProt peptide/propeptide) -- counting a real peptide end as a
## false positive would punish the model for being right. See the scan-PR note.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus")
TOL <- 2L

tst <- read.csv(file.path(S2, "test_set_tiered.csv"), stringsAsFactors = FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
cat(sprintf("test peptides: %d over %d genes\n", nrow(tst), n_distinct(tst$gene)))

## every OTHER known end, to hold out
kn <- fread(file.path(L, "knowns_near.csv")) %>% transmute(accession, anchor = as.integer(anchor), term)
ll <- fread(file.path(L, "ligand_list.csv")) %>% filter(!is.na(start), !is.na(end))
up <- fread(file.path(L, "uniprot_peptides.csv")) %>% filter(type %in% c("peptide","propeptide"))
other <- bind_rows(
  kn %>% transmute(accession, anchor, term),
  ll %>% transmute(accession, anchor = as.integer(start), term = "N"),
  ll %>% transmute(accession, anchor = as.integer(end),   term = "C"),
  up %>% transmute(accession, anchor = as.integer(start), term = "N"),
  up %>% transmute(accession, anchor = as.integer(end),   term = "C")) %>%
  filter(!is.na(anchor)) %>% distinct()
cat(sprintf("other known ends held out: %d N, %d C\n",
            sum(other$term=="N"), sum(other$term=="C")))

arms <- list(old = "lf_pepend_scan_%s_near.rds", rebuilt = "lf_pepend_scan_%s_allres.rds")
res <- list(); curves <- list(); off <- list(); wg <- list()
for (arm in names(arms)) for (tm in c("N","C")) {
  p <- path.expand(file.path("~/AF2_analysis", sprintf(arms[[arm]], tm)))
  if (!file.exists(p)) { cat("missing:", basename(p), "\n"); next }
  sc <- as.data.table(readRDS(p)$scan)[, .(accession, gene, anchor, score)]
  pos_tab <- tst %>% transmute(accession, anchor = if (tm=="N") start else end, gene, name)
  oth <- other %>% filter(term == tm)
  ## exact positives present in the scan
  sc[, key := paste(accession, anchor)]
  pkey <- paste(pos_tab$accession, pos_tab$anchor)
  sc[, is_pos := key %in% pkey]
  ## exclusion: within TOL of a test end (not the end itself) or of any other end
  near <- rbind(
    CJ(d = setdiff(-TOL:TOL, 0L), i = seq_len(nrow(pos_tab)))[, .(accession = pos_tab$accession[i], anchor = pos_tab$anchor[i] + d)],
    CJ(d = -TOL:TOL, i = seq_len(nrow(oth)))[, .(accession = oth$accession[i], anchor = oth$anchor[i] + d)])
  sc[, excl := key %in% paste(near$accession, near$anchor) & !is_pos]
  d <- sc[excl == FALSE]
  np <- sum(d$is_pos); nn <- sum(!d$is_pos)
  ## ROC + PR from the ranking
  o <- order(-d$score); y <- d$is_pos[o]
  tp <- cumsum(y); fp <- cumsum(!y)
  prec <- tp/(tp+fp); rec <- tp/np
  ap <- sum(prec[y])/np                       # AP = mean precision at each positive
  fpr <- fp/nn
  auc <- sum(diff(c(0,fpr)) * (head(c(0,rec),-1) + c(0,rec)[-1])/2)
  res[[length(res)+1]] <- tibble(arm, term = tm, n_pos = np, n_neg = nn,
    prevalence = np/(np+nn), roc_auc = auc, pr_auc = ap, lift = ap/(np/(np+nn)))
  k <- unique(round(seq(1, length(y), length.out = 1500)))
  curves[[length(curves)+1]] <- tibble(arm, term = tm, rec = rec[k], prec = prec[k], fpr = fpr[k])
  ## genome-wide rank of each test end + distance from the gene's best anchor
  pr <- d[is_pos == TRUE][, .(accession, gene, anchor, score)]
  pr[, rank_global := sapply(score, function(v) sum(d$score > v) + 1L)]
  best <- sc[, .(best_anchor = anchor[which.max(score)], best_score = max(score)), by = accession]
  pr <- merge(pr, best, by = "accession", all.x = TRUE)
  pr[, offset := best_anchor - anchor]
  ## within-gene rank
  pr[, rank_in_gene := mapply(function(a, v) sum(sc$score[sc$accession == a] > v) + 1L, accession, score)]
  pr[, n_anchors := sapply(accession, function(a) sum(sc$accession == a))]
  off[[length(off)+1]] <- pr[, .(arm, term = tm, gene, accession, anchor, score,
                                 rank_global, rank_in_gene, n_anchors, offset)]
}
R <- bind_rows(res); C <- bind_rows(curves); O <- bind_rows(off)
cat("\n=== ROC / PR, per anchor, genome-wide ===\n")
print(R %>% mutate(prevalence = signif(prevalence,3), roc_auc = round(roc_auc,4),
                   pr_auc = round(pr_auc,4), lift = round(lift)) %>% as.data.frame(), row.names = FALSE)
cat("\n=== distance from the precursor's best-scoring anchor to the true end ===\n")
print(O %>% group_by(arm, term) %>% summarise(n = n(),
  at_0 = sum(offset == 0), within_2 = sum(abs(offset) <= 2), within_5 = sum(abs(offset) <= 5),
  median_abs = median(abs(offset)), .groups="drop") %>% as.data.frame(), row.names = FALSE)
cat("\n=== rank of the true end among its own precursor's anchors ===\n")
print(O %>% group_by(arm, term) %>% summarise(n = n(), median_rank = median(rank_in_gene),
  top1 = sum(rank_in_gene == 1), top5 = sum(rank_in_gene <= 5), top10 = sum(rank_in_gene <= 10),
  .groups="drop") %>% as.data.frame(), row.names = FALSE)
cat("\n=== genome-wide rank of the true end (of ~2.98M anchors) ===\n")
print(O %>% group_by(arm, term) %>% summarise(n = n(), median = median(rank_global),
  top1k = sum(rank_global <= 1000), top10k = sum(rank_global <= 10000),
  top100k = sum(rank_global <= 100000), .groups="drop") %>% as.data.frame(), row.names = FALSE)
saveRDS(list(R=R, C=C, O=O), file.path(S2, "testset_roc.rds"))
write.csv(O, file.path(S2, "testset_per_peptide.csv"), row.names = FALSE)
