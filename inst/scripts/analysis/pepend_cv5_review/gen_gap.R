## Training vs validation held-out vs never-seen, on the CURRENT trunks.
## One model, one negative pool, three sets of positives -- so the three curves
## differ only in how much the model has seen of each.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus"); TOL <- 2L
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")
kn <- x$knowns %>% distinct(peps, term, .keep_all = TRUE)
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
ll <- fread(file.path(L,"ligand_list.csv")); up <- fread(file.path(L,"uniprot_peptides.csv")) %>%
  filter(type %in% c("peptide","propeptide"))

res <- list(); curves <- list()
for (tm in c("N","C")) {
  sc <- as.data.table(readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_scan_%s_near.rds", tm)))$scan)[, .(accession, gene, anchor, score)]
  sc[, key := paste(accession, anchor)]
  grp <- list(
    training   = kn %>% filter(term==tm, split=="train") %>% transmute(accession, anchor),
    validation = kn %>% filter(term==tm, split=="val")   %>% transmute(accession, anchor),
    `never seen` = tst %>% transmute(accession, anchor = if (tm=="N") start else end))
  ## one negative pool for all three: anchors not within TOL of ANY known end,
  ## from any group, plus ligand_list and UniProt peptides
  allends <- bind_rows(bind_rows(grp),
    ll %>% transmute(accession, anchor = as.integer(if (tm=="N") start else end)),
    up %>% transmute(accession, anchor = as.integer(if (tm=="N") start else end))) %>%
    filter(!is.na(anchor)) %>% distinct()
  near <- CJ(d=-TOL:TOL, i=seq_len(nrow(allends)))[, .(k = paste(allends$accession[i], allends$anchor[i]+d))]
  sc[, is_end := key %in% near$k]
  negs <- sc[is_end == FALSE]
  cat(sprintf("%s: %s negatives shared by all three curves\n", tm, format(nrow(negs), big.mark=",")))
  for (g in names(grp)) {
    pk <- paste(grp[[g]]$accession, grp[[g]]$anchor)
    pos <- sc[key %in% pk]
    if (!nrow(pos)) next
    d <- rbind(pos[, .(score, y=TRUE)], negs[, .(score, y=FALSE)])
    o <- order(-d$score); y <- d$y[o]; np <- sum(y); nn <- sum(!y)
    tp <- cumsum(y); fp <- cumsum(!y)
    prec <- tp/(tp+fp); rec <- tp/np; fpr <- fp/nn
    ap <- sum(prec[y])/np
    auc <- sum(diff(c(0,fpr)) * (head(c(0,rec),-1)+c(0,rec)[-1])/2)
    rk <- sapply(pos$score, function(v) sum(negs$score > v) + 1L)
    res[[length(res)+1]] <- tibble(term=tm, set=g, n=np, n_neg=nn, roc_auc=auc, ap=ap,
      pct_rank_median = median(rk/nn), rank_median = median(rk),
      top1k = mean(rk <= 1000), top10k = mean(rk <= 10000))
    k <- unique(round(seq(1, length(y), length.out=1500)))
    curves[[length(curves)+1]] <- tibble(term=tm, set=g, rec=rec[k], prec=prec[k], fpr=fpr[k])
  }
}
R <- bind_rows(res) %>% mutate(set = factor(set, levels=c("training","validation","never seen")))
C <- bind_rows(curves) %>% mutate(set = factor(set, levels=c("training","validation","never seen")))
cat("\n=== current trunks: seen vs held out vs never seen ===\n")
print(R %>% transmute(term, set, n, roc_auc=round(roc_auc,4), ap=round(ap,4),
  pct_rank=signif(pct_rank_median,3), rank_median, `in top 1k`=round(top1k,2),
  `in top 10k`=round(top10k,2)) %>% arrange(term, set) %>% as.data.frame(), row.names=FALSE)
saveRDS(list(R=R, C=C), file.path(S2, "gen_gap.rds"))
