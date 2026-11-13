suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus"); TOL <- 2L
source("R/pepend_windows.R")
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")
kn <- x$knowns %>% distinct(peps, term, .keep_all = TRUE)
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")
seqs <- setNames(fc$prec$seq, fc$prec$accession)
npro <- setNames(fc$prec$n_prot, fc$prec$accession)
ll <- fread(file.path(L,"ligand_list.csv")); up <- fread(file.path(L,"uniprot_peptides.csv")) %>%
  filter(type %in% c("peptide","propeptide"))
mot <- function(acc, anch, tm) mapply(function(a, n) {
  s <- seqs[[a]]; if (is.null(s)) return(NA_character_)
  lf_pepend_motif(list(s), tm, list(n), list(npro[[a]]))
}, acc, anch)

res <- list()
for (tm in c("N","C")) {
  sc <- as.data.table(readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_scan_%s_near.rds", tm)))$scan)[, .(accession, anchor, score)]
  sc[, key := paste(accession, anchor)]
  grp <- list(
    training   = kn %>% filter(term==tm, split=="train") %>% transmute(accession, anchor),
    validation = kn %>% filter(term==tm, split=="val")   %>% transmute(accession, anchor),
    `never seen` = tst %>% transmute(accession, anchor = if (tm=="N") start else end))
  allends <- bind_rows(bind_rows(grp),
    ll %>% transmute(accession, anchor=as.integer(if (tm=="N") start else end)),
    up %>% transmute(accession, anchor=as.integer(if (tm=="N") start else end))) %>%
    filter(!is.na(anchor)) %>% distinct()
  near <- CJ(d=-TOL:TOL, i=seq_len(nrow(allends)))[, .(k=paste(allends$accession[i], allends$anchor[i]+d))]
  sc[, is_end := key %in% near$k]
  negs <- sc[is_end==FALSE]$score; nn <- length(negs); negs_sorted <- sort(negs, decreasing=TRUE)
  for (g in names(grp)) {
    d <- grp[[g]] %>% mutate(motif = mot(accession, anchor, tm))
    pos <- merge(as.data.table(d), sc[, .(accession, anchor, score)], by=c("accession","anchor"))
    if (!nrow(pos)) next
    pos[, rk := sapply(score, function(v) sum(negs > v) + 1L)]
    agg <- pos[, .(n=.N, roc=mean((nn-(rk-1))/nn), pct_rank=median(rk/nn),
                   rank_median=as.numeric(median(rk)), top10k=mean(rk<=10000)), by=motif]
    agg[, `:=`(term=tm, set=g)]
    res[[length(res)+1]] <- agg
  }
}
R <- rbindlist(res)
R[, set := factor(set, levels=c("training","validation","never seen"))]
for (tmm in c("C","N")) {
  cat(sprintf("\n=== %s terminus ===\n", tmm))
  print(R[term==tmm][order(motif, set)][, .(motif, set, n, roc=round(roc,4),
    pct_rank=sprintf("%.2f%%", 100*pct_rank), rank_median, `top10k`=sprintf("%.0f%%", 100*top10k))],
    row.names=FALSE)
}
saveRDS(R, file.path(S2,"gap_by_motif.rds")); write.csv(R, file.path(S2,"gap_by_motif.csv"), row.names=FALSE)
