## PR WITHIN each C-terminal motif class, on the 37 unseen peptides.
## Within-group: positives and negatives both restricted to anchors of that
## class, so recall reaches 100% inside every group and the groups are
## comparable to each other.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus"); TOL <- 2L
source("R/pepend_windows.R")
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")

tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
kn <- fread(file.path(L,"knowns_near.csv")); ll <- fread(file.path(L,"ligand_list.csv"))
up <- fread(file.path(L,"uniprot_peptides.csv")) %>% filter(type %in% c("peptide","propeptide"))
oth <- bind_rows(kn %>% filter(term=="C") %>% transmute(accession, anchor=as.integer(anchor)),
                 ll %>% transmute(accession, anchor=as.integer(end)),
                 up %>% transmute(accession, anchor=as.integer(end))) %>%
  filter(!is.na(anchor)) %>% distinct()

fcO <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")       # old cache -> coverage flag
prec <- fcO$prec; seqs <- setNames(prec$seq, prec$accession)
nd <- rbindlist(lapply(seq_len(nrow(prec)), function(i) {
  m <- fcO$feats[[prec$accession[i]]]
  data.table(accession=prec$accession[i], anchor=seq_len(nrow(m)),
             nodata=as.integer(rowSums(abs(m[, ch4, drop=FALSE]))==0))}))
setkey(nd, accession, anchor)
## Share of the window's REAL residues with no data. A C window spans
## anchor-27 .. anchor+8, clipped to the mature range -- a rolling mean over a
## fixed 36 returns NA at the protein end, which silently deleted every
## `terminus` negative and made that group's AP meaningless.
np_ <- setNames(prec$n_prot, prec$accession); cp_ <- setNames(prec$c_prot, prec$accession)
nd[, miss := {
  cs <- c(0, cumsum(nodata)); a <- anchor
  lo <- pmax(a - 27L, np_[[.BY$accession]]); hi <- pmin(a + 8L, cp_[[.BY$accession]])
  lo <- pmax(lo, 1L); hi <- pmin(hi, .N)
  ok <- hi >= lo
  out <- rep(NA_real_, .N)
  out[ok] <- (cs[hi[ok] + 1L] - cs[lo[ok]]) / (hi[ok] - lo[ok] + 1L)
  out
}, by = accession]

res <- list()
for (arm in c("old","rebuilt")) {
  f <- sprintf("~/AF2_analysis/lf_pepend_scan_C_%s.rds", ifelse(arm=="old","near","allres"))
  sc <- as.data.table(readRDS(path.expand(f))$scan)[, .(accession, gene, anchor, score)]
  ## motif of every anchor, per precursor
  sc[, motif := {
    s <- seqs[[.BY$accession]]
    if (is.null(s)) rep(NA_character_, .N) else lf_pepend_motif_vec(s, "C", anchor, 1L)
  }, by = accession]
  sc <- merge(sc, nd[, .(accession, anchor, miss)], by=c("accession","anchor"), all.x=TRUE)
  sc[, key := paste(accession, anchor)]
  sc[, is_pos := key %in% paste(tst$accession, tst$end)]
  near <- rbind(
    CJ(d=setdiff(-TOL:TOL,0L), i=seq_len(nrow(tst)))[, .(accession=tst$accession[i], anchor=tst$end[i]+d)],
    CJ(d=-TOL:TOL, i=seq_len(nrow(oth)))[, .(accession=oth$accession[i], anchor=oth$anchor[i]+d)])
  sc[, excl := key %in% paste(near$accession, near$anchor) & !is_pos]
  for (g in c("terminus","dibasic","G + dibasic","monobasic","other")) {
    for (pool in c("all","complete")) {
      d <- sc[excl==FALSE & motif==g]
      if (pool=="complete") d <- d[is_pos | (!is.na(miss) & miss==0)]
      np <- sum(d$is_pos); if (np==0) next
      o <- order(-d$score); y <- d$is_pos[o]
      tp <- cumsum(y); prec_ <- tp/seq_along(y)
      ap <- sum(prec_[y])/np
      nn <- sum(!y)
      ## PREVALENCE-FREE measures. AP and precision both scale with how common
      ## positives are in the group, and prevalence here spans 1000x between
      ## groups, so AP cannot be read across them. These can:
      ##   roc_auc  = P(a positive outranks a random negative)
      ##   pct_rank = where the positive sits as a FRACTION of its own group
      ##   rec@k%   = recall when you follow the top k% of that group
      rk <- sapply(d$score[d$is_pos], function(v) sum(d$score > v) + 1L)
      roc <- mean(sapply(rk, function(r) (nn - (r - 1)) / nn))   # = 1 - normalised rank over negatives
      pctr <- rk / nrow(d)
      recat <- function(f) mean(rk <= max(1, ceiling(f * nrow(d))))
      res[[length(res)+1]] <- tibble(arm, motif=g, pool, n_pos=np, n_anchors=nrow(d),
        prevalence=np/nrow(d), ap=ap, lift=ap/(np/nrow(d)),
        roc_auc=roc, pct_rank_median=median(pctr),
        rec_top0.1pct=recat(0.001), rec_top1pct=recat(0.01), rec_top10pct=recat(0.10),
        rank_median=median(rk), top10=sum(rk <= 10))
    }
  }
}
R <- bind_rows(res)
for (p in c("all","complete")) {
  cat(sprintf("\n=== within-group PR, negatives = %s ===\n",
      ifelse(p=="all","every anchor of that motif class","only complete-data anchors of that class")))
  print(R %>% filter(pool==p) %>%
    transmute(motif, arm, n_pos, n_anchors, prev=signif(prevalence,2),
              roc_auc=round(roc_auc,4), pct_rank=signif(pct_rank_median,3),
              `rec@0.1%`=round(rec_top0.1pct,2), `rec@1%`=round(rec_top1pct,2),
              `rec@10%`=round(rec_top10pct,2), ap=round(ap,4)) %>%
    arrange(motif, arm) %>% as.data.frame(), row.names=FALSE)
}
saveRDS(R, file.path(S2,"motif_pr_testset.rds")); write.csv(R, file.path(S2,"motif_pr_testset.csv"), row.names=FALSE)
