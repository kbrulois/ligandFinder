## Same test set, but negatives restricted to anchors whose 36-residue window had
## FULL feature data in the OLD cache -- the pool where neither arm can win on
## coverage. This is the apples-to-apples version.
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus")
TOL <- 2L; ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")   # the OLD cache
nd <- rbindlist(lapply(seq_len(nrow(fc$prec)), function(i) {
  m <- fc$feats[[fc$prec$accession[i]]]
  data.table(accession = fc$prec$accession[i], anchor = seq_len(nrow(m)),
             nodata = as.integer(rowSums(abs(m[, ch4, drop=FALSE])) == 0))
}))
setkey(nd, accession, anchor)
## share of each C window's residues (anchor-27 .. anchor+8) with no data
nd[, miss_C := frollmean(nodata, 36, align = "right", fill = NA), by = accession]
nd[, miss_C := shift(miss_C, -8, type = "lag"), by = accession]   # window ends anchor+8
nd[, miss_N := frollmean(nodata, 36, align = "left", fill = NA), by = accession]
nd[, miss_N := shift(miss_N, 7, type = "lag"), by = accession]    # window starts anchor-7
cat(sprintf("per-residue no-data: %.1f%%\n", 100*mean(nd$nodata)))

tst <- read.csv(file.path(S2, "test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
kn <- fread(file.path(L,"knowns_near.csv")); ll <- fread(file.path(L,"ligand_list.csv"))
up <- fread(file.path(L,"uniprot_peptides.csv")) %>% filter(type %in% c("peptide","propeptide"))
other <- bind_rows(kn %>% transmute(accession, anchor=as.integer(anchor), term),
  ll %>% transmute(accession, anchor=as.integer(start), term="N"),
  ll %>% transmute(accession, anchor=as.integer(end),   term="C"),
  up %>% transmute(accession, anchor=as.integer(start), term="N"),
  up %>% transmute(accession, anchor=as.integer(end),   term="C")) %>%
  filter(!is.na(anchor)) %>% distinct()

arms <- list(old="lf_pepend_scan_%s_near.rds", rebuilt="lf_pepend_scan_%s_allres.rds")
res <- list(); curves <- list()
for (arm in names(arms)) for (tm in c("N","C")) {
  sc <- as.data.table(readRDS(path.expand(file.path("~/AF2_analysis", sprintf(arms[[arm]], tm))))$scan)[, .(accession, gene, anchor, score)]
  sc <- merge(sc, nd[, .(accession, anchor, m = get(paste0("miss_", tm)))], by=c("accession","anchor"), all.x=TRUE)
  pos <- tst %>% transmute(accession, anchor = if (tm=="N") start else end)
  oth <- other %>% filter(term==tm)
  sc[, key := paste(accession, anchor)]
  sc[, is_pos := key %in% paste(pos$accession, pos$anchor)]
  near <- rbind(
    CJ(d=setdiff(-TOL:TOL,0L), i=seq_len(nrow(pos)))[, .(accession=pos$accession[i], anchor=pos$anchor[i]+d)],
    CJ(d=-TOL:TOL, i=seq_len(nrow(oth)))[, .(accession=oth$accession[i], anchor=oth$anchor[i]+d)])
  sc[, excl := key %in% paste(near$accession, near$anchor) & !is_pos]
  for (pool in c("all negatives","complete-data negatives only")) {
    d <- if (pool == "all negatives") sc[excl==FALSE] else sc[excl==FALSE & (is_pos | (!is.na(m) & m == 0))]
    np <- sum(d$is_pos); nn <- sum(!d$is_pos); if (!np) next
    o <- order(-d$score); y <- d$is_pos[o]
    tp <- cumsum(y); fp <- cumsum(!y); prec <- tp/(tp+fp); rec <- tp/np; fpr <- fp/nn
    ap <- sum(prec[y])/np
    auc <- sum(diff(c(0,fpr)) * (head(c(0,rec),-1)+c(0,rec)[-1])/2)
    res[[length(res)+1]] <- tibble(arm, term=tm, pool, n_pos=np, n_neg=nn,
      prevalence=np/(np+nn), roc_auc=auc, pr_auc=ap, lift=ap/(np/(np+nn)))
    k <- unique(round(seq(1, length(y), length.out=1200)))
    curves[[length(curves)+1]] <- tibble(arm, term=tm, pool, rec=rec[k], prec=prec[k], fpr=fpr[k])
  }
}
R <- bind_rows(res); C <- bind_rows(curves)
cat("\n=== the 37 unseen peptides: both negative pools ===\n")
print(R %>% mutate(prevalence=signif(prevalence,3), roc_auc=round(roc_auc,4),
        pr_auc=round(pr_auc,4), lift=round(lift)) %>% as.data.frame(), row.names=FALSE)
saveRDS(list(R=R, C=C), file.path(S2, "testset_fair.rds"))
