## Same three sets, both heads: the CNN attention head and the per-member xgb
## head, from lf_winxgb_fullscan_{C,N}_near.npz (same 2,980,852 rows, so the two
## heads are directly comparable anchor for anchor).
suppressMessages({library(dplyr); library(data.table); library(reticulate)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus"); TOL <- 2L; PW <- 5L
np <- import("numpy", convert = TRUE)
prec <- fread(file.path(L, "prec_order.csv"))
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")
kn <- x$knowns %>% distinct(peps, term, .keep_all=TRUE)
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
ll <- fread(file.path(L,"ligand_list.csv")); up <- fread(file.path(L,"uniprot_peptides.csv")) %>%
  filter(type %in% c("peptide","propeptide"))

roc_res <- list(); off_res <- list()
for (tm in c("N","C")) {
  z <- np$load(path.expand(sprintf("~/AF2_analysis/lf_winxgb_fullscan_%s_near.npz", tm)), allow_pickle=TRUE)
  d <- data.table(accession = prec$accession[z[["prot_idx"]] + 1L],
                  anchor = as.integer(z[["anchor"]]),
                  cnn = as.numeric(z[["cnn"]]), xgb = as.numeric(z[["xgb_pm"]]))
  ## provenance: the npz's cnn column must match the scan rds we used before
  sc <- as.data.table(readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_scan_%s_near.rds", tm)))$scan)[, .(accession, anchor, score)]
  chk <- merge(d, sc, by=c("accession","anchor"))
  cat(sprintf("%s: %s rows; max |npz cnn - scan score| = %.3g\n", tm,
              format(nrow(chk), big.mark=","), max(abs(chk$cnn - chk$score))))
  setorder(d, accession, anchor)
  grp <- list(training = kn %>% filter(term==tm, split=="train") %>% transmute(accession, anchor),
              validation = kn %>% filter(term==tm, split=="val") %>% transmute(accession, anchor),
              `never seen` = tst %>% transmute(accession, anchor = if (tm=="N") start else end))
  allends <- bind_rows(bind_rows(grp),
    ll %>% transmute(accession, anchor=as.integer(if (tm=="N") start else end)),
    up %>% transmute(accession, anchor=as.integer(if (tm=="N") start else end))) %>%
    filter(!is.na(anchor)) %>% distinct()
  near <- CJ(k=-TOL:TOL, i=seq_len(nrow(allends)))[, .(key=paste(allends$accession[i], allends$anchor[i]+k))]
  d[, key := paste(accession, anchor)]
  d[, is_end := key %in% near$key]
  shifts <- c(lapply(seq_len(PW), \(k) list(k,"lag")), lapply(seq_len(PW), \(k) list(k,"lead")))
  for (hd in c("cnn","xgb")) {
    v <- d[[hd]]
    d[, rmax := do.call(pmax, c(list(get(hd)), lapply(shifts, \(z2) shift(get(hd), z2[[1]], fill=-Inf, type=z2[[2]])))), by=accession]
    pk_by <- split(d[get(hd) >= rmax]$anchor, d[get(hd) >= rmax]$accession)
    best  <- d[, .(top_anchor = anchor[which.max(get(hd))]), by=accession]
    negs  <- d[is_end==FALSE][[hd]]; nn <- length(negs)
    for (g in names(grp)) {
      q <- merge(as.data.table(grp[[g]]), d[, .(accession, anchor, s=get(hd))], by=c("accession","anchor"))
      if (!nrow(q)) next
      q[, rk := sapply(s, \(w) sum(negs > w) + 1L)]
      sall <- c(q$s, negs); y <- c(rep(TRUE, nrow(q)), rep(FALSE, nn))
      o <- order(-sall); y <- y[o]
      tp <- cumsum(y); fp <- cumsum(!y); prec_ <- tp/(tp+fp); rec <- tp/nrow(q); fpr <- fp/nn
      roc_res[[length(roc_res)+1]] <- tibble(term=tm, head=hd, set=g, n=nrow(q),
        roc_auc = sum(diff(c(0,fpr))*(head(c(0,rec),-1)+c(0,rec)[-1])/2),
        ap = sum(prec_[y])/nrow(q),
        pct_rank = median(q$rk/nn), rank_median = as.numeric(median(q$rk)),
        top10k = mean(q$rk <= 10000))
      q2 <- merge(q, best, by="accession", all.x=TRUE)
      q2[, top_in_prec := top_anchor - anchor]
      q2[, near_peak := mapply(function(a, an) { p <- pk_by[[a]]
        if (is.null(p) || !length(p)) NA_integer_ else p[which.min(abs(p-an))] - an }, accession, anchor)]
      off_res[[length(off_res)+1]] <- q2[, .(term=tm, head=hd, set=g, accession, anchor, top_in_prec, near_peak)]
    }
    d[, rmax := NULL]
  }
}
R <- bind_rows(roc_res) %>% mutate(set=factor(set, levels=c("training","validation","never seen")))
O <- rbindlist(off_res); O[, set := factor(set, levels=c("training","validation","never seen"))]
cat("\n=== retrieval, both heads ===\n")
print(R %>% transmute(term, head, set, n, roc_auc=round(roc_auc,4), ap=round(ap,4),
  pct_rank=sprintf("%.2f%%",100*pct_rank), rank_median, top10k=sprintf("%.0f%%",100*top10k)) %>%
  arrange(term, head, set) %>% as.data.frame(), row.names=FALSE)
f <- function(v) list(n=sum(!is.na(v)), exact=sum(v==0,na.rm=TRUE), w2=sum(abs(v)<=2,na.rm=TRUE), med=as.numeric(median(abs(v),na.rm=TRUE)))
for (meas in c("top_in_prec","near_peak")) {
  cat(sprintf("\n=== %s ===\n", meas))
  r <- O[, f(get(meas)), by=.(term, head, set)]
  r[, `:=`(`exact %`=sprintf("%.0f%%",100*exact/n), `w2 %`=sprintf("%.0f%%",100*w2/n))]
  print(r[order(term, head, set), .(term, head, set, n, exact, `exact %`, within2=w2, `w2 %`, median_abs=med)], row.names=FALSE)
}
saveRDS(list(R=R, O=O), file.path(S2,"xgb_redo.rds"))
