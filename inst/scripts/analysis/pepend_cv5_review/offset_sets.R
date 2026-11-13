## Distance to the peptide end, for training / validation / never-seen.
## Two measures, because they answer different questions:
##   top_in_prec  offset of the precursor's single best-scoring anchor to this end
##                -- "if I looked at this protein, would my top call be the end?"
##   near_peak    signed distance to the NEAREST local maximum (the scan's peaks)
##                -- "is there a peak on the end at all, wherever it ranks?"
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]
x  <- readRDS("~/AF2_analysis/lf_pepend_nn_input_near.rds")
kn <- x$knowns %>% distinct(peps, term, .keep_all=TRUE)
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
out <- list()
for (tm in c("N","C")) {
  S  <- readRDS(path.expand(sprintf("~/AF2_analysis/lf_pepend_scan_%s_near.rds", tm)))
  sc <- as.data.table(S$scan)[, .(accession, anchor, score)]
  pk <- as.data.table(S$peaks)[, .(accession, anchor, score)]
  best <- sc[, .(top_anchor = anchor[which.max(score)]), by=accession]
  grp <- list(training = kn %>% filter(term==tm, split=="train") %>% transmute(accession, anchor),
              validation = kn %>% filter(term==tm, split=="val") %>% transmute(accession, anchor),
              `never seen` = tst %>% transmute(accession, anchor = if (tm=="N") start else end))
  for (g in names(grp)) {
    d <- as.data.table(grp[[g]]); d <- merge(d, best, by="accession", all.x=TRUE)
    d[, top_in_prec := top_anchor - anchor]
    ## pk[accession == accession[i]] would resolve BOTH sides to pk's own
    ## column -- split pk once and look the block up by name instead
    pk_by <- split(pk$anchor, pk$accession)
    d[, near_peak := mapply(function(a, an) {
        p <- pk_by[[a]]
        if (is.null(p) || !length(p)) NA_integer_ else p[which.min(abs(p - an))] - an
      }, accession, anchor)]
    out[[length(out)+1]] <- d[, .(term=tm, set=g, accession, anchor, top_in_prec, near_peak)]
  }
}
O <- rbindlist(out)
O[, set := factor(set, levels=c("training","validation","never seen"))]
f <- function(v) c(n=sum(!is.na(v)), at0=sum(v==0, na.rm=TRUE), w2=sum(abs(v)<=2, na.rm=TRUE),
                   w5=sum(abs(v)<=5, na.rm=TRUE), med=median(abs(v), na.rm=TRUE))
for (meas in c("top_in_prec","near_peak")) {
  cat(sprintf("\n=== %s ===\n", meas))
  r <- O[, as.list(f(get(meas))), by=.(term, set)]
  r[, `:=`(`exact %`=sprintf("%.0f%%",100*at0/n), `within2 %`=sprintf("%.0f%%",100*w2/n),
           `within5 %`=sprintf("%.0f%%",100*w5/n))]
  print(r[order(term, set), .(term, set, n, exact=at0, `exact %`, within2=w2, `within2 %`,
                              within5=w5, `within5 %`, median_abs=med)], row.names=FALSE)
}
saveRDS(O, file.path(S2,"offset_sets.rds"))
