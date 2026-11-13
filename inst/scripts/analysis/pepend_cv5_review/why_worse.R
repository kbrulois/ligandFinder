## WHY is the rebuilt arm worse on the full anchor pool?
## Hypothesis: anchors that were all-zero in the OLD cache were suppressed by
## the old model for free. In the rebuilt cache they carry real features, so the
## rebuilt model scores them on their merits -- and some outrank the true ends.
## Test: of the negatives that outrank each test end, what share were
## old-cache-uncovered, and what topology are they in?
suppressMessages({library(dplyr); library(data.table)})
S2 <- commandArgs(TRUE)[1]; L <- path.expand("~/AF2_analysis/pepend_locus")
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
fcO <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fcO$prec
nd <- rbindlist(lapply(seq_len(nrow(prec)), function(i) {
  m <- fcO$feats[[prec$accession[i]]]
  data.table(accession=prec$accession[i], anchor=seq_len(nrow(m)),
             nodata=as.integer(rowSums(abs(m[, ch4, drop=FALSE]))==0))}))
topo <- readRDS(file.path(S2,"sa1_slim.rds"))[, c("accession","index","topo2")]
setDT(topo); setnames(topo, "index", "anchor")
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
sc <- list()
for (arm in c("old","rebuilt")) {
  f <- sprintf("~/AF2_analysis/lf_pepend_scan_C_%s.rds", ifelse(arm=="old","near","allres"))
  d <- as.data.table(readRDS(path.expand(f))$scan)[, .(accession, anchor, score)]
  setnames(d, "score", arm); sc[[arm]] <- d
}
d <- merge(sc$old, sc$rebuilt, by=c("accession","anchor"))
d <- merge(d, nd, by=c("accession","anchor"), all.x=TRUE)
d <- merge(d, topo, by=c("accession","anchor"), all.x=TRUE)
d[, is_pos := paste(accession, anchor) %in% paste(tst$accession, tst$end)]

cat("=== score given to anchors the OLD cache had no data for ===\n")
cat(sprintf("such anchors: %s of %s (%.1f%%)\n", format(sum(d$nodata==1, na.rm=TRUE), big.mark=","),
            format(nrow(d), big.mark=","), 100*mean(d$nodata==1, na.rm=TRUE)))
print(d[!is.na(nodata), .(mean_old=round(mean(old),4), mean_new=round(mean(rebuilt),4),
     median_old=round(median(old),4), median_new=round(median(rebuilt),4), n=.N),
   by=.(`no data in old cache`=nodata==1)])

cat("\n=== of the negatives that OUTRANK the median test end, what are they? ===\n")
for (arm in c("old","rebuilt")) {
  thr <- median(d[[arm]][d$is_pos])
  up <- d[get(arm) > thr & is_pos == FALSE]
  cat(sprintf("\n%s arm: %s negatives score above the median test end\n", arm, format(nrow(up), big.mark=",")))
  cat(sprintf("  of those, uncovered in the old cache: %.1f%%\n", 100*mean(up$nodata==1, na.rm=TRUE)))
  tt <- sort(table(up$topo2), decreasing=TRUE)
  cat("  topology: ", paste(sprintf("%s %.0f%%", names(tt), 100*tt/sum(tt)), collapse="  "), "\n")
}
cat("\n(topology key: e extracellular, i intracellular, t transmembrane, l lumenal, s signal peptide)\n")
