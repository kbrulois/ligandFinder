## The never-seen sets, scored with all 50 models -- honest for them, because no
## fold ever trained on their genes. Same pool rule as the CV table: every other
## annotated end is removed rather than counted as a false positive.
suppressMessages({library(data.table); library(dplyr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
u  <- fread(file.path(S2, "cv5_unseen_ends.csv"))
ll <- as.data.frame(readRDS("inst/extdata/ligand_list.rds"))
fast_auc <- function(s, y) { r <- frank(s, ties.method="average"); n1 <- sum(y); n0 <- length(y)-n1
  (sum(r[y==1]) - n1*(n1+1)/2) / (n1*n0) }
fast_ap <- function(s, y) { o <- order(-s); yy <- y[o]
  sum((cumsum(yy)/seq_along(yy))[yy==1]) / sum(y) }
## AP re-derived at a COMMON prevalence. Recall is prevalence-free; precision is
## not, so the negatives are reweighted until the positive rate equals `target`
## and precision is recomputed along the same ranking. This makes AP comparable
## across sets whose prevalence differs 17-fold, which raw AP is not.
TARGET <- 0.01
fast_ap_at <- function(s, y, target = TARGET) {
  o <- order(-s); yy <- y[o]
  n1 <- sum(y); n0 <- length(y) - n1
  w  <- n1 * (1 - target) / (target * n0)      # weight on each negative
  tp <- cumsum(yy); fp <- cumsum(1 - yy)
  sum((tp / (tp + w * fp))[yy == 1]) / n1
}
SETS <- c("never seen (peptide, < 60 aa)" = "unseen_pep",
          "large cleavage product ends"   = "unseen_big")
res <- list()
for (tm in c("N","C")) {
  b <- as.data.table(readRDS(file.path(A2, sprintf("lf_pepend_cv5_tracks_%s.rds", tm))))
  kn <- fread(file.path(A2, sprintf("lf_pepend_cv5_knowns_%s.csv", tm)))
  b[, key := paste(accession, anchor)]; setkey(b, accession)
  ann <- ll %>% filter(!is.na(start), !is.na(end)) %>%
    transmute(accession, anchor = if (tm=="N") start else end) %>% distinct() %>%
    mutate(key = paste(accession, anchor))
  all_ann <- union(ann$key, paste(kn$accession, kn$anchor))
  for (nm in names(SETS)) {
    p <- u[set == nm & term == tm]
    if (!nrow(p)) next
    pk <- paste(p$accession, p$anchor)
    d <- b[!(key %in% setdiff(all_ann, pk))]
    y <- as.integer(d$key %in% pk)
    for (hd in c("cnn","xgb_pm")) {
      col <- paste0(hd, "__all50"); s <- d[[col]]; ok <- !is.na(s)
      pos <- d[y == 1, .(accession, s = get(col))]
      ranks <- vapply(seq_len(nrow(pos)), function(i) {
        g <- d[.(pos$accession[i])]; sum(g[[col]] > pos$s[i], na.rm=TRUE) + 1L }, 1L)
      res[[length(res)+1L]] <- data.table(set = SETS[[nm]], term = tm, head = hd,
        npos = sum(y[ok]), pool = sum(ok), auc = fast_auc(s[ok], y[ok]),
        ap = fast_ap(s[ok], y[ok]), ap1 = fast_ap_at(s[ok], y[ok]), med_rank = median(ranks), top5 = mean(ranks <= 5))
    }
  }
  rm(b); invisible(gc()); message(tm, " done")
}
R <- rbindlist(res); R[, prev := npos/pool]
fwrite(R, file.path(S2, "win_stats_unseen.csv"))
print(as.data.frame(R %>% mutate(across(c(auc,ap,top5), ~round(.x,4)), prev = signif(prev,3))),
      row.names = FALSE)
