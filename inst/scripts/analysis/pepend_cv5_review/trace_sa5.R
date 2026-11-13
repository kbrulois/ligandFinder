suppressMessages({library(dplyr); library(purrr)})
S2 <- commandArgs(TRUE)[1]
## (1) how well does "share of residues with no data" alone separate known from the other negatives?
res <- readRDS(file.path(S2, "window_f0.rds"))
auc <- function(pos, neg) { r <- rank(c(pos, neg)); (sum(r[seq_along(pos)]) - length(pos)*(length(pos)+1)/2) / (length(pos)*length(neg)) }
for (tm in c("N","C")) for (sp in c("train","val")) {
  d <- res %>% filter(term == tm, split == sp, !is.na(f0))
  p <- -d$f0[d$win_type == "known"]
  cat(sprintf("%s %-5s  AUC of (1 - missing share), known vs all negatives: %.3f | vs non-near negatives: %.3f\n",
      tm, sp, auc(p, -d$f0[d$win_type != "known"]), auc(p, -d$f0[!d$win_type %in% c("known","near_end")])))
}
## (2) the 1,400 cap: where does `index` come from?
sa <- readRDS(file.path(S2, "sa_slim.rds"))
na_acc <- unique(sa$accession[is.na(sa$index)])
mx <- sapply(na_acc, function(a) suppressWarnings(max(sa$index[sa$accession == a], na.rm = TRUE)))
cat("\nmax non-NA index among the 125 NA-index accessions:\n"); print(head(sort(table(mx), decreasing = TRUE), 8))
sec <- readRDS("~/AF2_analysis/secretome_latest.rds")
for (a in c("Q02817","Q07954","Q9HCU4","P01871")) {
  i <- match(a, sec$accession)
  for (src in c("dssp","af_missense","cons")) {
    d <- sec[[src]][[i]]
    if (is.null(d) || !is.data.frame(d)) { cat(sprintf("%s %s %-11s: %s\n", a, sec$gene[i], src, class(d)[1])); next }
    ix <- if ("index" %in% names(d)) d$index else NA
    cat(sprintf("%s %-5s %-11s rows %5d (seq %5d)  index non-NA %5d  max %s\n", a, sec$gene[i], src, nrow(d),
        nchar(sec$sequence_uni[i]), sum(!is.na(ix)), suppressWarnings(max(ix, na.rm = TRUE))))
  }
}
