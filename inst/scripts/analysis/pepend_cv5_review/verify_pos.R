suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
x <- readRDS(file.path(S2, "sa1_slim.rds"))
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec; ap <- fc$all_params3
scale_cols <- setdiff(ap, c("AA_hydro","AA_charge","AA_mw","AA_pI","padding"))

## (1) is the data really there past row 1400 in the long proteins?
cat("--- data beyond row 1400, proteins whose old index capped at 1400 ---\n")
for (a in c("Q02817","Q07954","Q9HCU4")) {
  r <- x[x$accession == a, ]
  lo <- r[1:1400, ]; hi <- r[1401:nrow(r), ]
  cat(sprintf("%s %-7s rows %5d | relASA non-NA: rows1-1400 %.2f, rows1401+ %.2f | mean_afm %.2f / %.2f | cons_rs_n %.2f / %.2f\n",
    a, r$gene[1], nrow(r), mean(!is.na(lo$relASA)), mean(!is.na(hi$relASA)),
    mean(!is.na(lo$mean_afm)), mean(!is.na(hi$mean_afm)),
    mean(!is.na(lo$cons_rs_n)), mean(!is.na(hi$cons_rs_n))))
}

## (2) REGRESSION TEST: where the old index was valid, do positional features
## reproduce the old cache exactly? (the replacement for 10_7's 9.2 check)
old <- readRDS(file.path(S2, "sa_slim.rds"))     # the Jan table, extracellular-only
chk <- readRDS(file.path(S2, "sa1_chk.rds"))
good <- chk$acc[chk$idx_eq_pos %in% TRUE]
set.seed(1); samp <- sample(good, 400)
rows_new <- split(seq_len(nrow(x)), x$accession)
rows_old <- split(seq_len(nrow(old)), old$accession)
dev <- sapply(samp, function(a) {
  o <- old[rows_old[[a]], ]; o <- o[!is.na(o$index), ]
  if (!nrow(o)) return(NA_real_)
  n <- x[rows_new[[a]], ]
  if (!all(o$AA == n$AA[o$index])) return(Inf)        # letters must agree
  d <- as.matrix(o[, scale_cols]) - as.matrix(n[o$index, scale_cols])
  max(abs(d), na.rm = TRUE)
})
cat(sprintf("\n--- regression: %d precursors with a valid old index\n", sum(!is.na(dev))))
cat(sprintf("max |new_positional - old| over shared residues: %.3g  (Inf = letter disagreement: %d)\n",
            max(dev[is.finite(dev)], na.rm = TRUE), sum(!is.finite(dev) & !is.na(dev))))
