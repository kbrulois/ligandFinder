suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
x <- readRDS(file.path(S2, "sa1_slim.rds")); old <- readRDS(file.path(S2, "sa_slim.rds"))
ch <- c("cons_rs_n","min_afm","mean_afm","relASA")
chk <- readRDS(file.path(S2, "sa1_chk.rds"))
good <- chk$acc[chk$idx_eq_pos %in% TRUE]
set.seed(1); samp <- sample(good, 500)
rn <- split(seq_len(nrow(x)), x$accession); ro <- split(seq_len(nrow(old)), old$accession)
dev <- sapply(samp, function(a) {
  if (is.null(ro[[a]])) return(NA_real_)
  o <- old[ro[[a]], ]; o <- o[!is.na(o$index), ]
  if (!nrow(o)) return(NA_real_)
  n <- x[rn[[a]], ]
  if (!all(o$AA == n$AA[o$index])) return(Inf)
  max(abs(as.matrix(o[, ch]) - as.matrix(n[o$index, ch])), na.rm = TRUE)
})
cat(sprintf("regression over %d precursors with a valid old index:\n", sum(!is.na(dev))))
cat(sprintf("  max |positional(new) - old| on the 4 channels: %.3g\n", max(dev[is.finite(dev)], na.rm=TRUE)))
cat(sprintf("  letter disagreements: %d\n", sum(!is.finite(dev) & !is.na(dev))))

## how much data does positional indexing RECOVER, by channel source?
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec
sub <- x[x$accession %in% prec$accession, ]
cat(sprintf("\ntotal cached residues: %d; rows now reachable: %d\n", sum(nchar(prec$seq)), nrow(sub)))
cat("non-NA share per source, over ALL residues of cached precursors:\n")
cat(sprintf("  cons_rs_n %.3f   min/mean_afm %.3f   relASA(DSSP) %.3f\n",
    mean(!is.na(sub$cons_rs_n)), mean(!is.na(sub$mean_afm)), mean(!is.na(sub$relASA))))
cat("old cache had data for 54.0% of residues (46.0% all-four-zero)\n")
