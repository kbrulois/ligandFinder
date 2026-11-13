suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
fc  <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec
new <- readRDS(file.path(S2, "sa1_slim.rds")); old <- readRDS(file.path(S2, "sa_slim.rds"))
rn <- split(seq_len(nrow(new)), new$accession); ro <- split(seq_len(nrow(old)), old$accession)
detail <- bind_rows(lapply(c("MUC2","LRP1","CELSR2","ANO8","VWF","FBN1"), function(g) {
  i <- which(prec$gene == g)[1]; a <- prec$accession[i]; m <- fc$feats[[a]]
  z <- which(rowSums(abs(m[, ch4, drop=FALSE])) == 0)
  s <- z[round(seq(1, length(z), length.out = 3))]
  n <- new[rn[[a]], ]
  oi <- if (is.null(ro[[a]])) integer(0) else old$index[ro[[a]]]
  tibble(gene = g, res = s, AA = n$AA[s],
         cache_cons = m[s,"cons_rs_n"], cache_afm = m[s,"mean_afm"], cache_asa = m[s,"relASA"],
         SRC_cons = round(n$cons_rs_n[s],3), SRC_afm = round(n$mean_afm[s],3),
         SRC_asa = round(n$relASA[s],3), topo = n$topo2[s],
         jan_row = s %in% oi)
}))
print(as.data.frame(detail), row.names = FALSE)
write.csv(detail, file.path(S2, "missing_data_examples.csv"), row.names = FALSE)

cat("\n=== of residues the SHIPPED cache calls no-data, how many are recoverable from the source? ===\n")
set.seed(1); samp <- sample(prec$accession, 800)
rec <- bind_rows(lapply(samp, function(a) {
  m <- fc$feats[[a]]; if (is.null(rn[[a]])) return(NULL)
  n <- new[rn[[a]], ]; if (nrow(n) != nrow(m)) return(NULL)
  z <- rowSums(abs(m[, ch4, drop=FALSE])) == 0
  if (!any(z)) return(NULL)
  tibble(n_zero = sum(z), cons = sum(z & !is.na(n$cons_rs_n)),
         afm = sum(z & !is.na(n$mean_afm)), asa = sum(z & !is.na(n$relASA)))
}))
t <- colSums(rec)
cat(sprintf("  sampled 800 precursors; residues called no-data: %s\n", format(t[["n_zero"]], big.mark=",")))
for (k in c("cons","afm","asa")) cat(sprintf("    %s has a real value upstream: %s (%.1f%%)\n",
    k, format(t[[k]], big.mark=","), 100*t[[k]]/t[["n_zero"]]))
