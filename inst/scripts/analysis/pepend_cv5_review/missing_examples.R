suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
ch4 <- c("cons_rs_n","min_afm","mean_afm","relASA")
fc  <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")   # the SHIPPED cache
prec <- fc$prec
new <- readRDS(file.path(S2, "sa1_slim.rds"))   # secretome_aa_1 (Aug), all topologies
old <- readRDS(file.path(S2, "sa_slim.rds"))    # secretome_aa  (Jan), the one 10_7 read
rn <- split(seq_len(nrow(new)), new$accession); ro <- split(seq_len(nrow(old)), old$accession)

## how much of each named gene has NO data in the shipped cache
genes <- c("MUC2","LRP1","CELSR2","ANO8","CXCL14","NPY","ADM","MT-RNR2","VWF","FBN1")
summ <- bind_rows(lapply(genes, function(g) {
  i <- which(prec$gene == g)[1]; if (is.na(i)) return(NULL)
  a <- prec$accession[i]; m <- fc$feats[[a]]
  z <- rowSums(abs(m[, ch4, drop=FALSE])) == 0
  tibble(gene = g, accession = a, length = nrow(m),
         no_data = sum(z), pct = round(100*mean(z), 1),
         first_gap = if (any(z)) min(which(z)) else NA_integer_,
         old_tbl_rows = if (is.null(ro[[a]])) 0L else length(ro[[a]]),
         old_max_index = if (is.null(ro[[a]])) NA_integer_ else suppressWarnings(max(old$index[ro[[a]]], na.rm=TRUE)),
         new_tbl_rows = if (is.null(rn[[a]])) 0L else length(rn[[a]]))
}))
cat("=== how much of each gene has NO feature data in the SHIPPED cache ===\n")
print(as.data.frame(summ), row.names = FALSE)

## residue-level proof: cache says 0, the source table has a real value
detail <- bind_rows(lapply(genes, function(g) {
  i <- which(prec$gene == g)[1]; if (is.na(i)) return(NULL)
  a <- prec$accession[i]; m <- fc$feats[[a]]
  z <- which(rowSums(abs(m[, ch4, drop=FALSE])) == 0)
  if (!length(z)) return(NULL)
  s <- z[round(seq(1, length(z), length.out = min(3, length(z))))]
  n <- new[rn[[a]], ]
  tibble(gene = g, residue = s, AA = n$AA[s],
         cache_cons = m[s, "cons_rs_n"], cache_mean_afm = m[s, "mean_afm"],
         cache_relASA = m[s, "relASA"],
         src_cons = round(n$cons_rs_n[s], 4), src_mean_afm = round(n$mean_afm[s], 4),
         src_relASA = round(n$relASA[s], 4), src_topo = n$topo2[s],
         old_row_exists = !is.null(ro[[a]]) && s %in% old$index[ro[[a]]])
}))
cat("\n=== residue-level: SHIPPED cache value vs the SOURCE table's value ===\n")
cat("(src_* NA = genuinely absent upstream; a NUMBER next to cache 0 = data that exists but never reached the cache)\n")
print(as.data.frame(detail), row.names = FALSE)
write.csv(detail, file.path(S2, "missing_data_examples.csv"), row.names = FALSE)
write.csv(summ, file.path(S2, "missing_data_by_gene.csv"), row.names = FALSE)

## how often is a cache zero RECOVERABLE (source has a value) vs genuinely absent?
cat("\n=== over ALL 5,616 precursors: of residues the cache calls no-data, how many are recoverable? ===\n")
set.seed(1); samp <- sample(prec$accession, 800)
rec <- bind_rows(lapply(samp, function(a) {
  m <- fc$feats[[a]]; if (is.null(rn[[a]])) return(NULL)
  n <- new[rn[[a]], ]; if (nrow(n) != nrow(m)) return(NULL)
  z <- rowSums(abs(m[, ch4, drop=FALSE])) == 0
  if (!any(z)) return(NULL)
  tibble(n_zero = sum(z),
         cons_recoverable  = sum(z & !is.na(n$cons_rs_n)),
         afm_recoverable   = sum(z & !is.na(n$mean_afm)),
         asa_recoverable   = sum(z & !is.na(n$relASA)))
}))
t <- colSums(rec)
cat(sprintf("  residues the cache calls no-data: %s\n", format(t[["n_zero"]], big.mark=",")))
for (k in c("cons_recoverable","afm_recoverable","asa_recoverable"))
  cat(sprintf("  %-18s %s (%.1f%%)\n", sub("_recoverable","", k),
      format(t[[k]], big.mark=","), 100*t[[k]]/t[["n_zero"]]))
