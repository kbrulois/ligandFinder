suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); prec <- fc$prec; ap <- fc$all_params3
scale_cols <- setdiff(ap, c("AA_hydro","AA_charge","AA_mw","AA_pI","padding"))
x <- readRDS("~/peptide_alg/build_residue_db/processed/secretome_aa_1.rds")
x <- as.data.frame(x)[, c("accession","gene","index","AA","topo2", scale_cols, "cons_rs","species_limit")]
saveRDS(x, file.path(S2, "sa1_slim.rds"))
rows <- split(seq_len(nrow(x)), x$accession)
chk <- bind_rows(lapply(prec$accession, function(a) {
  r <- rows[[a]]; s <- prec$seq[prec$accession == a]
  if (is.null(r)) return(data.frame(acc = a, n = nchar(s), nrow = 0L))
  aa <- paste(x$AA[r], collapse = ""); ix <- x$index[r]; pos <- seq_along(r)
  data.frame(acc = a, n = nchar(s), nrow = length(r), same_seq = identical(aa, s),
             n_na_idx = sum(is.na(ix)), idx_eq_pos = all(ix[!is.na(ix)] == pos[!is.na(ix)]))
}))
saveRDS(chk, file.path(S2, "sa1_chk.rds"))
cat("cached precursors:", nrow(chk), "\n")
cat("rows == sequence length:", sum(chk$nrow == chk$n), "\n")
cat("row letters spell the cached sequence exactly:", sum(chk$same_seq %in% TRUE), "\n")
cat("non-NA index always equals row position:", sum(chk$idx_eq_pos %in% TRUE), "\n")
cat("precursors with some NA index:", sum(chk$n_na_idx > 0, na.rm=TRUE), " rows:", sum(chk$n_na_idx, na.rm=TRUE), "\n")
bad <- chk %>% filter(!(same_seq %in% TRUE) | !(idx_eq_pos %in% TRUE))
cat("\nexceptions:", nrow(bad), "\n"); print(head(bad, 15))
## channel-level NA after positional indexing, cached precursors only
sub <- x[x$accession %in% prec$accession, ]
cat("\nNA share per channel (all residues of cached precursors):\n")
print(round(100*colMeans(is.na(sub[, scale_cols])), 2))
g <- data.frame(cons = is.na(sub$cons_rs_n), afm = is.na(sub$min_afm) | is.na(sub$mean_afm), struct = is.na(sub$relASA))
cat("\njoint missingness pattern (cons, afm, struct):\n"); print(round(100*prop.table(table(g)), 2))
