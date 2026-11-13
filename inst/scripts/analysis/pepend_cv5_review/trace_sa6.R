suppressMessages(library(dplyr))
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds"); ap <- fc$all_params3
need <- setdiff(ap, c("AA_hydro","AA_charge","AA_mw","AA_pI","padding"))
for (f in c("secretome_aa_1.rds")) {
  x <- readRDS(file.path("~/peptide_alg/build_residue_db/processed", f))
  cat(f, ":", nrow(x), "rows,", ncol(x), "cols,", length(unique(x$accession)), "accessions\n")
  cat("  channels the cache needs that are missing:", paste(setdiff(need, names(x)), collapse = ", "), "\n")
  cat("  topo2:\n"); print(table(x$topo2, useNA = "ifany"))
  cat("  NA index:", sum(is.na(x$index)), "\n")
  cat("  cached precursors covered:", sum(fc$prec$accession %in% x$accession), "of", nrow(fc$prec), "\n")
  cat("  residues in cache:", sum(nchar(fc$prec$seq)), "\n")
}
