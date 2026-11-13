## per-residue min_afm vs relASA for every residue of every known peptide
suppressMessages({library(dplyr); library(tidyr); library(readr); library(purrr)})
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
fc  <- readRDS(file.path(A2, "lf_pepend_residue_cache.rds"))
pk  <- readRDS(file.path(A2, "lf_pocket_by_peptide.rds"))
cat("peptides in the pocket cache:", nrow(pk), "\n")

d <- pmap_dfr(pk, function(accession, gene, pep_id, pep_start, pep_end, pocket) {
  m <- fc$feats[[accession]]
  if (is.null(m)) return(NULL)
  r <- seq.int(pep_start, min(pep_end, nrow(m)))
  tibble(accession, gene, pep_id, resno = r,
         min_afm = m[r, "min_afm"], relASA = m[r, "relASA"],
         class = ifelse(r %in% unlist(pocket), "pep_pocket", "pep_other"))
})
## relASA is exactly 0 at EVERY residue of 716 precursors -- missing DSSP, not
## burial. Flag those peptides rather than letting them pile on the y = 0 line.
dead <- names(which(vapply(fc$feats, function(m) all(m[, "relASA"] == 0), logical(1))))
d <- d %>% mutate(no_dssp = accession %in% dead)
cat(sprintf("residues: %d from %d peptides / %d precursors\n", nrow(d),
            n_distinct(d$pep_id), n_distinct(d$accession)))
cat(sprintf("  pep_pocket %d, pep_other %d\n", sum(d$class == "pep_pocket"),
            sum(d$class == "pep_other")))
cat(sprintf("  in a no-DSSP precursor: %d residues (%d peptides)\n",
            sum(d$no_dssp), n_distinct(d$pep_id[d$no_dssp])))
cat("\nmin_afm:"); print(summary(d$min_afm))
cat("relASA:");  print(summary(d$relASA))
cat("\nexactly-zero relASA among the rest:",
    sum(d$relASA == 0 & !d$no_dssp), "of", sum(!d$no_dssp), "\n")
write_csv(d, file.path(S2, "pep_afm_asa.csv"))
