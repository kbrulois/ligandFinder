suppressMessages({library(dplyr)})
ROOT <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
A2 <- path.expand("~/AF2_analysis"); S2 <- commandArgs(TRUE)[1]
is_chemo <- function(g) grepl("^(CCL|CXCL|XCL|CX3CL)", g)
K <- readRDS(file.path(S2, "splits_K.rds"))
trained_genes <- unique(K$gene)
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache_allres.rds")); prec <- fc$prec
ll <- readRDS(file.path(ROOT, "inst", "extdata", "ligand_list.rds")) %>%
  filter(!is.na(accession), accession %in% prec$accession) %>%
  mutate(gene = prec$gene[match(accession, prec$accession)]) %>%
  filter(!gene %in% trained_genes)          # gene-level clean
cat("=== ligand_list peptides whose GENE has no training positive at all ===\n")
cat(sprintf("peptides %d over %d precursors | chemokine %d | non-chemokine %d\n",
    nrow(ll), n_distinct(ll$accession), sum(is_chemo(ll$gene)), sum(!is_chemo(ll$gene))))
nc <- ll %>% filter(!is_chemo(gene))
cat(sprintf("  of the non-chemokine ones, median length %d aa (range %d-%d)\n",
    median(nc$end - nc$start + 1), min(nc$end-nc$start+1), max(nc$end-nc$start+1)))
cat(sprintf("  short (<= 60 aa, i.e. peptide-like): %d\n", sum((nc$end-nc$start+1) <= 60)))
print(head(nc %>% filter(end-start+1 <= 60) %>%
  transmute(gene, pep = paste0(gene,"_",start,"x",end), len = end-start+1,
            name = if ("final_name" %in% names(.)) final_name else NA) %>% arrange(gene), 20) %>%
  as.data.frame(), row.names = FALSE)
write.csv(nc, file.path(S2, "test_candidates_clean.csv"), row.names = FALSE)

## per-residue split, for the tables
z <- reticulate::import("numpy", convert = TRUE)$load(
  file.path(A2, "lf_resid_labels_allres.npz"), allow_pickle = TRUE)
cat(sprintf("\n=== lf_resid (per-residue) split ===\nlabelled precursors %d in %d similarity groups; val_frac 0.4 on ROWS\n",
    sum(z[["labelled"]]), length(unique(z[["prec_group"]][z[["labelled"]]==1]))))
