suppressMessages({library(dplyr); library(purrr)})
ROOT <- "/Users/kbrulois/R_projects/ligandFinder/.claude/worktrees/modest-satoshi-87b2c8"
A2 <- path.expand("~/AF2_analysis")
id_map <- readRDS(file.path(ROOT, "data", "id_mapping.rds"))
e2s <- setNames(id_map[["Gene Names (primary)"]], id_map[["Entry Name"]])
is_chemo <- function(g) grepl("^(CCL|CXCL|XCL|CX3CL)", g)
K <- readRDS(file.path(commandArgs(TRUE)[1], "splits_K.rds"))
trained <- unique(paste(K$gene, K$pep_id))

## ---- 1. CHEMOKINES: docked, but dropped from the set entirely --------------
con <- readRDS(file.path(A2, "knowns.rds")) %>%
  filter(location == "relevant") %>% filter(rank == 1, .by = afpd_dir_name) %>%
  mutate(gene = unname(e2s[p2_name]), accession = p2_id,
         pep_id = paste0(gene, "_", p2_range),
         pep_start = as.integer(stringr::str_extract(p2_range, "^\\d+")),
         pep_end   = as.integer(stringr::str_extract(p2_range, "\\d+$"))) %>%
  group_by(pep_id) %>% slice(1) %>% ungroup()
ch <- con %>% filter(is_chemo(gene)) %>% distinct(gene, accession, pep_id, pep_start, pep_end)
cat(sprintf("=== CHEMOKINES: docked peptides dropped from the training set ===\n"))
cat(sprintf("peptides: %d over %d precursors; ends = %d (2 per peptide)\n",
            nrow(ch), n_distinct(ch$accession), 2*nrow(ch)))
cat("  in the negative pool?  10_7 does `prec %>% filter(!is_chemokine(gene))` -> NO\n")
cat("  in the positives?      filtered out at step 1 -> NO\n")
print(ch %>% count(gene, name = "peptides") %>% as.data.frame(), row.names = FALSE)

## ---- 2. ligand_list / UniProt peptide ends never used as positives ---------
ll <- readRDS(file.path(ROOT, "inst", "extdata", "ligand_list.rds"))
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache_allres.rds")); prec <- fc$prec
ll2 <- ll %>% filter(!is.na(accession), accession %in% prec$accession) %>%
  mutate(gene = prec$gene[match(accession, prec$accession)],
         pid = paste0(gene, "_", start, "x", end)) %>%
  filter(!paste(gene, pid) %in% trained)
cat(sprintf("\n=== ligand_list peptides in the secretome but NEVER a training positive ===\n"))
cat(sprintf("peptides: %d over %d precursors (%d chemokine)\n",
            nrow(ll2), n_distinct(ll2$accession), sum(is_chemo(ll2$gene))))
cat("  BUT: their genes are in `banned` (250 genes), so they were kept OUT of the\n")
cat("  negatives too -- never seen as either class.\n")
print(head(ll2 %>% transmute(gene, pep = pid, start, end) %>% arrange(gene), 15) %>% as.data.frame(), row.names = FALSE)
cat(sprintf("  ... %d more\n", max(0, nrow(ll2) - 15)))
