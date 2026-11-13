suppressMessages(library(dplyr))
A2 <- path.expand("~/AF2_analysis"); S2 <- commandArgs(TRUE)[1]
is_chemo <- function(g) grepl("^(CCL|CXCL|XCL|CX3CL)", g)
K <- readRDS(file.path(S2, "splits_K.rds")); tg <- unique(K$gene)
fc <- readRDS(file.path(A2, "lf_pepend_residue_cache_allres.rds")); prec <- fc$prec
ll <- readRDS("inst/extdata/ligand_list.rds")
d <- data.frame(accession = ll$accession, start = ll$start, end = ll$end,
                name = ll$final_name, stringsAsFactors = FALSE) %>%
  filter(!is.na(accession), accession %in% prec$accession, !is.na(start), !is.na(end)) %>%
  mutate(gene = prec$gene[match(accession, prec$accession)], len = end - start + 1L) %>%
  filter(!gene %in% tg) %>% distinct(gene, start, end, .keep_all = TRUE)
cat(sprintf("gene-level-clean ligand_list peptides: %d (chemokine %d, non-chemokine %d)\n",
            nrow(d), sum(is_chemo(d$gene)), sum(!is_chemo(d$gene))))
nc <- d %>% filter(!is_chemo(gene))
cat(sprintf("non-chemokine: median len %d aa (range %d-%d); <=60 aa: %d; <=40 aa: %d\n",
            median(nc$len), min(nc$len), max(nc$len), sum(nc$len <= 60), sum(nc$len <= 40)))
write.csv(nc %>% arrange(len) %>% select(gene, accession, name, start, end, len),
          file.path(S2, "test_set_candidates.csv"), row.names = FALSE)
cat("\n=== short named peptides (<= 40 aa), gene never a training positive ===\n")
print(nc %>% filter(len <= 40) %>% arrange(gene) %>%
        transmute(gene, name, span = paste0(start, "-", end), len) %>% as.data.frame(), row.names = FALSE)
