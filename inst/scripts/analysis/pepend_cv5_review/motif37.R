suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]
source("R/pepend_windows.R")
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
fc <- readRDS("~/AF2_analysis/lf_pepend_residue_cache_allres.rds")
seqs <- setNames(fc$prec$seq, fc$prec$accession)
out <- tst %>% rowwise() %>%
  mutate(s = seqs[[accession]],
         next3 = substring(s, end + 1L, end + 3L),
         motif = lf_pepend_motif(list(s), "C", list(end), list(1L))) %>%
  ungroup() %>%
  transmute(gene, name, end, len, after_end = ifelse(nchar(next3) > 0, next3, "<protein end>"), motif)
cat("=== what follows each test peptide's C-terminal end ===\n")
print(as.data.frame(out %>% arrange(motif, gene)), row.names = FALSE)
cat("\n=== counts ===\n"); print(table(out$motif))
cat(sprintf("\nG + dibasic (amidation motif): %d of %d\n", sum(out$motif=="G + dibasic"), nrow(out)))
write.csv(out, file.path(S2, "testset_cterm_motif.csv"), row.names = FALSE)
