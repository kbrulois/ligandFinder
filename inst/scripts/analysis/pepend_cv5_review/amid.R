suppressMessages({library(dplyr)})
S2 <- commandArgs(TRUE)[1]
tst <- read.csv(file.path(S2,"test_set_tiered.csv"), stringsAsFactors=FALSE) %>%
  filter(tier != "contaminated (negatives drawn)")
sec <- readRDS("~/AF2_analysis/secretome_latest.rds")
fc  <- readRDS("~/AF2_analysis/lf_pepend_residue_cache_allres.rds")
ROOT <- getwd(); source(file.path(ROOT,"R","pepend_windows.R"))
out <- bind_rows(lapply(seq_len(nrow(tst)), function(i) {
  a <- tst$accession[i]; e <- tst$end[i]
  j <- match(a, sec$accession); s <- fc$prec$seq[match(a, fc$prec$accession)]
  f <- sec$features[[j]]
  am <- NA_character_
  if (is.data.frame(f) && nrow(f)) {
    mr <- f %>% filter(type == "modified residue", !is.na(start), start == e,
                       grepl("amide", description, ignore.case = TRUE))
    if (nrow(mr)) am <- mr$description[1]
  }
  nxt <- substring(s, e+1L, e+3L)
  tibble(gene = tst$gene[i], name = tst$name[i], end = e, len = tst$len[i],
         next3 = ifelse(nchar(nxt)>0, nxt, "<end of protein>"),
         motif = lf_pepend_motif(list(s), "C", list(e), list(1L)),
         G_at_plus1 = substring(s, e+1L, e+1L) == "G",
         uniprot_amide = am)
}))
out$amidated <- !is.na(out$uniprot_amide)
cat(sprintf("=== C-terminal amidation among the %d test peptides ===\n", nrow(out)))
cat(sprintf("UniProt 'amide' annotated at the peptide's own C-terminal residue: %d of %d (%.0f%%)\n",
            sum(out$amidated), nrow(out), 100*mean(out$amidated)))
cat(sprintf("glycine immediately after the end (the PAM donor): %d\n", sum(out$G_at_plus1)))
cat(sprintf("agreement between the two: %d\n\n", sum(out$amidated == out$G_at_plus1)))
print(out %>% filter(amidated) %>% arrange(gene) %>%
  transmute(gene, name, end, next3, motif, uniprot_amide) %>% as.data.frame(), row.names=FALSE)
cat("\n--- G at +1 but NOT annotated amidated ---\n")
print(out %>% filter(G_at_plus1, !amidated) %>% transmute(gene, name, next3, motif) %>% as.data.frame(), row.names=FALSE)
cat("\n--- motif class of all 37 C ends ---\n"); print(table(out$motif))
cat("\n--- amidated x motif ---\n"); print(table(amidated = out$amidated, motif = out$motif))
write.csv(out, file.path(S2, "testset_amidation.csv"), row.names=FALSE)
