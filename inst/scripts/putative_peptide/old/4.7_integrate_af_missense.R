


secretome <- readRDS("~/peptide_alg/secretome_4.6.rds")

af_missense <- data.table::fread("~/peptide_alg/alpha_missense/AlphaMissense_aa_substitutions.tsv.gz")

af_missense %>%
   filter(uniprot_id %in% secretome[["accession"]]) %>%
   as_tibble %>%
   saveRDS(., "~/peptide_alg/af_missense_secretome.rds")

af_missense <- readRDS("~/peptide_alg/af_missense_secretome.rds")


af_missense <- af_missense %>%
  as_tibble %>%
  group_by(uniprot_id) %>%
  summarize(af_missense = list(tibble(protein_variant = protein_variant,
                                      am_pathogenicity = am_pathogenicity,
                                      am_class = am_class))) %>%
  dplyr::rename(accession = uniprot_id)



format_af_missense <- function(x) {
  
  x %>%
    mutate(AA = str_extract(protein_variant, "^."),
           variant = str_extract(protein_variant, ".$")) %>%
    mutate(group = (row_number() - 1) %/% 19 + 1) %>%
    {.[["AA"]][1] ->> miss_aa; .} %>%
    add_row(protein_variant = paste0(miss_aa, 1, miss_aa), 
            am_pathogenicity = NA,
            am_class = NA,
            AA = miss_aa,
            variant = miss_aa,
            group = 1,
            .before = 1) %>%
    pivot_wider(names_from = variant, id_cols = c("group", "AA"), values_from = am_pathogenicity) %>%
    dplyr::select(any_of(c("AA", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))) %>%
    mutate(mean_af_missense = rowMeans(across(-AA), na.rm = TRUE), .after = AA)
}


af_missense <- af_missense %>% 
                  mutate(af_missense = map(af_missense, format_af_missense))

af_missense <- af_missense %>%
  mutate(sequence_afm = map_chr(af_missense, ~paste(.[["AA"]], collapse = "")))

secretome <- left_join(secretome, af_missense, by = "accession")

secretome %>%
  mutate(same_seq = sequence_uni == sequence_afm) %>% pull(same_seq) %>% sum(., na.rm = TRUE)

library(Biostrings)

start <- Sys.time()

secretome <- secretome %>%
  mutate(af_missense_mapped = pmap(.l = list(seq1 = sequence_afm, seq2 = sequence_uni, to_map = af_missense),
                          .f = map_table))

end <- Sys.time()
end - start

saveRDS(secretome, "~/peptide_alg/secretome_4.7.rds")

