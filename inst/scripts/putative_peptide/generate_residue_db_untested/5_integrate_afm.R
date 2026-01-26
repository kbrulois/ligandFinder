

library(tidyverse)
library(ligandFinder)

uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_2.rds"))

af_missense <- data.table::fread(paste0(s_localDir, "/raw/alpha_missense/AlphaMissense_aa_substitutions.tsv.gz"))

af_missense <- af_missense %>%
  as_tibble %>%
  group_by(uniprot_id) %>%
  summarize(af_missense = list(tibble(protein_variant = protein_variant,
                                      am_pathogenicity = am_pathogenicity,
                                      am_class = am_class))) %>%
  dplyr::rename(accession = uniprot_id)

afm_columns <- c("AA", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

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
    tibble::add_column(bind_cols(Map(\(x) as.numeric(rep(NA, nrow(.))), afm_columns[!afm_columns %in% names(.)]))) %>%
    dplyr::select(any_of(afm_columns)) %>%
    mutate(mean_afm = rowMeans(across(-AA), na.rm = TRUE), .after = AA)
}


af_missense <- af_missense %>%
  mutate(af_missense = map(af_missense, format_af_missense))

af_missense <- af_missense %>%
  mutate(sequence_afm = map_chr(af_missense, ~paste(.[["AA"]], collapse = "")))

uniprot_t <- left_join(uniprot_t, af_missense, by = "accession")

uniprot_t %>%
  mutate(same_seq = sequence_uni == sequence_afm) %>% pull(same_seq) %>% sum(., na.rm = TRUE)

library(Biostrings)

start <- Sys.time()

uniprot_t <- uniprot_t %>%
  mutate(af_missense = pmap(.l = list(seq1 = sequence_afm,
                                      seq2 = sequence_uni,
                                      to_map = af_missense),
                            .f = map_table))
end <- Sys.time()
end - start

uniprot_t <- uniprot_t %>%
  mutate(afm_map_score = map_dbl(af_missense, `[[`, "score")) %>%
  mutate(af_missense = map(af_missense, `[[`, "ms"))


summarize_afm <- function(x) {
  if(is.data.frame(x)) {
    x %>%
      rowwise() %>%
      mutate(min_afm = min(c_across(!any_of(c("AA", "mean_afm", "mean_af_missense"))), na.rm = TRUE), .after = "AA") %>%
      mutate(max_afm = max(c_across(!any_of(c("AA", "mean_afm", "mean_af_missense"))), na.rm = TRUE), .after = "AA") %>%
      mutate(max_AA = names(.)[!names(.) %in% c("AA", "mean_afm", "mean_af_missense", "max_afm")][which.max(c_across(!any_of(c("AA", "mean_afm", "mean_af_missense", "max_afm"))))], .after = "AA") %>%
      mutate(min_AA = names(.)[!names(.) %in% c("AA", "mean_afm", "mean_af_missense", "min_afm", "max_AA")][which.min(c_across(!any_of(c("AA", "mean_afm", "mean_af_missense", "min_afm", "max_AA"))))], .after = "AA") %>%
      mutate(across(all_of(c("min_afm", "max_afm", "max_AA", "min_AA")), ~if_else(is.na(AA), NA, .)))
  } else {
    x
  }
}

start <- Sys.time()

options(future.globals.maxSize = Inf)
future::plan(strategy = future::multisession(workers = 8))
#future::plan(strategy = future::sequential())


uniprot_t$af_missense <- furrr::future_map(uniprot_t$af_missense,
                                           summarize_afm)

end <- Sys.time()
end - start


saveRDS(uniprot_t, paste0(s_localDir, "/processed/uniprot_3.rds"))

