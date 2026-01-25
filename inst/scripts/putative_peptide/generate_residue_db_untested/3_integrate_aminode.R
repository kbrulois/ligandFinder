###combine uniprot and aminode ~20min

uniprot_f <- data.table::fread(fs::path(s_localDir, "processed/uniprot_features.csv")) %>% as_tibble
uniprot_a <- data.table::fread(fs::path(s_localDir, "processed/uniprot_annotations.csv")) %>% as_tibble

uniprot_f <- uniprot_f %>%
              nest(features = !any_of(c("accession", "gene", "full_name", "sequence")))

uniprot_a <- uniprot_a %>%
                nest(annotations = !any_of(c("accession"))) %>%
                filter(accession %in% uniprot_f[["accession"]])

uniprot_t <- left_join(uniprot_f, uniprot_a, by = "accession")



aminode <- readRDS(paste0(s_localDir, "/processed/aminode.rds"))

aminode <- aminode %>%
            select(-alignment, -sim_mat)

aminode[aminode$gene == "C10ORF99", "gene"] <- "GPR15LG"

aminode <- aminode %>%
  mutate(cons = map2(.x = cons, .y = alignment_AA, .f = \(x, y) {

    x <- bind_cols(x, y)

    x %>%
      filter(AA != "N/A", index != "N/A") %>%
      mutate(frequency = cons_og/mean(cons_og)) %>%
      mutate(window = slider::slide_dbl(frequency, mean, .before = 5, .after = 5)) %>%
      mutate(smooth = slider::slide_dbl(window, mean, .before = 3, .after = 3)) %>%
      mutate(doubleSmooth = slider::slide_dbl(smooth, mean, .before = 3, .after = 3)) %>%
      mutate(across(all_of(c("frequency", "doubleSmooth")), ~scales::rescale(log(. + 0.01), c(1,0)), .names = "{.col}_lrs")) %>%
      mutate(across(all_of(c("frequency", "doubleSmooth")), ~scales::rescale(., c(1,0)), .names = "{.col}_rs"))

  })) %>%
  mutate(sequence = map_chr(.x = cons, .f = ~ paste(.[["AA"]], collapse = "")))

uniprot_t <- left_join(uniprot_t, aminode, by = "gene", suffix = c("_uni", "_ami"))


library(Biostrings)

start <- Sys.time()

future::plan(strategy = future::multisession(workers = 8))
#future::plan(strategy = future::sequential())

uniprot_t <- uniprot_t %>%
  mutate(cons = furrr::future_pmap(.l = list(seq1 = sequence_ami, seq2 = sequence_uni, to_map = cons),
                                   .f = map_table))

gc()

end <- Sys.time()

end - start

uniprot_t <- uniprot_t %>%
  mutate(ami_map_score = map_dbl(cons, `[[`, "score"), .before = "features") %>%
  mutate(cons = map(cons, `[[`, "ms")) %>%
  mutate(uni_ami_seq_len_dif = map2_dbl(cons, alignment_AA, \(x, y) {
    if(is.null(y)) {return(NA)} else {
    return(abs(nrow(x) - nrow(y)))
    }
  }), .after = "ami_map_score")

uniprot_t <- uniprot_t %>%
  mutate(map2(cons, alignment_AA, \(x, y) {
      if(is.data.frame(x) & is.data.frame(y)) {
        species_cols <- colnames(y)
        tibble(cons = list(x %>%
                  select(!all_of(species_cols))),
               alignment_AA = list(x %>%
                                     select(all_of(species_cols))))
      } else {
        tibble(cons = NA,
               alignment_AA = NA)

      }

  }) %>% bind_rows)


saveRDS(uniprot_t, paste0(s_localDir, "/processed/uniprot_1.rds"))
















secretome_og <- data.table::fread(paste0(s_localDir, "/raw/sa_location_Secreted HPA 2793 genes-1.tsv"))

secretome <- uniprot_t %>%
  mutate(secreted_any = map_lgl(.x = annotations, .f = \(x) {
    x %>%
      filter(annotation_name == "comment" &
               annotation_type == "subcellular location" &
               name_2 == "location") %>%
      {any(.[["annotation"]] == "Secreted")}
  })) %>%
  mutate(secreted_all = map_lgl(.x = annotations, .f = \(x) {
    x %>%
      filter(annotation_name == "comment" &
               annotation_type == "subcellular location" &
               name_2 == "location") %>%
      {all(.[["annotation"]] == "Secreted")}
  })) %>%
  mutate(secreted_HPA = gene %in% secretome_og[["Gene"]]) %>%
  mutate(secreted_final = secreted_any | secreted_HPA)

rm(uniprot_t, secretome_og)


