###combine uniprot and aminode ~20min


uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot.rds"))

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


aminode <- readRDS(paste0(s_localDir, "/processed/aminode.rds"))

aminode[aminode$gene == "C10ORF99", "gene"] <- "GPR15LG"

aminode <- aminode %>%
  mutate(cons = map(.x = cons, .f = \(x) {
    x %>%
      filter(AA != "N/A" | index != "N/A") %>%
      mutate(index = as.integer(index),
             cons = as.numeric(cons)) %>%
      mutate(frequency = cons/mean(cons)) %>%
      mutate(window = slider::slide_dbl(frequency, mean, .before = 5, .after = 5)) %>%
      mutate(smooth = slider::slide_dbl(window, mean, .before = 3, .after = 3)) %>%
      mutate(doubleSmooth = slider::slide_dbl(smooth, mean, .before = 3, .after = 3)) %>%
      mutate(across(all_of(c("frequency", "doubleSmooth")), ~scales::rescale(log(. + 0.01), c(1,0)), .names = "{.col}_scaled"))
    
  })) %>%
  mutate(sequence = map_chr(.x = cons, .f = ~ paste(.[["AA"]], collapse = "")))

secretome <- left_join(secretome, aminode, by = "gene", suffix = c("_uni", "_ami"))


library(Biostrings)

start <- Sys.time()

future::plan(strategy = future::sequential())

secretome <- secretome %>%
  mutate(cons = furrr::future_pmap(.l = list(seq1 = sequence_ami, seq2 = sequence_uni, to_map = cons),
                                   .f = map_table))

gc()

end <- Sys.time()

end - start

saveRDS(secretome, paste0(s_localDir, "/processed/secretome_1.rds"))

