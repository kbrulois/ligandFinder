



id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))


dat <- tibble(files = fs::dir_ls()) %>%
          filter(stringr::str_detect(files, ".pkl.xz$")) %>%
          mutate(uniprot_id = stringr::str_remove(files, ".pkl.xz$")) %>%
          mutate(uniprot_id = setNames(id_map$Entry, id_map$`Entry Name`)[uniprot_id])

dat2 <- data.table::fread("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/n_c_term_candidates_May_2025.csv")

dat2 <- dat2 %>%
          as_tibble %>%
          mutate(parse_proteins(file_name = model_id, delim_proteins= ";", delim_ranges = ",", delim_start_end="-"))

dat2 <- dat2 %>%
  mutate(proteins_in_msa_db = case_when(!p1_id %in% dat$uniprot_id ~ FALSE,
                                        !p2_id %in% dat$uniprot_id ~ FALSE,
                                        TRUE ~ TRUE))
