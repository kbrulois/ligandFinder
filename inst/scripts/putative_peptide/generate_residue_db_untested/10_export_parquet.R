


pq_path <- paste0(s_localDir, "/processed/secretome_parquet3")

dir.create(pq_path)

lists_to_unpack <- c("cons", "dssp", "af_missense", "af_xyz")

secretome <- secretome %>%
  mutate(across(all_of(lists_to_unpack), 
                .fns = list(tagtoremove = ~map(.x, .f = ~`[[`(., "ms")),
                            score = ~map_dbl(.x, .f = ~`[[`(., "score"))),
                .unpack = TRUE)) %>%
  select(-all_of(lists_to_unpack)) %>%
  rename_with(.cols = ends_with("_tagtoremove"), .fn = ~sub("_tagtoremove", "", .))


clean_list_cols <- \(x) {
  
  not_tibble <- !sapply(x, is_tibble)
  
  if(any(not_tibble)) {
    dummy_data <- x[!not_tibble][[1]] %>% filter(FALSE)
    x[not_tibble] <- map(sum(not_tibble), \(x) dummy_data)
  }
  
  x
}

secretome %>%
  mutate(across(lists_to_unpack, .fns = clean_list_cols)) %>%
  mutate(gene_grp = stringr::str_sub(gene, 1, 1)) %>%
  group_by(gene_grp) %>%
  select(-annotations) %>%
  arrow::write_dataset(path = pq_path, format = "parquet")

