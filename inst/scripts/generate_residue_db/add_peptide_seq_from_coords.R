


ligs <- openxlsx::read.xlsx("~/Desktop/things_to_test.xlsx")




ligs <- ligs %>% as_tibble




pq_path <- "~/ligandFinder_data/residue_db"

res_db <- arrow::open_dataset(source = pq_path)

residue_data <- res_db %>%
  filter(gene %in% ligs$ligand) %>%
  select(gene, sequence_uni) %>%
  collect()

ligs <- left_join(ligs, residue_data, by = join_by(ligand == gene))

get_pep_seq <- function(x, sequence_uni) {

  coords <- stringr::str_split(x, "-", simplify = TRUE) %>%
    as_tibble() %>%
    rename(start = V1, end = V2) %>%
    mutate(across(everything(), as.numeric))

  stringr::str_sub(sequence_uni, start = coords[[1]], end = coords[[2]])

}


ligs <- ligs %>%
  mutate(pep_sequence = get_pep_seq(ligand.coordinates, sequence_uni), .after = ligand.coordinates) %>%
  select(-sequence_uni)

openxlsx::write.xlsx(ligs, file = "~/Desktop/to_test_table.xlsx")









