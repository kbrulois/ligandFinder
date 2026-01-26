
library(ligandFinder)
library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"
uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_6.rds"))



secretome <- uniprot_t %>%
  mutate(location = case_when(hpa_sloc == 2 & uniprot_loc_min == 2 ~ "4l",
                              hpa_sloc == 2 & uniprot_loc_max == 2 ~ "3l",
                              hpa_sloc == 2 | uniprot_loc_max == 2 ~ "2l",
                              hpa_sloc == 0 & uniprot_loc_max == 1 ~ "1l",
                              uniprot_topo_max == 2 & uniprot_loc_max == 2 ~ "4t",
                              (uniprot_topo_max == 2 | hpa_sloc == 2 ) & uniprot_loc_max == 1 ~ "3t",
                              uniprot_topo_max == 2 & uniprot_loc_max == 0 ~ "2t",
                              uniprot_topo_max == 1 ~ "1t",
                              TRUE ~ "IC"))

table(secretome$location)

secretome <- secretome %>%
  filter(location != "IC" & !is.na(sequence_uni))

rm(uniprot_t)
gc()





features <- uniprot_t %>% filter(accession == "P01275") %>% pull(features) %>% `[[`(1)

sequence_uni <- uniprot_t %>% filter(accession == "P01275") %>% pull(sequence_uni)







expand_by_residue <- function(x, dat_to_expand = c("topo", "features_expanded", "dssp", "af_missense", "cons", "alignment_AA")) {

x <- x %>%
  mutate(features_expanded = map2(features, sequence_uni, expand_features))

x <- x %>%
  mutate(topo = map2(sequence_uni, topo, \(x, y) tibble(AA = str_split(x, "", simplify = TRUE) %>% c,
                                                        topo = str_split(y, "", simplify = TRUE) %>% c))) %>%
  mutate(to_expand = pmap(rlang::syms(dat_to_expand),
                          bind_cols, .name_repair = "minimal")) %>%
  mutate(to_expand = map(to_expand, \(x) x[, !duplicated(colnames(x))])) %>%
  mutate(to_expand = map(to_expand, \(x) x[, colnames(x) != ""]))

to_return <- x %>%
  select(-where(is.list), to_expand) %>%
  unnest(to_expand)

to_return <- to_return %>%
  mutate(topo2 = if_else(has_topo & ("e" %in% topo), topo, "e"))

return(to_return)

}


secretome_aa <- expand_by_residue(secretome)


secretome_aa <- secretome_aa %>%
  filter(topo2 == "e" & topo != "s")





