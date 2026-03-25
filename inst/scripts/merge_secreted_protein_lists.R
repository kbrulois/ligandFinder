



####check secreted proteins list

uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_6.rds"))




sg1 <- data.table::fread("~/AF2_analysis/ECB_hSecreted_GN_wPep.csv") %>% as_tibble %>% mutate(source = "ecb_irina")

sg2 <- data.table::fread("~/AF2_analysis/not_in_ECB_hSecreted_GN.csv") %>% as_tibble %>% mutate(source = "irina")

sg2 <- bind_rows(sg1, sg2)


genes1 <- sg2 %>% filter(source == "ecb_irina") %>% pull(GN)

genes2 <- sg2 %>% filter(source == "irina") %>% pull(GN)

genes3 <- sg2 %>% pull(GN)

test <- uniprot_t %>%
  mutate(location = case_when(hpa_sloc == 2 & uniprot_loc_min == 2 ~ "4l",
                              hpa_sloc == 2 & uniprot_loc_max == 2 ~ "3l",
                              hpa_sloc == 2 | uniprot_loc_max == 2 ~ "2l",
                              hpa_sloc == 0 & uniprot_loc_max == 1 ~ "1l",
                              uniprot_topo_max == 2 & uniprot_loc_max == 2 ~ "4t",
                              (uniprot_topo_max == 2 | hpa_sloc == 2 ) & uniprot_loc_max == 1 ~ "3t",
                              uniprot_topo_max == 2 & uniprot_loc_max == 0 ~ "2t",
                              uniprot_topo_max == 1 ~ "1t",
                              TRUE ~ "IC")) %>%
  select(gene, accession, uniprot_loc, `Subcellular location`, `Secretome location`, `Secretome function`,
         `Subcellular main location`, `Subcellular additional location`,
         location, uniprot_loc_max, uniprot_loc_mean,
         uniprot_loc_min, uniprot_topo_max, uniprot_topo_mean, uniprot_topo_min, has_TM,

         hpa_loc, hpa_sloc) %>%
  filter(gene %in% c(secretome$gene, sg2$GN)) %>%
  mutate(source_k = case_when(location %in% c("4l", "3l", "2l", "1l") ~ "kevin",
                            location %in% c("4t", "3t", "2t", "1t") ~ "kevin_tm",
                            TRUE ~ "not")) %>%
  mutate(source = case_when(gene %in%  genes1 & source_k %in% c("kevin", "kevin_tm") ~ "everyone",
                            gene %in%  genes2 & source_k %in% c("kevin", "kevin_tm") ~ "irina_kevin",
                            gene %in% genes3 & source_k == "not" ~ "ecb_irina",
                            TRUE ~ source_k))


not_in <- tibble(gene = sg2$GN[!sg2$GN %in% uniprot_t$gene],
       source = "ecb_irina")


test <- bind_rows(test, not_in)


table(test$source)

test <- test %>%
          arrange(source) %>%
          relocate(source, source_k, .after = accession)

data.table::fwrite(test, "~/Desktop/extracellular_proteins_v2.csv")




