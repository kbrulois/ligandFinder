
library(ligandFinder)
library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"
uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_6.rds"))

secretome_genes <- data.table::fread("~/Desktop/Peptides/extracellular_proteins_v2.csv") %>% as_tibble



secretome <- uniprot_t %>%
  mutate(location = case_when(gene == "BRINP1" ~ "4l",
                              hpa_sloc == 2 & uniprot_loc_min == 2 ~ "4l",
                              hpa_sloc == 2 & uniprot_loc_max == 2 ~ "3l",
                              hpa_sloc == 2 | uniprot_loc_max == 2 ~ "2l",
                              hpa_sloc == 0 & uniprot_loc_max == 1 ~ "1l",
                              uniprot_topo_max == 2 & uniprot_loc_max == 2 ~ "4t",
                              (uniprot_topo_max == 2 | hpa_sloc == 2 ) & uniprot_loc_max == 1 ~ "3t",
                              uniprot_topo_max == 2 & uniprot_loc_max == 0 ~ "2t",
                              uniprot_topo_max == 1 ~ "1t",
                              gene %in% secretome_genes[["gene"]] ~ "sec",
                              TRUE ~ "IC"))

table(secretome$location)

secretome <- secretome %>%
  filter(location != "IC" & !is.na(sequence_uni))

rm(uniprot_t)
gc()







expand_by_residue <- function(x, dat_to_expand = c("topo", "features_expanded", "dssp", "af_missense", "cons", "alignment_AA")) {

x <- secretome %>%
  mutate(features_expanded = map2(features, sequence_uni, expand_features))

x <- x %>%
  mutate(topo = map2(sequence_uni, topo, \(x, y) tibble(AA = str_split(x, "", simplify = TRUE) %>% c,
                                                        topo = str_split(y, "", simplify = TRUE) %>% c)))

x <- x %>%
  mutate(to_expand = pmap(pick(any_of(dat_to_expand)),
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


#secretome_aa <- secretome_aa %>%
#  filter(topo2 == "e" & topo2 != "s")




offlimits_features <- c("signal peptide", "E")

black_balled <- c("thrombin light chain", "epiregulin", "INSL5 (A chain)")


one_hot_encode <- function(x, tag = "") {
  x[is.na(x)] <- "unknown"
  unique_aa <- unique(x)
  to_return <- t(sapply(x, function(y) as.integer(unique_aa == y)))
  colnames(to_return) <- paste0(tag, "_", unique_aa)
  return(as_tibble(to_return))
}



AAs <- c("AA_Q", "AA_P", "AA_V", "AA_L", "AA_T", "AA_S", "AA_A",
         "AA_G", "AA_R", "AA_F", "AA_C", "AA_I", "AA_N", "AA_Y", "AA_W",
         "AA_K", "AA_D", "AA_E", "AA_M", "AA_H", "AA_U", "minAA_E", "minAA_A", "minAA_I",
         "minAA_V", "minAA_L", "minAA_P", "minAA_T", "minAA_S", "minAA_Q",
         "minAA_G", "minAA_D", "minAA_F", "minAA_N", "minAA_R", "minAA_H",
         "minAA_Y", "minAA_K", "minAA_M", "minAA_C", "minAA_W", "maxAA_C",
         "maxAA_D", "maxAA_K", "maxAA_W", "maxAA_P", "maxAA_F", "maxAA_M",
         "maxAA_E", "maxAA_H", "maxAA_L", "maxAA_I", "maxAA_G", "maxAA_Q",
         "maxAA_N", "maxAA_T", "maxAA_Y", "maxAA_R", "maxAA_S", "maxAA_A",
         "maxAA_V")

cons_metrics <- c("cons_rs", "cons_lrs", "blos_wt_all",
                  "blos_uw_all", "blos_wt_mam", "blos_uw_mam", "gran_wt_all", "gran_uw_all",
                  "gran_wt_mam", "gran_uw_mam", "blos_nr_all", "blos_nr_mam", "gran_nr_all",
                  "gran_nr_mam")

basic <- c(cons_metrics, paste0(cons_metrics, "_n"), c("max_afm", "min_afm", "mean_afm"), "relASA")

SS <- c("SS_P", "SS_S",
        "SS_E", "SS_-", "SS_T", "SS_G", "SS_B", "SS_H", "SS_I")

angles <- c("Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin")

energy <- c("NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy")

#sites <- c("DB_Dibasic",  "SV_sequence variant")






nn <- list()
nn[["nn1"]] <- c(basic, SS)
nn[["nn2"]] <- c(nn[["nn1"]], angles, energy)
nn[["nn4"]] <- c(nn[["nn2"]], AAs)


nn <- tibble(neural_net = names(nn),
             parameters = nn)

all_params <- unique(do.call(c, nn[["parameters"]]))





context <- c("_lag1", "_lag2", "_lead1", "_lead2")

nn2 <- nn %>%
  mutate(parameters = map(parameters, .f = ~c(., do.call(c, map(., ~paste0(., context)))))) %>%
  mutate(neural_net = paste0(neural_net, "c"))

names(nn2[["parameters"]]) <- paste0(names(nn2[["parameters"]]), "c")

nn <- bind_rows(nn, nn2)

nn <- nn %>%
  mutate(parameters = map(parameters, .f = ~c(., "known")))


if(FALSE) {
  to_remove <- "^SS_.*_(lead|lag)\\d+$"

  nn3 <- nn %>%
    filter(neural_net %in% c("nn5c", "nn6c", "nn7c")) %>%
    mutate(parameters = map(parameters, .f = ~.[!str_detect(., to_remove)]))

  names(nn3$parameters) <- nn3$neural_net <- paste0(nn3$neural_net, "noSS")

  nn <- bind_rows(nn, nn3)

  nn_og <- readRDS("~/peptide_alg/nn_models_256_64_4_leaky_relu_lion_opt_large_val_better_angles_w_context.rds")
  nn <- bind_cols(nn, nn_og[,3:5])
}


secretome_aa <- secretome_aa %>%
  mutate(SS = replace_na(SS, "-")) %>%
  mutate(one_hot_encode(AA, tag = "AA")) %>%
  mutate(one_hot_encode(SS, tag = "SS")) %>%
  mutate(one_hot_encode(sites_Dibasic_type, tag = "DB")) %>%
  mutate(one_hot_encode(`uniprot_sequence variant_type`, tag = "SV")) %>%
  mutate(one_hot_encode(gpcrdb_gtp_type, tag = "known")) %>%
  mutate(one_hot_encode(min_AA, tag = "minAA")) %>%
  mutate(one_hot_encode(max_AA, tag = "maxAA")) %>%
  mutate(across(c("Phi", "Psi"), .fns = ~tibble(cos = cos(. * pi/180),
                                                sin = sin(. * pi/180)), .unpack = TRUE)) %>%
  dplyr::group_by(accession) %>%
  mutate(across(ends_with("_energy"), .fns = \(x) {
    x[is.na(x)] <- 0
    scales::rescale(x, to = c(0, 1))})) %>%
  ungroup() %>%
  select(!ends_with("_unknown")) %>%
  mutate(known = rowSums(across(starts_with("known_gpcr_pep"))) > 0) %>%
  mutate(across(all_of(cons_metrics), ~ .x * species_limit, .names = "{.col}_n"))



















offlimits_features <- c("signal peptide", "E")

black_balled <- c("thrombin light chain", "epiregulin", "INSL5 (A chain)")


one_hot_encode <- function(x, tag = "") {
  x[is.na(x)] <- "unknown"
  unique_aa <- unique(x)
  to_return <- t(sapply(x, function(y) as.integer(unique_aa == y)))
  colnames(to_return) <- paste0(tag, "_", unique_aa)
  return(as_tibble(to_return))
}



AAs <- c("AA_Q", "AA_P", "AA_V", "AA_L", "AA_T", "AA_S", "AA_A",
         "AA_G", "AA_R", "AA_F", "AA_C", "AA_I", "AA_N", "AA_Y", "AA_W",
         "AA_K", "AA_D", "AA_E", "AA_M", "AA_H", "AA_U", "minAA_E", "minAA_A", "minAA_I",
         "minAA_V", "minAA_L", "minAA_P", "minAA_T", "minAA_S", "minAA_Q",
         "minAA_G", "minAA_D", "minAA_F", "minAA_N", "minAA_R", "minAA_H",
         "minAA_Y", "minAA_K", "minAA_M", "minAA_C", "minAA_W", "maxAA_C",
         "maxAA_D", "maxAA_K", "maxAA_W", "maxAA_P", "maxAA_F", "maxAA_M",
         "maxAA_E", "maxAA_H", "maxAA_L", "maxAA_I", "maxAA_G", "maxAA_Q",
         "maxAA_N", "maxAA_T", "maxAA_Y", "maxAA_R", "maxAA_S", "maxAA_A",
         "maxAA_V")

cons_metrics <- c("cons_rs", "cons_lrs", "blos_wt_all",
                  "blos_uw_all", "blos_wt_mam", "blos_uw_mam", "gran_wt_all", "gran_uw_all",
                  "gran_wt_mam", "gran_uw_mam", "blos_nr_all", "blos_nr_mam", "gran_nr_all",
                  "gran_nr_mam")

basic <- c(cons_metrics, paste0(cons_metrics, "_n"), c("max_afm", "min_afm", "mean_afm"), "relASA")

SS <- c("SS_P", "SS_S",
        "SS_E", "SS_-", "SS_T", "SS_G", "SS_B", "SS_H", "SS_I")

angles <- c("Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin")

energy <- c("NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy")

#sites <- c("DB_Dibasic",  "SV_sequence variant")






nn <- list()
nn[["nn1"]] <- c(basic, SS)
nn[["nn2"]] <- c(nn[["nn1"]], angles, energy)
nn[["nn3"]] <- c(nn[["nn2"]], AAs)


nn <- tibble(neural_net = names(nn),
             parameters = nn)

all_params <- unique(do.call(c, nn[["parameters"]]))





context <- c("_lag1", "_lag2", "_lead1", "_lead2")

nn2 <- nn %>%
  mutate(parameters = map(parameters, .f = ~c(., do.call(c, map(., ~paste0(., context)))))) %>%
  mutate(neural_net = paste0(neural_net, "c"))

names(nn2[["parameters"]]) <- paste0(names(nn2[["parameters"]]), "c")

nn <- bind_rows(nn, nn2)

nn <- nn %>%
  mutate(parameters = map(parameters, .f = ~c(., "known")))


if(FALSE) {
  to_remove <- "^SS_.*_(lead|lag)\\d+$"

  nn3 <- nn %>%
    filter(neural_net %in% c("nn5c", "nn6c", "nn7c")) %>%
    mutate(parameters = map(parameters, .f = ~.[!str_detect(., to_remove)]))

  names(nn3$parameters) <- nn3$neural_net <- paste0(nn3$neural_net, "noSS")

  nn <- bind_rows(nn, nn3)

  nn_og <- readRDS("~/peptide_alg/nn_models_256_64_4_leaky_relu_lion_opt_large_val_better_angles_w_context.rds")
  nn <- bind_cols(nn, nn_og[,3:5])
}


secretome_aa <- secretome_aa %>%
  mutate(SS = replace_na(SS, "-")) %>%
  mutate(one_hot_encode(AA, tag = "AA")) %>%
  mutate(one_hot_encode(SS, tag = "SS")) %>%
  mutate(one_hot_encode(gpcrdb_gtp_type, tag = "known")) %>%
  mutate(one_hot_encode(min_AA, tag = "minAA")) %>%
  mutate(one_hot_encode(max_AA, tag = "maxAA")) %>%
  mutate(across(c("Phi", "Psi"), .fns = ~tibble(cos = cos(. * pi/180),
                                                sin = sin(. * pi/180)), .unpack = TRUE)) %>%
  dplyr::group_by(accession) %>%
  mutate(across(ends_with("_energy"), .fns = \(x) {
    x[is.na(x)] <- 0
    scales::rescale(x, to = c(0, 1))})) %>%
  ungroup() %>%
  select(!ends_with("_unknown")) %>%
  mutate(known = rowSums(across(starts_with("known_gpcr_pep"))) > 0) %>%
  mutate(across(all_of(cons_metrics), ~ .x * species_limit, .names = "{.col}_n"))


