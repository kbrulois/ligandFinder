



res_pairs <- readRDS("/oak/stanford/groups/ebutcher/kevin/res_pairs.rds")

pair_lut <- setNames(res_pairs[["name"]], res_pairs[["pairs"]])

runs_c <- runs_c %>%
  mutate(AA_pair = paste0(AA_rec, AA_lig1)) %>%
  mutate(AA_pair_type = if_else(AA_pair %in% res_pairs[["pairs"]], pair_lut[AA_pair], NA))


metrics <- c("paeL", "paeR")


test <- runs_m %>%
  filter(location == "relevant") %>%
  mutate(map_dfr(contacts, \(x) {
    if(!is.null(x)) {
    x %>%
      filter(in_pocket & !grepl("ICL", protein_segment) & dist < 2) %>%
      filter(mean_af_missense_rec > 0.7 & mean_af_missense_lig1 > 0.7) %>%
      summarise(across(all_of(metrics), ~mean(., na.rm = TRUE), .names = "{.col}_mean_cons"),
                across(all_of(metrics), ~min(., na.rm = TRUE), .names = "{.col}_min_cons")) %>%
      pivot_wider(names_from = "name", values_from = "value")
    }

  }))




runs_c %>%
  filter(afpd_dir_name %in% c("hNPY2R_hNPYx29x64", "hAGTR1_hANGTx25x32", "hAGRA1_hANGTx25x32", "hNPY2R_hBRNP1x356x368", "hNPY2R_hUCN2x72x109") & rank == 0) %>%
  select(!where(is.list)) %>%
data.table::fwrite(., "/scratch/groups/ebutcher/deorphan/analysis/select_cons.csv")


test <- runs_c %>%
  filter(in_pocket & !grepl("ICL", protein_segment) & dist > 1.5 & run_name == "bm_sep28" & !is.na(AA_pair_type) & dist < 4.6 & area > 6) %>%
  filter(mean_af_missense_rec > 0.7 & mean_af_missense_lig1 > 0.7) %>%
  mutate(across(all_of(c("paeL", "paeR")), ~ (30 - .)/ 30)) %>%
  group_by(code) %>%
  summarise(across(all_of(metrics), ~mean(., na.rm = TRUE), .names = "{.col}_mean_cons"),
            across(all_of(metrics), ~max(., na.rm = TRUE), .names = "{.col}_max_cons"),
            total = n())


test3 <- runs_c %>%
  filter(in_pocket & !grepl("ICL", protein_segment) & dist > 1.5 & run_name == "bm_sep28" & !is.na(AA_pair_type) & dist < 4.6 & area > 6) %>%
  mutate(wts = case_when(mean_af_missense_rec >= 0.7 & mean_af_missense_lig1 >= 0.7 ~ 1,
                         mean_af_missense_rec >= 0.7 & mean_af_missense_lig1 < 0.7 ~ 0.5,
                         mean_af_missense_rec < 0.7 & mean_af_missense_lig1 >= 0.7 ~ 0.5,
                         mean_af_missense_rec < 0.7 & mean_af_missense_lig1 < 0.7 ~ 0.1)) %>%
  #filter(mean_af_missense_rec > 0.7 & mean_af_missense_lig1 > 0.7) %>%
  mutate(across(all_of(c("paeL", "paeR")), ~ (30 - .)/ 30)) %>%
  group_by(code) %>%
  summarise(across(all_of(metrics), ~mean(., na.rm = TRUE), .names = "{.col}_mean_alll"),
            across(all_of(metrics), ~weighted.mean(x = ., w = wts, na.rm = TRUE), .names = "{.col}_mean_alll_w"),
            across(all_of(metrics), ~max(., na.rm = TRUE), .names = "{.col}_max_alll"),
            rl_cor = cor(mean_af_missense_rec, mean_af_missense_lig1),
            rl_rat = sum(mean_af_missense_rec > 0.7 & mean_af_missense_lig1 > 0.7)/n())



test2 <- left_join(runs_m  %>%
                     filter(location == "relevant" & run_name == "bm_sep28"), test, by = "code")
test2 <- left_join(test2, test3, by = "code")

sum(!is.na(test2[["paeL_mean_cons"]]))
sum(!is.na(test2[["paeL_mean_alll"]]))

test2 %>%
  mutate(cons_avail = !is.na(paeL_mean_cons)) %>%
  {table(.[["cons_avail"]], .[["known_pair"]])}

test2 <- test2 %>%
  mutate(paeL_mean_cons = if_else(is.na(paeL_mean_cons), paeL_mean_alll, paeL_mean_cons)) %>%
  mutate(paeL_max_cons = if_else(is.na(paeL_max_cons), paeL_max_alll, paeL_max_cons))


yardstick::roc_auc_vec(factor(test2$known_pair), test2$paeL_mean_cons)

yardstick::roc_auc_vec(factor(test2$known_pair), test2$rl_rat)

yardstick::roc_auc_vec(factor(test2$known_pair), test2$paeL_mean_alll)
yardstick::roc_auc_vec(factor(test2$known_pair), test2$paeL_mean_alll_w)


yardstick::roc_auc_vec(factor(test2$known_pair), test2$paeL_max_cons)

yardstick::roc_auc_vec(factor(test2$known_pair), test2$total * test2$paeL_mean_cons)


