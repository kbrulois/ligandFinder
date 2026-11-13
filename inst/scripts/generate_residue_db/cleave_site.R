


peps <- data.table::fread("~/AF2_analysis/TermDBX_CPs_256_v3_Jun23_4AF2screen.csv") %>% as_tibble()


peps <- peps %>%
  mutate(stringr::str_extract(`Peptide ID`, "x\\d+x\\d+$") %>%
           stringr::str_remove(., "^x") %>%
           ligandFinder::tibble_split(., "x", names = c("start", "end")) %>%
           mutate(across(everything(), as.numeric)))


peps <- left_join(peps %>% dplyr::rename(gene = Gene), secretome %>% select(gene, accession, sequence_uni, sp_ind), by = "gene")

scores <- list()
wins <- c("C", "N")

for(win in wins) {

peps2 <- peps %>%
  rowwise() %>%
  mutate(windows_C = pmap(list(sequence_uni = sequence_uni,
                             start = start,
                             end = end,
                             target = win,
                             sp_ind = sp_ind,
                             pep_id = `Peptide ID`),
                        .f = get_adj_db_sites))


peps3 <- peps2 %>%
  rowwise() %>%
  reframe(get_pep_data(window = windows_C, p_id = accession, gene = gene),
          gene = gene)


c_dat <- bind_cols(peps[,4],peps3)

c_dat <- c_dat %>%
  filter(!is.na(peps)) %>%
  #filter(map_lgl(dat, filter_window_dat)) %>%
  filter(map_lgl(dat, ~nrow(.) == 36))


sep_meta_dat <- function(data, nn_params = all_params3) {
  list(
    meta_data = select(data, -any_of(nn_params)),
    data = select(data, any_of(nn_params))
  )
}

c_dat <- c_dat %>%
  mutate(split = map(dat, sep_meta_dat)) %>%
  select(-dat) %>%
  unnest_wider(split) %>%
  mutate(known = 0)

c_dat <- c_dat %>%
  mutate(known_idx = list(c(rep("none", 36))))

c_dat <- c_dat %>%
  mutate(data = map(data, \(x) { x %>% mutate(across(everything(), ~replace_na(., 0)))}))

nn_in_C <- generate_keras_input(list(all = c_dat))

tmp2 <- models[[win]] |> predict(nn_in_C[["all"]][["x"]])

tmp <- data.frame(tmp2[["global"]])

tmp$`Peptide ID` <- c_dat$`Peptide ID`
tmp[[paste0("win_type_", win)]] <- c_dat$win_type

colnames(tmp)[1] <- paste0("score_", win)

scores[[win]] <- tmp

}

for(x in names(scores)) {

  peps <- left_join(peps, scores[[x]], by = "Peptide ID")

}

peps <- peps %>%
  mutate(pred_end = if_else(score_C > score_N, "C", "N")) %>%
  rowwise() %>%
  mutate(max_score = max(c_across(c(score_C, score_N)), na.rm = TRUE)) %>%
  ungroup()


data.table::fwrite(peps, file = "~/AF2_analysis/TermDBX_CPs_256_v3_Jun23_4AF2screen_KB.csv")


c_dat2 <- c_dat2 %>% arrange(target)

nn_in_C <- generate_keras_input(list(all = c_dat2 %>% filter(target == "C")))
nn_in_N <- generate_keras_input(list(all = c_dat2 %>% filter(target == "N")))

scores_N <- models$N |> predict(nn_in_N[["all"]][["x"]])
scores_C <- models$C |> predict(nn_in_C[["all"]][["x"]])

val_pred_comb <- c(scores_N[["global"]][, 1], scores_C[["global"]][, 1])

nn_input_comb <- c_dat2

nn_input_comb$pred <- val_pred_comb

nn_input_comb$per_index <- c(
  extract_per_index(scores_N[["per_index_cat"]], names(classes)),
  extract_per_index(scores_C[["per_index_cat"]], names(classes))
)



pred_to_plot <- nn_input_comb


