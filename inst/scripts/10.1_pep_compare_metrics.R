



cols_to_sum <- c("cons_rs", "cons_lrs", "blos_wt_all",
                 "blos_uw_all", "blos_wt_mam", "blos_uw_mam", "gran_wt_all", "gran_uw_all",
                 "gran_wt_mam", "gran_uw_mam", "blos_nr_all", "blos_nr_mam", "gran_nr_all",
                 "gran_nr_mam")

cols_to_sum2 <- c(cols_to_sum, paste0(cols_to_sum, "_n"), c("max_afm", "min_afm", "mean_afm", "relASA"))

uniprot_sub <- secretome_aa %>%
  filter(is.na(`uniprot_signal peptide_type`)) %>%
  filter(!is.na(cons_og)) %>%
  filter(!gene %in% known_genes[["chem"]]) %>%
  #filter(location %in% paste0(1:4, "l")) %>%
  #filter(accession %in% known_peps$accession) %>%
  filter(relASA > 0.5 & SS != "E") %>%
  mutate(across(all_of(cols_to_sum), ~ .x * species_limit, .names = "{.col}_n"))

mean_dif <- uniprot_sub %>%
  group_by(known) %>%
  summarise(across(all_of(cols_to_sum2), ~mean(., na.rm = TRUE))) %>%
  summarise(across(all_of(cols_to_sum2), ~ (.x[2] - .x[1])))

sd_res <- uniprot_sub %>%
  summarise(across(all_of(cols_to_sum2), ~sd(., na.rm = TRUE)))


svglite::svglite(filename = "~/AF2_analysis/metric_comparison3.svg", width = 6, height = 6)
ggplot2::ggplot((mean_dif/sd_res) %>%
                  pivot_longer(everything()) %>%
                  mutate(col = if_else(stringr::str_detect(name, "_n$"), "myo_normalized", "not_normalized"))) +
  geom_col(aes(x = value, y = name, fill = col), color = "black", width = 0.8) +
  xlab("Cohen's D") +
  ylab("") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  ggtitle("residues of known ligands versus other residues\n(gated on relASA > 0.5 and non-strand non-sp)") +
  theme_bw()
dev.off()










ggplot2::ggplot(uniprot_sub %>% mutate(known = ifelse(known == 1, "known", "unknown"))) +
  ggplot2::geom_boxplot(aes(x = known, y = `cons_lrs`))
















ggplot2::ggplot(uniprot_t_e %>% filter(relASA > 0.5)) +
  ggplot2::geom_violin(mapping = ggplot2::aes(x = known, y = cons_rs))



uniprot_t_e %>%
  filter(relASA > 0.5 & SS != "E") %>%
  mutate(across(all_of(cols_to_sum), ~ .x * species_limit, .names = "{.col}_n")) %>%
  group_by(known) %>%
  summarise(across(all_of(cols_to_sum2), ~mean(., na.rm = TRUE))) %>%
  View






uniprot_t_e %>%
  nest(data = -known) %>%
  summarise(
    tidy(t.test(value ~ group, data = df))
  )


uniprot_t_e %>%
  mutate(blos_wt_all_n = blos_wt_all * species_limit) %>%
  mutate(index = row_number(), .by = accession) %>%
  select(accession, AA, index, blos_wt_mam, blos_wt_all_n) %>%
  {data.table::fwrite(., "~/AF2_analysis/conservation_metrics.csv")}
