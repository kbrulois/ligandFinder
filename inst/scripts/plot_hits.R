

library(tidyverse)


files <- c("~/AF2_analysis/bm_brinp_68W.csv", "~/AF2_analysis/CXCL14peptidesOneQuarter.csv")
files <- c("~/AF2_analysis/top200NC_Nov17.csv", "~/AF2_analysis/all_metrics_oct17.csv")

files <- c("~/AF2_analysis/new_peps_scg_tm_focused.csv")


metrics <- c("iptm", "paeL_mean_in", "paeL_mean_all", "paeL_mean_out")

dat <- map(files, ~as_tibble(data.table::fread(.))) %>%
        bind_rows


dat <- dat %>%
        mutate(run_name = if_else(run_name == "bm_sep28", "bm", run_name))



dat <- dat %>%
  mutate(lig1_end2 = case_when(lig1_end == "1C" ~ "C-terminal\ninsertion",
                               lig1_end == "1N" ~ "N-terminal\ninsertion",
                               TRUE ~ "loop\ninsertion")) %>%
  mutate(location = case_when(depth < 0.2 & depth > -0.5 & radius < 1 ~ lig1_end2,
                            TRUE ~ "not inserted")) %>%
  mutate(across(all_of(metrics), ~if_else(location == "not inserted", 0, .), .names = "{.col}_loc")) %>%
  mutate(pep = paste0(p2_name, "_", p2_range)) %>%
  mutate(depth2 = if_else(location != "not inserted", depth, NA))






####peptide versus receptors

dat$pep %>% unique() %>% {.[gtools::mixedorder(.)]}

peps <- c("BRNP3_369x380", "BRNP2_386x397", "BRNP2_372x397", "BRNP2_347x397", "BRNP1_356x368", "PYY_31x64", "PYY_29x64")
peps <- c("CXL17_22x119", "CXL17_64x119", "CXL11_22x94", "SG3A2_70x82", "SG1C1_24x33", "SG1C1_86x95", "SG1D1_75x90", "SG2A1_19x28", "SG3A1_78x87",
          "SG3A1_95x104")
peps <- c("CXL14_35x102")

peps <- c("GHRL_24x51", "GHRL_76x98" )

peps <- tibble(end = c(95:111),
                        uniprot_name = "CXL14",
                        start = 88) %>%
                  mutate(peps = paste0(uniprot_name, "_", start, "x", end)) %>%
                  pull(peps)

peps <- dat$pep %>% unique() %>% {.[gtools::mixedorder(.)]}



pq_path <- "~/ligandFinder_data/residue_db"

res_db <- arrow::open_dataset(source = pq_path)

pep_uniprot_genes <- peps %>% stringr::str_split(., "_", simplify = TRUE)
pep_uniprot_genes <- pep_uniprot_genes[ ,1] %>% unique


residue_data <- res_db %>%
  filter(uni_gene %in% pep_uniprot_genes) %>%
  select(uni_gene, sequence_uni) %>%
  collect()




file_name <- "~/AF2_analysis/new_peps_scg_tm_focused_paeL.svg"

dat2 <- dat %>%
  #filter(p1_name %in% !!common_receptors) %>%
  filter(pep %in% !!peps) %>%
  mutate(tibble_split(p2_range, "x", names = c("start", "end"))) %>%
  {left_join(., residue_data, by = join_by(p2_name == uni_gene))} %>%
  mutate(legacy_ind = map2_chr(end, sequence_uni, \(x, y) {
    paste0((as.numeric(x) - 34), stringr::str_sub(y, x, x))
  }), .after = "afpd_dir_name") %>%
  mutate(pep = paste0(pep, "(", legacy_ind, ")"))


uni_peps <- dat2$pep %>% unique

dat2 <- dat2 %>%
  mutate(pep = factor(pep, levels = uni_peps[gtools::mixedorder(uni_peps)]))






dat2 <- dat2 %>%
  group_by(p1_name, pep) %>%
  arrange(rank, .by_group = TRUE) %>%
  mutate(id = row_number())



p <- ggplot(data = dat2, aes(x = id, y = paeL_mean_in_loc, fill = depth2, shape = location)) +
  ggplot2::geom_point(size = 1.5, stroke = 0.1) +
  #ggiraph::geom_point_interactive(aes(tooltip = p1_name), size = 0.6) +
  ggplot2::facet_grid(rows = vars(pep), cols = vars(p1_name), switch = "y") +
  ggplot2::scale_shape_manual(values = c(`C-terminal\ninsertion` = 25,
                                         `N-terminal\ninsertion` = 24,
                                         `loop\ninsertion` = 21,
                                         `not inserted` = 4),
                              name = "ligand\norientation") +
  ggplot2::scale_fill_gradientn(limits = c(max(dat2 %>%
                                                  filter(location != "not inserted") %>%
                                                  pull(depth)),
                                            min(dat2 %>%
                                                  filter(location != "not inserted") %>%
                                                  pull(depth))),colors = viridis::turbo(n = 20),
                                 name = "insertion\ndepth", na.value = "black") +
  xlab("") +
  ylab("") +
  ggtitle("PAEL_mean_in") +
  theme(axis.text.x = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        strip.text.x = element_text(angle = 90, hjust = 0), strip.placement = "outside",
        panel.grid.major = element_line(colour = "black", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey98", color = NA),
        panel.border = element_rect(colour = "black", linewidth = 0.1),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "lines")) +
    coord_cartesian(clip = "off")


svglite::svglite(filename = file_name, width = 40, height = length(levels(dat2$pep)))
print(p)
dev.off()




htmlwidgets::saveWidget(ggiraph::girafe(ggobj = p, width_svg = 40, height_svg = 10),
                        file = "~/Desktop/brinp_pep_rec_spec_iptm_new.html",
                        selfcontained = TRUE)



####receptors versus ligands

receptors <- c("BKRB1", "BKRB2", "AGTR1", "AGTR2", "GALR3", "RL3R1", "RL3R2", "GPR15", "GPR25")

receptors <- c("GPR15", "GPR25")


file_name <- "~/Desktop/secroglob_v_ligands_iptm.svg"

dat2 <- dat %>%
  filter(run_name == "bm" | pep %in% !!peps) %>%
  filter(p1_name %in% !!receptors) %>%
  group_by(p1_name, pep) %>%
  arrange(rank, .by_group = TRUE) %>%
  mutate(id = row_number())


p <- ggplot(data = dat2, aes(x = id, y = iptm_loc, fill = depth2, shape = location)) +
  ggplot2::geom_point(size = 1.5, stroke = 0.1) +
  #ggiraph::geom_point_interactive(aes(tooltip = p1_name), size = 0.6) +
  ggplot2::facet_grid(rows = vars(p1_name), cols = vars(pep), switch = "y") +
  ggplot2::scale_shape_manual(values = c(`C-terminal\ninsertion` = 25,
                                         `N-terminal\ninsertion` = 24,
                                         `loop\ninsertion` = 21,
                                         `not inserted` = 4),
                              name = "ligand\norientation") +
  ggplot2::scale_fill_gradientn(limits = c(max(dat2 %>%
                                                 filter(location != "not inserted") %>%
                                                 pull(depth)),
                                           min(dat2 %>%
                                                 filter(location != "not inserted") %>%
                                                 pull(depth))),colors = viridis::turbo(n = 20),
                                name = "insertion\ndepth", na.value = "black") +
  xlab("") +
  ylab("") +
  ggtitle("IPTM") +
  theme(axis.text.x = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        strip.text.x = element_text(angle = 90, hjust = 0), strip.placement = "outside",
        panel.grid.major = element_line(colour = "black", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey98", color = NA),
        panel.border = element_rect(colour = "black", linewidth = 0.1),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "lines")) +
  coord_cartesian(clip = "off")

svglite::svglite(filename = file_name, width = 30, height = length(receptors) * 2)
print(p)
dev.off()













library(tidyverse)

ligand_list <- data.table::fread(system.file("extdata/GPCRdb_known_pairings_human_plus2more_unique.csv",
                                             package = "ligandFinder")) %>% as_tibble

ligand_list <- ligand_list %>%
  mutate(afpd_dir_name = paste0(rec, "_", lig),
         known_pair = "known")

dat <- data.table::fread("~/Desktop/brinp_benchmarking_kb_oct11.csv") %>% as_tibble

dat <- dat %>%
  select(-known_pair) %>%
  {left_join(., ligand_list, by = "afpd_dir_name")} %>%
  mutate(known_pair = if_else(is.na(known_pair), "unknown", known_pair))

dat3 <- dat %>%
  filter(location == "relevant") %>%
  filter(!(run_name == "GPCRvCXCL14_oct7" & !p2_range %in% c("35x102", "35x104"))) %>%
  mutate(pep = paste0(p2_name, "_", p2_range)) %>%
  group_by(p1_name, pep) %>%
  summarise(paeL_mean_in = if_else(all(is.na(paeL_mean_in)), NA, max(paeL_mean_in, na.rm = TRUE)),
            known_pair = unique(known_pair),
            popup = paste0(unique(p1_name), unique(pep), unique(known_pair), paste0(unique(lig1_end), collapse = ","), collapse = "\n")) %>%
  mutate(paeL_mean_in = ifelse(paeL_mean_in < 0, 0, paeL_mean_in))


p <- ggplot(data = dat3, aes(x = p1_name, y = pep, fill = paeL_mean_in)) +
      #ggiraph::geom_tile_interactive(aes(tooltip = popup)) +
      ggplot2::geom_tile(color = "black") +
      ggplot2::geom_tile(data = dat3 %>% filter(known_pair == "known"), aes(x = p1_name, y = pep),
                         fill = NA,
                         color = "black",
                         linewidth = 1,
                         linejoin = "round") %>%
          ggfx::with_outer_glow(., colour = "white", sigma = 2, expand = 5) %>%
          ggrastr::rasterise(., dev = "ragg", dpi = 300) +
      scale_fill_viridis_c(option = "H", name = "") +
      scale_x_discrete(position = "top") +
      xlab("") +
      ylab("") +
      theme(axis.text.x = element_text(angle = 45),
            panel.background = element_blank(),
            axis.text.x.top = element_text(hjust = 0, vjust = 0))




svglite::svglite(filename = "~/Desktop/brinp_pep_rec_spec_paeL_mean_in.svg", width = 40, height = 32)
print(p)
dev.off()









htmlwidgets::saveWidget(ggiraph::girafe(ggobj = p, width_svg = 80, height_svg = 62),
                        file = "~/Desktop/brinp_pep_rec_spec_iptm_neww.html",
                        selfcontained = TRUE)







