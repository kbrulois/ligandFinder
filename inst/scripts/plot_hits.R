

library(tidyverse)

data.table::fread("~/AF2_analysis/brinp_v5.csv") %>%
          as_tibble %>%
          #mutate(run_name = "bm_sep28") %>%
  mutate(depth = if_else(is.na(depth), depth_percent, depth)) %>%
  mutate(radius = if_else(is.na(radius), radius_percent, radius)) %>%
  {data.table::fwrite(., "~/AF2_analysis/brinp_v5.csv")}

files <- c("~/AF2_analysis/all_metrics_oct17.csv")

metrics <- c("iptm", "paeL_mean_in")

dat <- map(files, ~as_tibble(data.table::fread(.))) %>%
        bind_rows

dat <- dat %>%
  filter(run_name %in% c("bm_sep28", "cxc17_gp15l"))

dat <- dat %>%
        mutate(run_name = if_else(run_name == "bm_sep28", "bm", run_name))


common_receptors <- c(intersect(dat %>% filter(run_name == "bm") %>% pull(p1_name) %>% unique,
                              dat %>% filter(run_name != "bm") %>% pull(p1_name) %>% unique))


peps <- c("BRNP3_369x380", "BRNP2_386x397", "BRNP2_372x397", "BRNP2_347x397", "BRNP1_356x368", "PYY_31x64", "PYY_29x64")
peps <- c("CXL17_22x119", "CXL17_64x119", "GP15L_1x81", "GP15L_25x81", "GP15L_71x81")


dat2 <- dat %>%
  mutate(lig1_end2 = case_when(lig1_end == "1C" ~ "C-terminal\ninsertion",
                               lig1_end == "1N" ~ "N-terminal\ninsertion",
                               TRUE ~ "loop\ninsertion")) %>%
mutate(location = case_when(depth < 0.2 & depth > -0.5 & radius < 1 ~ lig1_end2,
                            TRUE ~ "not inserted")) %>%
  mutate(across(all_of(metrics), ~if_else(location == "not inserted", 0, .), .names = "{.col}_loc")) %>%
  filter(run_name != "bm" | (run_name == "bm" & p2_name == "PYY")) %>%
  filter(p1_name %in% !!common_receptors) %>%
  mutate(pep = paste0(p2_name, "_", p2_range)) %>%
  group_by(p1_name, pep) %>%
  mutate(id = row_number()) %>%
  mutate(depth2 = if_else(location != "not inserted", depth, NA)) %>%
  filter(pep %in% !!peps)


dat2 <- dat %>%
  mutate(lig1_end2 = case_when(lig1_end == "1C" ~ "C-terminal\ninsertion",
                               lig1_end == "1N" ~ "N-terminal\ninsertion",
                               TRUE ~ "loop\ninsertion")) %>%
  mutate(location = case_when(depth < 0.2 & depth > -0.5 & radius < 1 ~ lig1_end2,
                              TRUE ~ "not inserted")) %>%
  mutate(across(all_of(metrics), ~if_else(location == "not inserted", 0, .), .names = "{.col}_loc")) %>%
  filter(run_name != "bm" | (run_name == "bm" & p2_name %in% c("KNG1"))) %>%
  filter(p1_name %in% !!common_receptors) %>%
  mutate(pep = paste0(p2_name, "_", p2_range)) %>%
  group_by(p1_name, pep) %>%
  mutate(id = row_number()) %>%
  mutate(depth2 = if_else(location != "not inserted", depth, NA))



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


svglite::svglite(filename = "~/Desktop/cxcl17_gp15l_lig_spec_iptm_new6.svg", width = 40, height = 5)
print(p)
dev.off()




htmlwidgets::saveWidget(ggiraph::girafe(ggobj = p, width_svg = 20, height_svg = 10),
                        file = "~/Desktop/brinp_pep_rec_spec_iptm_new.html",
                        selfcontained = TRUE)






file_name <- "~/Desktop/ACKR5_CXCL17_receptor_spec.svg"

dat2 <- dat %>%
  mutate(lig1_end2 = case_when(lig1_end == "1C" ~ "C-terminal\ninsertion",
                               lig1_end == "1N" ~ "N-terminal\ninsertion",
                               TRUE ~ "loop\ninsertion")) %>%
  mutate(location = case_when(depth < 0.2 & depth > -0.5 & radius < 1 ~ lig1_end2,
                              TRUE ~ "not inserted")) %>%
  mutate(across(all_of(metrics), ~if_else(location == "not inserted", 0, .), .names = "{.col}_loc")) %>%
  filter(p1_name %in% c("GP182", "GPR25")) %>%
  mutate(pep = paste0(p2_name, "_", p2_range)) %>%
  group_by(p1_name, pep) %>%
  mutate(id = row_number()) %>%
  mutate(depth2 = if_else(location != "not inserted", depth, NA))



p <- ggplot(data = dat2, aes(x = id, y = iptm_loc, color = depth2, shape = location)) +
  ggplot2::geom_point(size = 1) +
  #ggiraph::geom_point_interactive(aes(tooltip = p1_name), size = 0.6) +
  ggplot2::facet_grid(rows = vars(p1_name), cols = vars(pep), switch = "y") +
  ggplot2::scale_shape_manual(values = c(`relevant` = 19, `irrelevant` = 4)) +
  ggplot2::scale_color_gradientn(limits = c(max(dat2 %>%
                                                  filter(location == "relevant") %>%
                                                  pull(depth)),
                                            min(dat2 %>%
                                                  filter(location == "relevant") %>%
                                                  pull(depth))),colors = viridis::turbo(n = 20),
                                 name = "insertion\ndepth") +
  xlab("") +
  ylab("") +
  ggtitle("IPTM") +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90), strip.placement = "outside",
        panel.grid.major = element_line(colour = "black", linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.1),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "lines")) +
  coord_cartesian(clip = "off")


svglite::svglite(filename = file_name, width = 20, height = 4)
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







