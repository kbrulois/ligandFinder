


bm <- data.table::fread("~/Desktop/most_recent/bm_update_3_subset_lig_features_coexpression_lig_type_clusters_receptor_features_depcod_new_metrics.csv")

bm2 <- data.table::fread('~/Desktop/most_recent/bm_contact_score_latest3.csv')


common_cols <- intersect(setdiff(names(bm), "code"), names(bm2))

# Drop them from right-hand side before the join
bm2_clean <- bm2 %>% select(-all_of(common_cols))

# Perform left join
result <- left_join(bm2_clean, bm, by = "code")









dat <- data.table::fread("~/Desktop/most_recent/shap_final2.csv")

dat <- as_tibble(dat)






names(all_colors$subset.colors) <- unique(dat$feature)

dat <- dat %>%
  mutate(feature = factor(feature, levels = rev(c("paeL", "paeR", "pLDDT_rec", "pLDDT_lig1", "frequency_scaled_lig1",
                                                  "frequency_scaled_rec", "mean_af_missense_rec", "mean_af_missense_lig1",
                                                  "favorability", "CP", "sb", "ds", "area_scaled")))) %>%
  filter(ligand_index != "") %>%
  mutate(ligand_index = factor(ligand_index, levels = c("L1-5", "L6-10", "L11-20", "L21-50", "L51-277")))

bm_dat <- data.table::fread("~/Desktop/most_recent/bm_contact_score_latest2.csv")

bm_dat <- bm_dat %>%
  mutate(ligand_type = if_else(lig1_NorC_position == 1, lig1_NorC, "loop"))

dat <- left_join(dat, bm_dat %>% select(code, ligand_type), by = "code")



to_plot <- dat %>%
  group_by(shap_name, ligand_index) %>%
  summarise(shap_value = sum(shap_value))

to_plot <- left_join(to_plot, dat %>%
                       select(shap_name, known_pair, gpcr_family, data_split_6, ligand_index, feature) %>%
                       distinct(ligand_index, shap_name, .keep_all = TRUE), by = join_by(shap_name, ligand_index))

p <- ggplot2::ggplot(data = to_plot,
                       mapping = ggplot2::aes(x = ligand_index, y = shap_value)) +

    ggplot2::geom_violin(
      scale = "width",
      #position = ggplot2::position_dodge(0.9),
      trim = T,
      linewidth = 0.1,
      width = 0.5,
      alpha = 0.5,
      show.legend = F) +

    #ggplot2::facet_wrap(facets = vars(ligand_index), nrow = 2) +

    ggbeeswarm::geom_quasirandom(mapping = aes(color = known_pair,
                                             fill = known_pair),
                               pch = 21,
                               size = 0.7,
                               width = 0.5,
                               stroke = 0.1,
                               inherit.aes = TRUE,
                               show.legend = TRUE) +

  ggplot2::scale_fill_manual(values = c(`known` = "#709AE1FF", `unknown` = "#FED43980"),
                               name = "known receptor\nligand pair",
                               na.value = "white",
                               guide = ggplot2::guide_legend(keyheight = 1, override.aes = list(size = 3,
                                                                                                alpha = 0.6, stroke = 0.2,
                                                                                                color = c(`known` = "#709AE1FF", `unknown` = "#FED439FF"),
                                                                                                fill = c(`known` = "#709AE1FF", `unknown` = "#FED439FF"),
                                                                                                pch = 21, linetype = 0))) +

    ggplot2::scale_color_manual(values = c(`known` = "#709AE1FF", `unknown` = "#FED439FF"),
                                name = "known receptor\nligand pair",
                                na.value = "white",
                                guide = NULL) +

    ggplot2::geom_hline(yintercept = 0) +

    ylab("shap value") +

    ggtitle("SHAP analysis of ligand residues") +

    ggplot2::theme(
      #strip.text = ggplot2::element_text(face = "bold", margin = margin(t = 3, b = 3)),
      #strip.background.x = elementalist::element_rect_round(radius = unit(4, "pt")),
      #strip.background.y = elementalist::element_rect_round(inherit.blank = TRUE, radius = unit(4, "pt")),
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.position = "right",
      strip.switch.pad.grid = ggplot2::unit(0.2, "cm"),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(),
      axis.title.y = ggplot2::element_text(size = 18, face = 2),
      axis.text.y = ggplot2::element_text(size = 12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.grid.major = element_line(colour = "grey94"),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor.x = element_blank(),
      panel.grid.minor = element_line(linewidth = rel(0.2)),
      legend.text = ggplot2::element_text(size = 10, face = 2),
      legend.title = ggplot2::element_text(size = 10, face = 2),
      legend.title.position = "top",
      #panel.spacing=ggplot2::unit(0.2,"cm"),
      plot.title = ggplot2::element_text(size = 14, face = 2, hjust = 0.5)
    )

  svglite::svglite(filename = "~/Desktop/most_recent/new/ligand_residues.svg", width = 7, height = 5)
  print(p)

  dev.off()


  dat2 <- dat %>%
    group_by(shap_name, bw_index) %>%
    summarise(shap_value = sum(shap_value, na.rm = TRUE))

  bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

  dat2 <- left_join(dat2, bw_align %>% rename(bw_index = BW), by = "bw_index")

  dat2 <- left_join(dat2, dat %>%
                      select(shap_name, known_pair, gpcr_family, data_split_6, ligand_index, feature, bw_index, ligand_type) %>%
                      distinct(bw_index, shap_name, .keep_all = TRUE), by = join_by(shap_name, bw_index))

  dat2 <- dat2 %>%
          filter(bw_index != "")


dat2 <-  dat2 %>%
    filter(bw_index != "") %>%
    mutate(CP = if_else(CP == "CP", "CP", if_else(is.na(CP), "other", "other")))
dat2$CP[is.na(dat2$CP)] <- "other"

bw_align2 <- c(paste0("N", 18:1), bw_align$BW[c(1:120, 225:227, 121:212)])

dat2 <- dat2 %>%
          filter(bw_index %in% bw_align2) %>%
          mutate(bw_index = factor(bw_index, levels = bw_align2)) %>%
          mutate(region2 = as.character(bw_index) %>% str_extract(., "^.")) %>%
          mutate(region2 = if_else(region2 == "N", "N-term", paste0("TM", region2)))

uni_bw <- tibble(bw_index = levels(dat2[["bw_index"]]))

uni_bw <- left_join(uni_bw, dat2 %>% ungroup %>% distinct(bw_index, .keep_all = TRUE)  %>% select(bw_index, region2), by = "bw_index")

dat3 <- factor_to_uniprotFeature(uni_bw$region2)

label_dat <- dat2 %>%
              filter(shap_value > 0) %>%
              pull(bw_index) %>%
              as.character %>%
              unique

cp_bw <- c("1.28", "1.35", "1.39", "2.53", "2.56", "2.57", "2.59", "2.60", "2.61", "2.63", "2.64", "3.21", "3.26", "3.28", "3.29", "3.30", "3.32", "3.33", "3.34", "3.36", "3.37", "3.40", "4.56", "4.57", "4.58", "4.59", "4.60", "4.61", "45.51", "45.52", "5.33", "5.35", "5.36", "5.37", "5.39", "5.40", "5.43", "5.44", "5.46", "5.461", "5.47", "6.44", "6.48", "6.51", "6.52", "6.54", "6.55", "6.58", "6.59", "6.61", "6.62", "7.24", "7.27", "7.30", "7.31", "7.33", "7.34", "7.35", "7.37", "7.38", "7.39", "7.41", "7.42")
both <- c("1.35", "1.39", "2.53", "2.57", "2.61", "2.64", "3.28", "3.29", "3.32", "3.33", "3.36", "3.37", "3.40", "4.56", "4.57", "4.60", "4.61", "5.39", "5.43", "5.46", "5.461", "6.44", "6.48", "6.51", "6.52", "6.55", "6.58", "6.59", "7.34", "7.38", "7.41", "7.42")

label_dat <- uni_bw %>%
              mutate(x_pos = 1:nrow(.)) %>%
              filter(bw_index %in% label_dat) %>%
              mutate(y_pos = rep(c(-1.4, -1.6, -1.8, -2, -2.2, -2.4), 20)[1:nrow(.)]) %>%
              mutate(known = case_when(bw_index %in% both ~ "Gloriam et al.\n& Ngo et al.",
                                       bw_index %in% cp_bw ~ "Ngo et al.",
                                       TRUE ~ "previously\nundescribed"))

dat2 <- dat2 %>%
          mutate(gpcr_family = stringr::str_remove(gpcr_family, "\\sreceptors$"))

uni_gpcr_fam <- as.character(unique(dat2$gpcr_family))

uni_gpcr_fam <- uni_gpcr_fam[gtools::mixedorder(uni_gpcr_fam)]

gpcr_fam_cols <- rep(c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF",
                   "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF",
                   "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF",
                   "#FD8CC1FF"), 2)[1:length(uni_gpcr_fam)]

names(gpcr_fam_cols) <- uni_gpcr_fam

color_scale <- list(`known receptor\nligand pair` = c(`known` = "#709AE1FF", `unknown` = "#FED43980"),
                  gpcr_family = gpcr_fam_cols)

col_toplot <- "gpcr_family"

dat2 <- dat2 %>%
  mutate(gpcr_family = factor(gpcr_family, levels = names(color_scale$gpcr_family)))

dat2 <- dat2 %>%
        mutate(chemokine = if_else(gpcr_family == "Chemokine", "chemokine", "non-chemokine"))



  library(ggtree)

test <- dat2 %>%
  filter(known_pair == "known" & shap_value >= 0) %>%
  group_by(gpcr_family, shap_name, bw_index) %>%
  summarize(shap_value = sum(shap_value, na.rm = TRUE))

test <- test %>%
  pivot_wider(names_from = c("shap_name", "bw_index"), values_from = "shap_value", values_fill = 0)



D <- proxy::dist(test %>% ungroup %>% select(-gpcr_family) %>% as.matrix, method = "euclidean")
tre <- ape::nj(D)
tre <- ape::ladderize(tre)

tre$tip.label <- as.character(test[["gpcr_family"]])


dist_matrix <- cophenetic(tre)
hc <- hclust(as.dist(dist_matrix), method = "average")
groups <- cutree(hc, k = 8)
groups

groups <- paste0("group", groups)



tree_plot <- ggtree::ggtree(tre, layout = "rectangular") +
            ggtree::geom_tiplab(aes(label = label), hjust = 0)

tree_plot <- tree_plot + xlim(0, max(tree_plot$data$x) + 10)

tree_plot

svglite::svglite(filename = "~/Desktop/most_recent/new/gpcr_family_shap_tree.svg", width = 5, height = 7)
print(tree_plot)

dev.off()

tree_plot <- tree_plot %<+% data.frame(node = 1:length(groups), group = as.character(unname(groups))) +
                aes(color = group) +
            scale_color_discrete()

print(tree_plot)



col_toplot <- names(color_scale)[1]

gpcr_fam_sub <- c("Chemokine", "Bradykinin", "Neuropeptide FF/neuropeptide AF", "Cholecystokinin", "Class A orphans", "Galanin", "Melanocortin")

shape_scale <- c(`N` = 21, `C` = 22, `loop` = 23)

residue_scale <- c("#FD7446FF", "#1A9993FF", "black")

names(residue_scale) <- c("Gloriam et al.\n& Ngo et al.", "Ngo et al.", "previously\nundescribed")



p <- ggplot2::ggplot(data = dat2 %>%
                       filter(!is.na(ligand_type)) %>%
                       mutate(ligand_type = factor(ligand_type, levels = c("N", "C", "loop"))),
                     mapping = ggplot2::aes(x = bw_index, y = shap_value)) +


  ggrastr::rasterise(ggbeeswarm::geom_quasirandom(mapping = aes(x = bw_index,
                                                             y = shap_value,
                                                             color = known_pair,
                                                             shape = ligand_type),
                                                  fill = "#FFFFFF00",
                                               corral.width = 1,
                                               corral = "random",
                                               size = 0.7, cex = 0.01,
                                               stroke = 0.1,
                                               inherit.aes = TRUE,
                                               show.legend = TRUE), dpi = 600) +

  ggplot2::geom_violin(mapping = ggplot2::aes(x = bw_index, y = shap_value),
                       scale = "width",
                       fill = NA,
                       #position = ggplot2::position_dodge(1),
                       trim = T,
                       linewidth = 0.1,
                       width = 1,
                       alpha = 1,
                       show.legend = F) +

  ggplot2::scale_color_manual(values = color_scale[[col_toplot]],
                              name = "receptor-ligand\npair",
                              na.value = "black",
                              guide = ggplot2::guide_legend(keyheight = 1, override.aes = list(size = 3, stroke = 0.5, pch = 21, fill = "#FFFFFF00"))) +

  ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_segment(data = dat3,
                                                                 aes(x = start,
                                                                     xend = end,
                                                                     y = -1.2),
                                                                 color = "white",
                                                                 linewidth = 4, inherit.aes = FALSE),
                                           colour = "black",
                                           sigma = 1,
                                           expand = 5),
                     dev = "ragg",
                     dpi = 600) +

  ggplot2::geom_text(data = dat3 %>% rowwise %>% mutate(center = mean(c_across(all_of(c('start', 'end'))))),
                     aes(x = center, y = -1.2, label = type), inherit.aes = FALSE) +


  #ggplot2::facet_grid(rows = vars(gpcr_family), switch = "y") +


  ggplot2::scale_shape_manual(values = unname(shape_scale),
                              name = "ligand insertion",
                              guide = ggplot2::guide_legend(keyheight = 1, override.aes = list(size = 3, stroke = 0.5))) +



  ggnewscale::new_scale_color() +

  ggplot2::scale_color_manual(values = residue_scale,
                              name = "involvement in\nligand binding",
                              na.value = "black", guide = ggplot2::guide_legend(keyheight = 1, override.aes = list(size = 3, stroke = 0.5))) +

  ggplot2::geom_label(data = label_dat, aes(x = bw_index, y = y_pos, label = bw_index, color = known), size = 2) +



  ggplot2::geom_hline(yintercept = 0) +

  ylab("SHAP value") +

  ggtitle("Receptor residues") +

  ggplot2::theme(
    #strip.text = ggplot2::element_text(face = "bold", margin = margin(t = 3, b = 3)),
    #strip.background.x = elementalist::element_rect_round(radius = unit(4, "pt")),
    #strip.background.y = elementalist::element_rect_round(inherit.blank = TRUE, radius = unit(4, "pt")),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.position = "right",
    strip.switch.pad.grid = ggplot2::unit(0.2, "cm"),
    axis.ticks = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    #axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.text.x = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_line(),
    axis.title.y = ggplot2::element_text(size = 14, face = 2),
    axis.text.y = ggplot2::element_text(size = 12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = NA),
    panel.grid.major = element_line(colour = "grey94"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor = element_line(linewidth = rel(0.2)),
    legend.text = ggplot2::element_text(size = 10, face = 2),
    legend.title = ggplot2::element_text(size = 10, face = 2),
    legend.title.position = "top",
    #panel.spacing=ggplot2::unit(0.2,"cm"),
    plot.title = ggplot2::element_text(size = 14, face = 2, hjust = 0.5)
  )


svglite::svglite(filename = "~/Desktop/most_recent/new/receptor_residues.svg", width = 18, height = 5)
print(p)

dev.off()





