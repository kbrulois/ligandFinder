



score_var <- "score_nn8_s8_entire"

scores <- secretome_roi %>% 
  filter(!grepl("^phs_", roi_name)) %>% 
  pull(!!sym(score_var))

cats <- quantile(scores, c(0.25, 0.4, 0.8), na.rm = TRUE)

names(cats) <- paste0("top ", 100 - as.numeric(sub("%$", "", names(cats))), "%\n(score > ", round(unname(cats), 2), ")")

cats <- c(setNames(0, "all"), cats)

cat_colors <- setNames(c("#8A9197FF", "#FD7446FF", "#FED439FF", "#709AE1FF"),
                       names(cats))


dat <- bind_rows(
  lapply(names(cats), \(x)
         secretome_roi %>%
           filter(!!sym(score_var) > cats[[x]]) %>%
           summarize(AA_known = sum(roi_length[roi_type == "gtp"], na.rm = TRUE),
                     AA_putative = sum(roi_length[roi_type != "gtp"], na.rm = TRUE),
                     peptides_known = sum(roi_type == "gtp", na.rm = TRUE),
                     peptides_putative = sum(roi_type != "gtp", na.rm = TRUE),
                     cutoff = cats[[x]],
                     color = x) 
  ))

dat$color <- factor(dat$color, levels = names(cats))

dat <- bind_rows(
  dat %>%
    dplyr::rename(known = peptides_known, putative = peptides_putative) %>%
    select(known, putative, cutoff, color) %>%
    mutate(type = "peptide"),
  dat %>%
    dplyr::rename(known = AA_known, putative = AA_putative) %>%
    select(known, putative, cutoff, color) %>%
    mutate(type = "AA")
)


library(ggforce)

d <- dat %>%
  pivot_longer(cols = c('known', 'putative'), names_to = 'known') %>%
  group_by(type, known) %>%
  mutate(value_resc = value/max(value)) %>%
  ungroup

d$r <- sqrt(d$value_resc / pi)
d$y0 <- max(d$r) / 2
d$x0 <- max(d$r) / 2


d <- d %>%
  group_by(type, known) %>%
  mutate(x0 = x0 - 0.5 * r) %>%
  mutate(y0 = y0 + 0.5 * r) %>%
  mutate(x_label = max(x0)) %>%
  mutate(y_edges = sqrt(r^2 - (x0 - x_label)^2) + y0) %>%
  arrange(y_edges) %>%
  mutate(y_label = ifelse(row_number() == 1, 
                          y_edges - r,
                          (y_edges + lag(y_edges)) /2)) %>%
  arrange(desc(y_edges)) %>%
  mutate(percent = 100 * value/max(value)) %>%
  dplyr::rename(score = cutoff) %>%
  ungroup %>%
  arrange(score)

venn <- ggplot(d, aes(x0 = x0, y0 = y0, r = r, fill = color)) +
  geom_circle() +
  geom_label(aes(x = x_label, 
                 y = y_label, 
                 label = paste0(value, " (", round(percent, 0), "%)")), inherit.aes = FALSE) +
  facet_grid(rows = vars(type), cols = vars(known), switch = "y") +
  scale_fill_manual(values = cat_colors, name = "Score Cutoff") +
  theme_void() +  theme(legend.key.height = unit(2, "lines"), strip.text.y = element_text(angle = 90)
  )






v_lines <- d %>%
  filter(type == "peptide" & known == "known") %>%
  dplyr::slice(-1)

v_lines2 <- d %>%
  filter(type == "peptide" & known == "putative") %>%
  dplyr::slice(-1)


violin <- ggplot(secretome_roi %>% arrange(!is.na(percent_ol_gtp)), 
                 mapping = aes(x = !!sym(score_var), y = roi_type, color = percent_ol_gtp)) +
  ggbeeswarm::geom_quasirandom(fill = "#80808000", position = ggplot2::position_dodge(0.9), size = 2, stroke = 0.8, pch = 21, orientation = "y", na.rm = TRUE) +
  scale_color_viridis_c(option = "H") +
  ggnewscale::new_scale_color() + 
  geom_vline(data = v_lines,  
             mapping = aes(xintercept = score, color = color), 
             lwd = 2, 
             lty = "dashed") + 
  ggplot2::annotate(geom = "label", x = v_lines[["score"]], y = 0.5, label = paste0(round(v_lines[["percent"]], 0), "%"), angle = 0, hjust = 0.5) +
  ggplot2::annotate(geom = "text", x = min(v_lines[["score"]]) - 0.05, y = 0.5, label = "percent of known: ", hjust = 1) +
  ggplot2::annotate(geom = "label", x = v_lines[["score"]], y = 2.5, label = paste0(round(v_lines2[["percent"]], 0), "%"), angle = 0, hjust = 0.5) +
  ggplot2::annotate(geom = "text", x = min(v_lines[["score"]]) - 0.05, y = 2.5, label = "percent of putative: ", hjust = 1) +
  scale_y_discrete(na.translate = FALSE) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = cat_colors[-1], name = "Score Cutoff)", na.value = "#80808000") +
  theme_bw()+
  theme(legend.key.height = unit(2, "lines"),
        panel.border = element_blank()
  ) +
  ylab("")




design <- "
12
#2
"

p_al <- violin + venn +
  
  patchwork::plot_layout(design = design)


ggsave("~/peptide_alg/gating5.svg", plot = p_al, device = svglite::svglite, width = 16, height = 6.75)












x_var <- sym("roi_cons")
y_var <- sym("roi_ASA")

color_var <- sym("score")

p <- ggplot2::ggplot(secretome_roi %>% drop_na(known_ligand3) %>% arrange(desc(is.na(!!color_var)),!!color_var)) + 
  ggrastr::rasterise(ggplot2::geom_point(aes(x = !!x_var, y = !!y_var),
                                         size = 1,
                                         pch = 1),
                     dev = "ragg",
                     dpi = 100) +
  xlab("conservation") +
  ylab("relASA")

p <- p +
  scale_color_viridis_c(option = "B",  guide = guide_colorbar(
    barwidth = 1, barheight = 8, frame.colour = "black", ticks.linewidth = 0.7, ticks.colour = "black", frame.linewidth = 0.5
  )) +
  ggnewscale::new_scale_color() + 
  facet_wrap(facets = "known_ligand3", ncol = 2, drop = TRUE)

p <- p + ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_density_2d(mapping = ggplot2::aes(x = !!x_var,
                                                                                                  y = !!y_var),
                                                                           color = "black",
                                                                           show.legend = FALSE),
                                                  colour = "white", sigma = 2, expand = 3),
                            dev = "ragg",
                            dpi = 600)



for(i in 2:length(cats)) {
  gate_dat <- secretome_roi %>% filter(score > cats[i])
  hull_indices <- chull(gate_dat[[x_var]], gate_dat[[y_var]]) 
  hull_points <- gate_dat[hull_indices, ]    
  hull_points <- rbind(hull_points, hull_points[1, ])
  
  hull_points2 <- hull_points
  hull_points2$known_ligand3 <- "known"
  hull_points$known_ligand3 <- "putative"
  hull_points <- bind_rows(hull_points, hull_points2)
  
  hull_points$color <- names(cats)[i]
  
  p <- p +                       
    ggrastr::rasterise(ggfx::with_outer_glow(geom_polygon(data = hull_points,     
                                                          aes(x = !!x_var, y = !!y_var, color = color),
                                                          fill = "#33333300", lwd = 1, lty = "dashed"),
                                             colour = "white", sigma = 2, expand = 3),
                       dev = "ragg",
                       dpi = 600) +
    scale_color_manual(values = cat_colors[-1], name = "Score Cutoff", na.value = "#80808000") +
    theme(legend.key.height = unit(6, "lines")  
    )
  
}

p <- p +  theme_bw() + theme(legend.position = "none",
                             legend.key.height = unit(6, "lines")  
)





design <- "
111333
222333
"

p_al <- violin + p + venn +
  
  patchwork::plot_layout(design = design)


ggsave("~/peptide_alg/gating2.svg", plot = p_al, device = svglite::svglite, width = 16, height = 8)


