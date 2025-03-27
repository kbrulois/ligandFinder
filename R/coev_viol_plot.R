


plot_coev_viol <- function(coev_res = coev_res,
                         params = c("MI", "COVinv", "COVnorm", "PWcons"),
                         plot_dir = "~/Desktop/tmp_plots",
                         pdb_names = coev_res[["pdb_names"]],
                         extra_residues = "CXCL14_W68",
                         extra_residues_name = "interacting_68W",
                         output_filename = paste0(out_path, "/viol_coev_")) {
  
  list2env(coev_res[!names(coev_res) == "pdb_names"], envir = environment())
  
  pairs_input <- clean_coev_data(co_evol_mat)
  
  coev_viol <- bind_rows(
  pairs_input %>% mutate(dPWcons = factor(rep("All", n()))),
  pairs_input)

lapply(params, \(var_tp) {
  
  p <- ggplot(data = coev_viol %>% 
                        arrange(desc(is.na(!!sym(paste0(var_tp, "_rank_top10")))), !!sym(paste0(var_tp, "_rank_top10"))) %>%
                        filter(!is.na(discDist)) %>%
                        mutate(ranks = paste0(interaction_status, pair_type, pdb_model, !!sym(paste0(var_tp, "_rank_top10")), sep = "_")), 
              mapping = aes(x = discDist, 
                            y = !!sym(var_tp), 
                            fill = discDist)) +
    geom_violin(scale = "width", fill = NA) + 
    ggbeeswarm::geom_quasirandom(aes(color = !!sym(paste0(var_tp, "_rank_top10")), 
                                     size = !!sym(paste0(var_tp, "_rank_top10"))), 
                                 fill = "#80808000", 
                                 position = ggplot2::position_dodge(0.9), 
                                 stroke = 0.8, 
                                 pch = 21, 
                                 orientation = "x", 
                                 na.rm = TRUE) +
    ggrastr::rasterise(
      ggfx::with_outer_glow(
        stat_summary(fun = mean, color = "#C80813FF", width = 0.25, fill = NA, pch = 23, outlier.shape = NA), 
        colour = "white", 
        sigma = 1, 
        expand = 1.5),
      device = "ragg", 
      dpi = 300) +
    ggrastr::rasterise(
      ggfx::with_outer_glow(
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                     geom = "errorbar", width = 0.2, color = "#C80813FF"), 
        colour = "white", 
        sigma = 1, 
        expand = 1.5),
      device = "ragg", 
      dpi = 300) +
    ggh4x::facet_nested(cols = vars(dPWcons), rows = vars(pdb_model, pair_type), scales = "free_y", independent = "y", switch = "y") +
    scale_color_manual(values = c(top10 = "#D2AF81FF"), na.value = "#197EC080") + 
    scale_size_manual(values = c(top10 = 3), na.value = 1) + 
    ggtitle("Coevolution versus residue distance, stratifed by pairwise conservation (columns) and pdb model and pairing type (rows)") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),  
          panel.grid.minor.x = element_blank())
  
  
  
  ggsave(filename = paste0(plot_dir, "/", var_tp, ".svg"), plot = p, device = svglite::svglite, width = 18, height = 24)
  
})

html_slide_show(svg_directory = plot_dir,
                output_file = paste0(output_filename, paste0(chain_names, collapse = "_"), ".html"),
                frames = vars_tp2,
                categories = NULL,
                title = paste0("coev_viol_", paste0(chain_names, collapse = "_")),
                columns = 1)

}
