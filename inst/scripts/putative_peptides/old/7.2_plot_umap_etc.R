



subsetter <- secretome$roi_length < 3 & !is.na(secretome$roi_length) & (is.na(secretome$dbr_terminal_type) | secretome$dbr_terminal_type == 0)

to_umap <- secretome[!subsetter, ] %>%
  dplyr::select(roi_afm, prox_afm_3, roi_cons, roi_ASA, prox_relASA_3, prox_cons_3, dbr_WY_cons, dbr_WY_asa, score) %>%
  mutate(across(contains("_cons"), .fns = \(x) scales::rescale(x, c(0,2)))) %>%
  mutate(across(everything(), .fns = ~replace_na(data = ., replace = 0))) %>%
  #mutate(across(everything(), .fns = \(x) {x[is.infinite(x)] <- 0; return(x)})) %>%
  #mutate(relative_cons = roi_cons - non_roi_cons) %>%
  #select(roi_cons, prox_relASA_6) %>%
  as.matrix

umap_config <- umap::umap.defaults

umap_config$min_dist <- 0.5
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200

umap_res <- umap::umap(d = to_umap, config = umap_config)

secretome[["UMAP1"]] <- NA
secretome[["UMAP2"]] <- NA

secretome[["UMAP1"]][!subsetter] <- umap_res[["layout"]][,1]
secretome[["UMAP2"]][!subsetter] <- umap_res[["layout"]][,2]






to_plot <- secretome[!subsetter,] %>%
  mutate(stream = str_detect(roi, "_us"))

features <- c("score",
              paste0("prox_relASA_", c(3)),
              paste0("prox_cons_", c(3)),
              "dbr_Dibasic_asa",
              "roi_ASA",
              "roi_cons",
              "roi_afm",
              "prox_afm_3",
              "dbr_WY_cons",
              "dbr_Dibasic_cons")

plot_type <- c("cons_v_asa", "umap")

tmp_dir <- "~/peptide_alg/tmp6/"
dir.create(tmp_dir)

for(y in plot_type) {
  
  if(y == "cons_v_asa") {
    x_var <- sym("roi_cons")
    y_var <- sym("roi_ASA")
  }
  if(y == "umap") {
    x_var <- sym("UMAP1")
    y_var <- sym("UMAP2")
  }
  
  for(x in features) {
    
    color_var <- sym(x)
    
    is_numeric <- is.numeric(to_plot[[x]])
    
    p <- ggplot2::ggplot(to_plot %>% arrange(desc(is.na(!!color_var)),!!color_var)) + 
      ggrastr::rasterise(ggplot2::geom_point(aes(x = !!x_var, y = !!y_var, color = !!color_var),
                                             size = 1,
                                             pch = 19),
                         dev = "ragg",
                         dpi = 100) +
      xlab(as.character(x_var)) +
      ylab(as.character(y_var))
    
    p <- p +
      scale_color_viridis_c(option = "B",  guide = guide_colorbar(
        barwidth = 1, barheight = 8, frame.colour = "black", ticks.linewidth = 0.7, ticks.colour = "black", frame.linewidth = 0.5
      )) +
      ggnewscale::new_scale_color()
    
    p <- p + ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_density_2d(mapping = ggplot2::aes(x = !!x_var,
                                                                                                      y = !!y_var),
                                                                               color = "black",
                                                                               show.legend = FALSE),
                                                      colour = "white", sigma = 2, expand = 3),
                                dev = "ragg",
                                dpi = 600)
    
    
    # to_plot2 <- to_plot %>% 
    #   filter(gene %in% to_annotate$gene & dbr_name %in% to_annotate$roi) %>%
    #   mutate(CCL_type = setNames(to_annotate$name, to_annotate$gene)[gene])
    
    to_plot2 <- to_plot %>% 
      filter(dbr_type == "gtp")
    
    
    # p <- p + ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_density_2d(data = to_plot2, mapping = ggplot2::aes(x = !!x_var,
    #                                                                                                                    y = !!y_var),
    #                                                                            color = "red",
    #                                                                            show.legend = FALSE),
    #                                                   colour = "white", sigma = 2, expand = 3),
    #                             dev = "ragg",
    #                             dpi = 600)
    
    p <- p + ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_point(data = to_plot2, 
                                                                          aes(x = !!x_var, y = !!y_var, color = score), 
                                                                          size = 2, 
                                                                          pch = 19), 
                                                      colour = "black",
                                                      sigma = 3, 
                                                      expand = 5),
                                dev = "ragg",
                                dpi = 600) +
      # scale_color_manual(values = setNames(to_annotate$color, to_annotate$name), name = NULL) + 
      scale_color_viridis_c(option = "D") +
      guides(color = guide_legend(override.aes = list(size = 6)))
    
    
    p <- p +  theme_bw()
    
    ggsave(filename = paste0(tmp_dir, x, ".svg"),
           plot = p, svglite::svglite, width = 6, height = 6)
    
    if(x == "score" & y == "umap") {
      p_UMAP_base <- p
    }
    if(x == "score" & y == "cons_v_asa") {
      p_cons_v_asa_base <- p
    }
    
  }
  
  files <- list.files(tmp_dir)
  
  files <- sub(".svg$", "", files)
  
  html_slide_show(svg_directory = tmp_dir,
                  output_file = paste0("~/peptide_alg/", y, ".html"),
                  frames = files,
                  categories = NULL,
                  title = y,
                  columns = 1)
  
}



