
to_visualize <- list(goi = c("RARRES2", "KNG1", "GCG", "CARTPT"))
to_visualize <- lapply(to_visualize, \(x) x[gtools::mixedorder(x)])


all_genes <- do.call(c, to_visualize)


names(all_genes) <- do.call(c, lapply(names(to_visualize), \(x) rep(x, length(to_visualize[[x]]))))

plot_dat <- secretome %>%
  filter(gene %in% all_genes) %>%
  group_by(gene) %>%
  filter(nchar(sequence_uni) == max(nchar(sequence_uni))) %>%
  #slice_head(n = 1) %>%
  ungroup() %>%
  split(., .[["gene"]])


ss_features <- c(setNames("Alpha helix (4-12)", "H"),
                 setNames("Isolated beta-bridge residue", "B"),
                 setNames("Strand", "E"),
                 setNames("3-10 helix", "G"),
                 setNames("Pi helix", "I"),
                 setNames("Turn", "T"),
                 setNames("Bend", "S"),
                 setNames("Kappa helix", "P"))


plot_features <- c("strand", "helix", "turn", "signal peptide", "Dibasic", "Cysteine", "W or Y", "c-term", "n-term", names(ss_features))

extra_color <- c("#FD7446FF","#FD8CC1FF","#197EC0FF","#1A9993FF",
                 "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF",
                 "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF",
                 "#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF")

color_scale <- c(setNames("#FED439FF", "H"),
                 setNames("#FED439FF", "helix"),
                 setNames("#709AE1FF", "B"),
                 setNames("#FD7446FF", "E"),
                 setNames("#FD7446FF", "strand"),
                 setNames("#D5E4A2FF", "G"),
                 setNames("#46732EFF", "I"),
                 setNames("#71D0F5FF", "T"),
                 setNames("#71D0F5FF", "turn"),
                 setNames("#C80813FF", "S"),
                 setNames("#D2AF81FF", "P"),
                 setNames("#370335FF", "signal peptide"),
                 setNames("#FD8CC1FF", "Dibasic"),
                 setNames("#91331FFF", "c-term"),
                 setNames("#1A9993FF", "n-term"),
                 setNames("#075149FF", "Cysteine"),
                 setNames("#197EC0FF", "W or Y"))

labels <- sapply(names(color_scale), \(x) if(x %in% names(ss_features)) return(ss_features[[x]]) else return(x))


plot_dir <- "~/peptide_alg/plots_tmp11"

dir.create(plot_dir)

future::plan(strategy = future::sequential())

furrr::future_map(names(plot_dat), \(x) {
  
  message("plotting ", x)
  
  pd <- plot_dat[[x]]
  
  
  feats <- pd[["features"]][[1]] %>%
    filter(type %in% plot_features | grepl("^phs|^hsr", source) | source %in% c("gtp", "ensemble"))
  
  feats[["source"]] <- factor(feats[["source"]], levels = rev(c("AA", 
                                                                "ensemble", 
                                                                "alpha fold", 
                                                                "sites", 
                                                                "phs_n-term-Cys", 
                                                                "phs_c-term-Cys",
                                                                "phs_n-term",
                                                                "phs_c-term",
                                                                "phs_dbr",
                                                                "phs", 
                                                                "phs_hsr",
                                                                "phs_hsr_N",
                                                                "phs_hsr_C",
                                                                "hsr",
                                                                "gtp", 
                                                                "uniprot")), ordered = TRUE)
  
  p_order <- levels(feats[["source"]])
  
  cons_dat <- pd[["cons_mapped"]][[1]][["ms"]]
  
  
  if("AA" %in% names(cons_dat)) {
    to_plot <- bind_cols(pd[["af_mapped"]][[1]][["ms"]], cons_dat %>% dplyr::select(-index, -AA))
  } else {
    to_plot <- pd[["af_mapped"]][[1]][["ms"]] %>%
      mutate(cons = NA,
             frequency = NA,
             window = NA,
             smooth = NA,
             doubleSmooth = NA)
  }
  
  p_size <- nrow(to_plot)
  
  to_plot[["index"]] <- 1:p_size
  
  buffer <- 1/p_size
  
  feats <- bind_rows(feats, tibble(type = NA, evidence = NA, start = NA, end = NA, source = "AA"))
  
  cons_cutoff <- 0.2
  is_conserved <- ifelse(to_plot[["frequency"]] < cons_cutoff, "conserved", "not conserved")
  is_conserved[is.na(is_conserved)] <- "not conserved"
  open_cutoff <- 0.5
  is_open <- ifelse(to_plot[["relASA"]] > open_cutoff, 2, 1)
  is_open[is.na(is_open)] <- 1
  
  not_beta_cutoff <- 0.2
  is_not_beta <- ifelse(to_plot[["frequency"]] < cons_cutoff, 2, 1)
  font_color_scale <- c("black", "red")
  

  
  p_feats <- ggplot2::ggplot(to_plot) 
  
  feats <- feats %>%
    drop_na(source) %>%
    split(., .[["source"]])
  
  for(feat in p_order[p_order %in% names(feats)]) {
    show_legend <- TRUE
    col_var <- "type"
    
    if(feat == "AA") {show_legend <- FALSE}
    if(feat == "gtp") {
      extra_feats <- feats[[feat]][["type"]]
      color_scale2 <- extra_color[1:length(extra_feats)]
      names(color_scale2) <- extra_feats
      color_scale <- c(color_scale, color_scale2)
    }
    if(feat == "ensemble") {
      feat_dat <- feats[[feat]] %>%
        filter(start != max(start)) %>%
        mutate(feat_len = end - start) %>%
        arrange(desc(feat_len))
      
    } else if(grepl("^phs|hsr", feat)) {
      feat_dat <- feats[[feat]] %>%
        mutate(feat_len = end - start) %>%
        arrange(desc(feat_len)) %>%
        mutate(score = map2_dbl(.x = start, .y = end, .f = \(start, end) {
          mean(scales::rescale(pd$aa_scores[[1]][["score_nn4_s8"]], to = c(0,1))[start:end], na.rm = TRUE)
        }))
      col_var <- "score"
      
    }
    else {
      feat_dat <- feats[[feat]] %>%
        mutate(feat_len = end - start) %>%
        arrange(desc(feat_len))
    }
    p_feats <- p_feats +
      ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_segment(data = feat_dat, 
                                                                     aes(x = start, 
                                                                         xend = end, 
                                                                         y = source, 
                                                                         color = !!sym(col_var)), 
                                                                     linewidth = 4),
                                               colour = "black",
                                               sigma = 1, 
                                               expand = 5),
                         dev = "ragg",
                         dpi = 600) + 
      ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_point(data = feat_dat %>% filter(start == end & source != "ensemble"), 
                                                                   aes(x = start, y = source, color = !!sym(col_var)), 
                                                                   pch = 19, 
                                                                   size = 3),
                                               colour = "black",
                                               sigma = 1, 
                                               expand = 5),
                         dev = "ragg",
                         dpi = 600) + 
      ggplot2::geom_point(data = feat_dat %>% filter(start == end & source == "ensemble"), 
                          aes(x = start, y = source), 
                          pch = 25, stroke = 1, color = "grey40",
                          fill = "#FF573300",
                          size = 2) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::expand_limits(x = c(0, p_size + 0.4)) +
      scale_x_continuous(expand = expansion(mult = buffer),
                         breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
                         minor_breaks = seq(0, p_size, by = 10))
    
    if(grepl("^phs|hsr", feat)) {
      p_feats <- p_feats +
        ggplot2::scale_color_viridis_c(name = "ROI score", option = "H", limits = c(0,1)) +
        ggnewscale::new_scale_color()
    } else if(!feat %in% "ensemble"){
      p_feats <- p_feats +
        ggplot2::scale_color_manual(name = feat, values = color_scale, labels = labels) +
        guides(color = guide_legend(nrow = 6)) +
        ggnewscale::new_scale_color()
    }
  }
  
  
  
  p_feats <- p_feats +
    ggplot2::theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),  # Remove major horizontal grid lines
      panel.grid.minor.y = element_blank(), 
      panel.grid.major.x = element_line(linetype = "solid"),
      panel.grid.minor.x = element_line(linetype = "dashed"),
      legend.position = "bottom",
      legend.justification = "left",
      legend.title.position = "top"
    )
  
  p_feats <- p_feats + ggplot2::geom_text(aes(x = 1:p_size, 
                                              y = 1, 
                                              label = to_plot[["AA"]], 
                                              color = is_conserved,
                                              fontface = ifelse(is_conserved == "conserved", 2, 1)), 
                                          size = 3) +
    ggplot2::scale_color_manual(name = "conservation", 
                                values = c(setNames("red", "conserved"), 
                                           setNames("black", "not conserved")))
  
  
  
  score_cats <- tibble(
    params = c("score_og", "score_og_a4", "score_og_a6", 
               "score_og_a8", "score_glm", "score_nn", "score_nn2", "score_nn3", 
               "score_nn4", "score_nn5", "score_nn6", "score_b4", "score_b6", "score_b8", "score_glm_s", 
               "score_nn_s", "conservation", "relASA", "pathogenicity", "conservation_s4", 
               "conservation_s6", "conservation_s8", "relASA_s4", "relASA_s6", 
               "relASA_s8", "pathogenicity_s4", "pathogenicity_s6", "pathogenicity_s8", "score_nn_s4", "score_nn_s6", "score_nn_s8", "score_nn2_s4", 
               "score_nn2_s6", "score_nn2_s8", "score_nn3_s4", "score_nn3_s6", 
               "score_nn3_s8", "score_nn4_s4", "score_nn4_s6", "score_nn4_s8", 
               "score_nn5_s4", "score_nn5_s6", "score_nn5_s8", "score_nn6_s4", 
               "score_nn6_s6", "score_nn6_s8", "score_nn_s_s4", "score_nn_s_s6", 
               "score_nn_s_s8"),
    raw_smooth = c("raw", 
                   "smooth",
                   "smooth",
                   "smooth",
                   "raw",
                   "raw",
                   "raw",
                   "raw",
                   "raw",
                   "raw",
                   "raw",
                   "smooth",
                   "smooth",
                   "smooth",
                   "smooth",
                   "smooth",
                   "raw",
                   "raw", 
                   "raw",
                   rep("smooth", 30)
    ),
    type = c(rep("score", 16), rep("parameter", 12), rep("score", 21))
  )
  
  score_cats <- score_cats %>%
                  filter(grepl("score_nn4_", params))
  
  p_scores_dat <- pd$aa_scores[[1]] %>%
    mutate(index = 1:nrow(.), .before = everything()) %>%
    mutate(across(contains("_nn") | contains("_glm"), .fns = \(x) x - min(x, na.rm = TRUE))) %>%
    mutate(across(contains("_nn") | contains("_glm"), .fns = \(x) x/max(x, na.rm = TRUE))) %>%
    dplyr::select(c("index", score_cats$params)) %>%
    reshape2::melt(id.vars = "index", variable.name = "params")
  p_scores_dat <- left_join(p_scores_dat, score_cats, by = "params") %>%
    split(., .[["type"]])
  
  derivative <- function(x) {
    subsetter <- !is.na(x) & !is.nan(x) & !is.infinite(x)
    to_return <- as.numeric(rep(NA, length(x)))
    to_return[subsetter] <- tryCatch({predict(smooth.spline(x[subsetter]), deriv = 1)[["y"]]}, 
                                     error = function(e) {as.numeric(rep(NA, sum(subsetter)))})
    return(to_return)
  }
  
  p_scores_dat[["rate_of_change"]] <- p_scores_dat$score %>%
    group_by(params) %>%
    mutate(value = derivative(value))
  
  
  p_scores <- lapply(names(p_scores_dat), \(x) {
    p <- ggplot2::ggplot(p_scores_dat[[x]]) +
      ggplot2::geom_tile(aes(x = index, y = params, fill = value), show.legend = x == "parameter") +
      geom_text(data = to_plot, aes(x = index, y = -0.2, label = AA, color = is_conserved,
                                    fontface = ifelse(is_conserved == "conserved", 2, 1)), 
                size = 3, show.legend = FALSE) +
      ggplot2::scale_color_manual(name = "conservation", 
                                  values = c(setNames("red", "conserved"), 
                                             setNames("black", "not conserved"))) +
      scale_x_continuous(expand = expansion(mult = buffer)) +
      scale_y_discrete(expand = expansion(mult = c(1, 0))) + 
      xlab("") +
      ylab(x) +
      ggplot2::theme_bw()
    
    if(x == "rate_of_change") {
      p <- p + scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", 
        midpoint = 0
      )
    } else {
      p <- p + ggplot2::scale_fill_viridis_c(option = "H")
      
    }
    p
    
  }
  )
  
  
  
  design <- "
1,
2,
3
"
  
  anno_dat <- secretome %>%
    filter(gene == x & roi_name %in% feats[["gtp"]][["type"]])
  
  # p_UMAP <-  p_UMAP_base + 
  #   ggnewscale::new_scale_color() +
  #   
  #   ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_point(data = anno_dat,
  #                                                                mapping = aes(x = UMAP1, 
  #                                                                              y = UMAP2, 
  #                                                                              color = dbr_name), 
  #                                                                size = 5, 
  #                                                                pch = 21, 
  #                                                                fill = "#FF573300",
  #                                                                stroke = 2, 
  #                                                                show.legend = FALSE),
  #                                            colour = "black",
  #                                            sigma = 3, 
  #                                            expand = 8),
  #                      dev = "ragg",
  #                      dpi = 600) + 
  #   ggplot2::scale_color_manual(values = color_scale) +
  #   theme(legend.position = "none")
  # 
  # p_cons_v_asa <- p_cons_v_asa_base +
  #   ggnewscale::new_scale_color() +
  #   ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_point(data = anno_dat,
  #                                                                mapping = aes(x = prox_cons_3, 
  #                                                                              y = prox_relASA_3, 
  #                                                                              color = dbr_name), 
  #                                                                size = 5, 
  #                                                                pch = 21, 
  #                                                                fill = "#FF573300",
  #                                                                stroke = 2, 
  #                                                                show.legend = FALSE),
  #                                            colour = "black",
  #                                            sigma = 3, 
  #                                            expand = 8),
  #                      dev = "ragg",
  #                      dpi = 600) +
  #   ggplot2::scale_color_manual(values = color_scale)
  # 
  
  
  plot_width <- max(p_size * (13/150), 10) + 6
  
  # x + y = plot_width
  # 
  # x/y * plot_width = 6
  # 
  # x + x * plot_width/6 = plot_width
  # 
  pw = plot_width/(1+ (plot_width/6))
  
  ph = plot_width - pw
  
  p_al <- 
    p_scores[[1]] +
    p_scores[[2]] +
    p_feats +
    patchwork::plot_layout(ncol = 1, axes = "collect_x", widths = c(pw, ph)) +
    patchwork::plot_annotation(title = x,
                               subtitle = paste(feats[["gtp"]][["type"]], collapse = "; "),
                               theme = theme(plot.title = element_text(face = 2,
                                                                       size = 22,
                                                                       hjust = 0),
                                             plot.subtitle = element_text(face = 2,
                                                                          size = 18,
                                                                          hjust = 0)))
  
  
  ggsave(filename = paste0(plot_dir, "/", x, ".svg"),
         plot = p_al, svglite::svglite, width = plot_width, height = 6, limitsize = FALSE)
  
  
})




files <- list.files(plot_dir)

files <- sub(".svg$", "", files)

all_genes2 <- all_genes[all_genes %in% files]


html_slide_show(svg_directory = plot_dir,
                output_file = "~/peptide_alg/known_peptides_latest3.html",
                frames = all_genes2,
                categories = names(all_genes2),
                title = "known_peptides_latest3",
                columns = 1)
