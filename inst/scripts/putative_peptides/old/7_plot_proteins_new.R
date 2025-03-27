##plot proteins

adhesionGPCRs <- list(`Group I` = paste0("ADGRG", 1:6),
                      `Group II` = paste0("ADGRL", 1:3),
                      `Group III` = paste0("ADGRB", 1:3),
                      `Group IV` = paste0("ADGRE", 1:5),
                      `Group V` = paste0("ADGRF", 1:5),
                      `Group VI` = paste0("ADGRA", 1:3),
                      `Group VII` = paste0("ADGRV1"),
                      PARS = paste0("PAR", 1:3))

chemokines <- unique(secretome_roi[["gene"]][str_detect(secretome_roi[["gene"]], "^CCL|^CXCL|^XCL|^CX3CL")])

chemokines <- lapply(c("^CCL", "^CXCL", "^XCL", "^CX3CL"), \(x) chemokines[str_detect(chemokines, x)])

names(chemokines) <- c("CCL", "CXCL", "XCL", "CX3CL")

known_ligands <- secretome_roi %>% 
  filter(roi_type == "gtp") %>%
  distinct(gene) %>%
  pull(gene)





#to_visualize <- list(examples = c("RLN2", "TGFA", "CORT", "GIP", "NPS", "GHRH", "KISS1", "PPY", "EREG", "F2", "NPY", "NPPC"))

#to_visualize <- list(examples = c("RARRES2", "CXCL17", "GPR15LG", "KNG1"))



to_visualize <- c(list(goi = c("RARRES2", "NICOL1", "DMKN", "CARTPT", "GPR15LG"),
                       `known ligands` = known_ligands), 
                  chemokines)

to_visualize <- lapply(to_visualize, \(x) x[gtools::mixedorder(x)])



to_visualize <- list(`top 100 hits` = secretome_roi %>%
                         filter(roi_length > 3) %>%
                         mutate(secreted_strict = map_lgl(annotations, .f = \(x) {
                           to_test <- x %>%
                             filter(annotation_type == "subcellular location") %>%
                             pull("annotation")
                           any(grepl("Secreted", to_test)) & !any(grepl("membrane|cytoplasm", to_test, ignore.case = TRUE))
                         })) %>%
                        filter(secreted_strict) %>%
                         group_by(gene) %>%
                         mutate(max_score_per_gene = max(score_nn_new_s8_entire)) %>%
                         filter(score_nn_new_s8_entire == max_score_per_gene) %>%
                         ungroup() %>%
                         arrange(desc(max_score_per_gene)) %>%
                         dplyr::slice(1:100) %>%
                         pull(gene))



to_visualize <-  list(`N-terminal` = secretome_roi %>%
                                        filter(roi_length > 3) %>%
                                        filter(!is.na(`percent_ol_phs_hsr:phs_N`)) %>%
                                        arrange(desc(score_nn10c__entire)) %>%
                                        dplyr::slice(1:30) %>%
                                        pull(gene),
                      `C-terminal` = secretome_roi %>%
                                        filter(roi_length > 3) %>%
                                        filter(!is.na(`percent_ol_phs_hsr:phs_C`)) %>%
                                        arrange(desc(score_nn10c__entire)) %>%
                                        dplyr::slice(1:30) %>%
                                        pull(gene)
                      )
  
                      





all_genes <- do.call(c, to_visualize)
nice_mode <- FALSE


names(all_genes) <- do.call(c, lapply(names(to_visualize), \(x) rep(x, length(to_visualize[[x]]))))

plot_dat <- secretome_roi %>%
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


plot_features <- c("strand", "helix", "turn", "signal peptide", "Dibasic", "Cysteine", "W or Y", "c-term", "n-term", "disulfide bond", names(ss_features))

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
                 setNames("#197EC0FF", "W or Y"),
                 setNames("#197EC000", "disulfide bond"))

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
  
  
  p_cons <- ggplot2::ggplot(to_plot) +
    ggplot2::geom_col(aes(x = index, y = frequency), width = 0.8) +
    ggplot2::geom_line(aes(x = index, y = doubleSmooth), 
                       color = "#8A919799", 
                       linewidth = 2) + 
    geom_hline(yintercept = cons_cutoff, linetype = "dashed", color = "red") +
    geom_text(aes(x = index, y = -0.3, label = AA, color = is_conserved,
                  fontface = ifelse(is_conserved == "conserved", 2, 1)), 
              size = 3) +
    ggplot2::scale_color_manual(name = "conservation", 
                                values = c(setNames("red", "conserved"), 
                                           setNames("black", "not conserved"))) +
    ggplot2::expand_limits(x = c(0, p_size)) +
    scale_x_continuous(expand = expansion(mult = buffer),
                       breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
                       minor_breaks = seq(0, p_size, by = 10)) +
    ggplot2::xlab("") +
    ggplot2::ylab("conservation (Aminode)") +
    ggnewscale::new_scale_color() +
    ggplot2::theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),  
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(linetype = "solid"),
      panel.grid.minor.x = element_line(linetype = "dashed")
    )
  
  
  
  
  
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
          mean(scales::rescale(pd$aa_scores[[1]][[roi_color]], to = c(0,1))[start:end], na.rm = TRUE)
        }))
      col_var <- "score"
      
    }
    else {
      feat_dat <- feats[[feat]] %>%
        mutate(feat_len = end - start) %>%
        arrange(desc(feat_len))
    }
    
    if(nice_mode) {
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
                         dpi = 600)
    } else {
      p_feats <- p_feats +
        ggplot2::geom_segment(data = feat_dat, 
                              aes(x = start, 
                                  xend = end, 
                                  y = source, 
                                  color = !!sym(col_var)), 
                                  linewidth = 4) + 
        ggplot2::geom_point(data = feat_dat %>% filter(start == end & source != "ensemble"), 
                            aes(x = start, y = source, color = !!sym(col_var)), 
                            pch = 19, 
                            size = 3)
      
    }
    p_feats <- p_feats +
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
    
    curve_dat <- feat_dat %>% filter(type == "disulfide bond") %>% drop_na(start, end)
    
    if(nrow(curve_dat) > 0) {
    p_feats <- p_feats +
      ggplot2::geom_curve(data = curve_dat,
                          aes(x = start, xend = end, y = "sites", yend = "sites"),
                          color = "grey40", curvature = 0.4) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::expand_limits(x = c(0, p_size + 0.4)) +
      scale_x_continuous(expand = expansion(mult = buffer),
                         breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
                         minor_breaks = seq(0, p_size, by = 10))
    }
    
    
    if(grepl("^phs|hsr", feat)) {
      p_feats <- p_feats +
        ggplot2::scale_color_viridis_c(name = paste("ROI score", roi_color, sep = "\n"), option = "H", limits = c(0,1)) +
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
  
  if("mean_af_missense" %in% colnames(pd$af_missense_mapped[[1]][["ms"]])) {
  p_afm_dat <- pd$af_missense_mapped[[1]][["ms"]] %>%
    mutate(index = 1:nrow(.), .before = everything()) %>%
    dplyr::rename(mean = mean_af_missense) %>%
    dplyr::select(-AA) %>%
    reshape2::melt(id.vars = "index", variable.name = "afm")
  } else {
    p_afm_dat <- tibble(index = 1:p_size,
                        afm = as.numeric(rep(NA, p_size)),
                        value = as.numeric(rep(NA, p_size)))
  }
  
  p_afm <- ggplot2::ggplot(p_afm_dat) +
    ggplot2::geom_tile(aes(x = index, y = afm, fill = value), show.legend = FALSE) + 
    ggplot2::scale_fill_viridis_c(option = "H") +
    scale_x_continuous(expand = expansion(mult = buffer)) +
    ylab("per AA AF missense") +
    ggplot2::theme_bw()
  
  
  p_feature_dat <-  pd$aa_scores[[1]] %>%
    mutate(index = 1:nrow(.), .before = everything()) %>%
    dplyr::select(starts_with("index"),
                  starts_with("conservation_og"), 
                  starts_with("pathogenicity"),
                  starts_with("relASA"),
                  all_of(params_tp)) %>%
  dplyr::select(-matches("_lead\\d|lag\\d")) %>%
    reshape2::melt(id.vars = "index", variable.name = "params") %>%
    mutate(param_type = str_extract(params, "^[^_]+"))
    
  
  
  
  p_scores_dat <- pd$aa_scores[[1]] %>%
    mutate(index = 1:nrow(.), .before = everything()) %>%
    dplyr::select(starts_with("index"), 
                  starts_with("score_")) %>%
    reshape2::melt(id.vars = "index", variable.name = "params") %>%
    filter(params %in% nn_tp[["params"]]) %>%
    left_join(., nn_tp, by = "params")
  
  derivative <- function(x) {
    subsetter <- !is.na(x) & !is.nan(x) & !is.infinite(x)
    to_return <- as.numeric(rep(NA, length(x)))
    to_return[subsetter] <- tryCatch({predict(smooth.spline(x[subsetter]), deriv = 1)[["y"]]}, 
                                     error = function(e) {as.numeric(rep(NA, sum(subsetter)))})
    return(to_return)
  }
  
  p_scores_roc <- p_scores_dat %>%
    group_by(params) %>%
    mutate(value = derivative(value))
  
  p_feat_dat <- ggplot2::ggplot(p_feature_dat) +
    ggplot2::geom_tile(aes(x = index, y = params, fill = value), show.legend = x == "parameter") +
    ggh4x::facet_nested(rows = vars(param_type), scales = "free_y", switch = "y") +
    # geom_text(data = to_plot %>% mutate(nn = last_facet), aes(x = index, y = 0, label = AA, color = is_conserved,
    #                               fontface = ifelse(is_conserved == "conserved", 2, 1)), 
    #           size = 3, show.legend = FALSE) +
    ggplot2::scale_color_manual(name = "conservation", 
                                values = c(setNames("red", "conserved"), 
                                           setNames("black", "not conserved"))) +
    scale_x_continuous(expand = expansion(mult = buffer)) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) + 
    xlab("") +
    ylab(x) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red")
  
  
  
  p_scores <- ggplot2::ggplot(p_scores_dat) +
    ggplot2::geom_tile(aes(x = index, y = params, fill = value), show.legend = x == "parameter") +
      ggh4x::facet_nested(rows = vars(model_type, nn_type2, nn_type), scales = "free_y", switch = "y") +
      # geom_text(data = to_plot %>% mutate(nn = last_facet), aes(x = index, y = 0, label = AA, color = is_conserved,
      #                               fontface = ifelse(is_conserved == "conserved", 2, 1)), 
      #           size = 3, show.legend = FALSE) +
      ggplot2::scale_color_manual(name = "conservation", 
                                  values = c(setNames("red", "conserved"), 
                                             setNames("black", "not conserved"))) +
    scale_x_continuous(expand = expansion(mult = buffer)) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) + 
    xlab("") +
    ylab(x) +
    ggplot2::theme_bw() +
    theme(axis.text.y = element_blank()) +
    ggplot2::scale_fill_viridis_c(option = "H")
      
  
  
  
  design <- "
  1
  1
  1
  2
  2
  2
  2
  2
  2
  2
  3
  3
 3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  3
  4
  4
  4
  4
  6
  6
  6
"
  
  anno_dat <- secretome_roi %>%
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
  
  anno1 <- pd$annotations[[1]] %>%
    filter(annotation_name == "comment") %>%
    filter(annotation_type %in% c("subcellular location", "tissue specificity", "disease")) %>%
    {paste0(.[["annotation_type"]], ": ", .[["annotation"]])}
  
  anno2 <- pd$annotations[[1]] %>%
    filter(annotation_name == "dbReference" & name_1 == "disease") %>% 
    {paste0(.[["annotation_type"]], ": ", .[["annotation"]])}
  
  
  really_long_subtitle <- c(paste0("full name: ", pd$full_name[1]),
                            paste0("known ligand: ", feats[["gtp"]][["type"]], collapse = "; "),
                            anno1,
                            anno2)
  
  really_long_subtitle <- paste(really_long_subtitle, collapse = "\n")
  
  p_al <- 
    p_cons + 
    p_feat_dat +
    p_scores +
    p_feats +
    p_afm +
    patchwork::plot_layout(design = design, ncol = 2, axes = "collect_x", widths = c(pw, ph)) +
    patchwork::plot_annotation(title = x,
                               subtitle = really_long_subtitle,
                               theme = theme(plot.title = element_text(face = 2,
                                                                       size = 22,
                                                                       hjust = 0),
                                             plot.subtitle = element_text(face = 1,
                                                                          size = 12,
                                                                          hjust = 0)))
  
  
  ggsave(filename = paste0(plot_dir, "/", x, ".svg"),
         plot = p_al, svglite::svglite, width = plot_width, height = 30, limitsize = FALSE)
  
  
})




files <- list.files(plot_dir)

files <- sub(".svg$", "", files)

all_genes2 <- all_genes[all_genes %in% files]


html_slide_show(svg_directory = plot_dir,
                output_file = "~/peptide_alg/chemokines_knowns.html",
                frames = all_genes2,
                categories = names(all_genes2),
                title = "chemokines_knowns",
                columns = 3)


