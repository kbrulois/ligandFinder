##plot proteins


to_visualize <- list(examples = c("EDN1", "EDN2", "EDN3", "KNG1", "CXCL12", "GCG", "AGT"),
                     top10_relative_prox = secretome %>%
                       arrange(desc(relative_prox_cons_4)) %>%
                       slice_head(n = 1) %>%
                       distinct(gene) %>%
                       pull(gene),
                     score = secretome %>%
                       arrange(desc(score)) %>%
                       slice_head(n = 1) %>%
                       distinct(gene) %>%
                       pull(gene)
)


to_visualize <- chemokines






to_visualize <- list(examples = c("CARTPT", "C17ORF67", "C15ORF61", "RARRES2", "NICOL1", "DMKN"),
                     `known ligands` = secretome %>%
                                          filter(dbr_type == "gtp") %>%
                                          distinct(gene) %>%
                                          pull(gene))


to_visualize <- lapply(to_visualize, \(x) x[gtools::mixedorder(x)])


all_genes <- do.call(c, to_visualize)


names(all_genes) <- do.call(c, lapply(names(to_visualize), \(x) rep(x, length(to_visualize[[x]]))))

plot_dat <- secretome %>%
  filter(gene %in% all_genes) %>%
  group_by(gene) %>%
  filter(roi_length == max(roi_length)) %>%
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


plot_dir <- "~/peptide_alg/plots_tmp10"

dir.create(plot_dir)

future::plan(strategy = future::sequential())

furrr::future_map(names(plot_dat)[36:93], \(x) {

  message("plotting ", x)

  pd <- plot_dat[[x]]

  score_dat <- secretome %>%
                filter(gene == x) %>%
                dplyr::select(dbr_name, score) %>%
                dplyr::rename(type = dbr_name)

  feats <- pd[["features"]][[1]] %>%
    filter(type %in% plot_features | grepl("^dbr_|^c-term_|^n-term_", type) | source %in% c("gtp", "ensemble"))

  feats[["source"]] <- factor(feats[["source"]], levels = rev(c("AA", "ensemble", "alpha fold", "sites", "dbr", "c-term", "n-term", "gtp", "uniprot")))

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



    # p_Phi <- ggplot2::ggplot(to_plot) +
    #   ggplot2::geom_col(data = to_plot,
    #                     aes(x = index, y = Phi),
    #                     position = position_dodge(),
    #                     width = 0.8) +
    #   ggplot2::expand_limits(x = c(0, p_size)) +
    #   scale_x_continuous(expand = expansion(mult = buffer),
    #                      breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
    #                      minor_breaks = seq(0, p_size, by = 10)) +
    #   ggplot2::xlab("") +
    #   ggplot2::ylab("Phi") +
    #   ggnewscale::new_scale_color() +
    #   ggplot2::theme_bw() +
    #   theme(
    #     panel.grid.major.y = element_blank(),
    #     panel.grid.minor.y = element_blank(),
    #     panel.grid.major.x = element_line(linetype = "solid"),
    #     panel.grid.minor.x = element_line(linetype = "dashed")
    #   )
    #
    # p_Psi <- ggplot2::ggplot(to_plot) +
    #   ggplot2::geom_col(aes(x = index, y = Psi),
    #                     position = position_dodge(),
    #                     width = 0.8) +
    #   ggplot2::expand_limits(x = c(0, p_size)) +
    #   scale_x_continuous(expand = expansion(mult = buffer),
    #                      breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
    #                      minor_breaks = seq(0, p_size, by = 10)) +
    #   ggplot2::xlab("") +
    #   ggplot2::ylab("Psi") +
    #   ggnewscale::new_scale_color() +
    #   ggplot2::theme_bw() +
    #   theme(
    #     panel.grid.major.y = element_blank(),
    #     panel.grid.minor.y = element_blank(),
    #     panel.grid.major.x = element_line(linetype = "solid"),
    #     panel.grid.minor.x = element_line(linetype = "dashed")
    #   )

    p_ASA <- ggplot2::ggplot(to_plot) +
      ggplot2::geom_col(aes(x = index, y = relASA), width = 0.8) +

      ggplot2::geom_line(aes(x = index, y = relASA_sss),
                         color = "#8A919799",
                         linewidth = 2) +
      geom_hline(yintercept = open_cutoff, linetype = "dashed", color = "red") +
      ggplot2::expand_limits(x = c(0, p_size)) +

      scale_x_continuous(expand = expansion(mult = buffer),
                         breaks = seq(0, p_size, by = 20),    # Major grid lines every 20
                         minor_breaks = seq(0, p_size, by = 10)) +
      ggplot2::xlab("") +
      ggplot2::ylab("relativeASA") +
      ggplot2::theme_bw() +
      ggnewscale::new_scale_color() +
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

      } else if(feat %in% c("c-term", "dbr", "n-term")) {
        feat_dat <- feats[[feat]] %>%
          mutate(feat_len = end - start) %>%
          arrange(desc(feat_len))
        feat_dat <- left_join(feat_dat, score_dat, by = "type")
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

    if(feat %in% c("c-term", "dbr", "n-term")) {
      p_feats <- p_feats +
        ggplot2::scale_color_viridis_c(name = "ROI score", option = "B", limits = c(0,1)) +
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


     p_afm_dat <- pd$af_missense_mapped[[1]][["ms"]] %>%
                    mutate(index = 1:nrow(.), .before = everything()) %>%
                    dplyr::rename(mean = mean_af_missense) %>%
                    dplyr::select(-AA) %>%
                    reshape2::melt(id.vars = "index", variable.name = "afm")

     p_afm <- ggplot2::ggplot(p_afm_dat) +
              ggplot2::geom_tile(aes(x = index, y = afm, fill = value)) +
              ggplot2::scale_fill_viridis_c(option = "B") +
       scale_x_continuous(expand = expansion(mult = buffer)) +
       ggplot2::theme_bw()







    design <- "
  13
  13
  13
  14
  24
  24
  25
  25
  #5
  #6
  #6
  #6
"

      anno_dat <- secretome %>%
                    filter(gene == x & dbr_name %in% feats[["gtp"]][["type"]])

      p_UMAP <-  p_UMAP_base +
                  ggnewscale::new_scale_color() +

        ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_point(data = anno_dat,
                                      mapping = aes(x = UMAP1,
                                                    y = UMAP2,
                                                    color = dbr_name),
                                      size = 5,
                                      pch = 21,
                                      fill = "#FF573300",
                                      stroke = 2,
                                      show.legend = FALSE),
                                      colour = "black",
                                      sigma = 3,
                                      expand = 8),
                           dev = "ragg",
                           dpi = 600) +
                ggplot2::scale_color_manual(values = color_scale) +
                theme(legend.position = "none")

    p_cons_v_asa <- p_cons_v_asa_base +
      ggnewscale::new_scale_color() +
      ggrastr::rasterise(ggfx::with_outer_glow(ggplot2::geom_point(data = anno_dat,
                          mapping = aes(x = prox_cons_3,
                                        y = prox_relASA_3,
                                        color = dbr_name),
                          size = 5,
                          pch = 21,
                          fill = "#FF573300",
                          stroke = 2,
                          show.legend = FALSE),
    colour = "black",
    sigma = 3,
    expand = 8),
dev = "ragg",
dpi = 600) +
      ggplot2::scale_color_manual(values = color_scale)



    plot_width <- max(p_size * (13/150), 10) + 6

    # x + y = plot_width
    #
    # x/y * plot_width = 6
    #
    # x + x * plot_width/6 = plot_width
    #
   pw = plot_width/(1+ (plot_width/6))

    ph = plot_width - pw

    p_al <- p_UMAP +
      p_cons_v_asa +
      p_ASA +
      #p_Phi +
      #p_Psi +
      p_cons +
      p_afm +
      p_feats +
      patchwork::plot_layout(design = design, ncol = 2, axes = "collect_x", widths = c(pw, ph)) +
      patchwork::plot_annotation(title = x,
                                 subtitle = paste(feats[["gtp"]][["type"]], collapse = "; "),
                                 theme = theme(plot.title = element_text(face = 2,
                                                                         size = 22,
                                                                         hjust = 0),
                                               plot.subtitle = element_text(face = 2,
                                                                            size = 18,
                                                                            hjust = 0)))


    ggsave(filename = paste0(plot_dir, "/", x, ".svg"),
           plot = p_al, svglite::svglite, width = plot_width, height = 13, limitsize = FALSE)


})




files <- list.files(plot_dir)

files <- sub(".svg$", "", files)

all_genes2 <- all_genes[all_genes %in% files]


html_slide_show(svg_directory = plot_dir,
                output_file = "~/peptide_alg/known_peptides_latest.html",
                frames = all_genes2,
                categories = names(all_genes2),
                title = "known_peptides_new",
                columns = 3)


