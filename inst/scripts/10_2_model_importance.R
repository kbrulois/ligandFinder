




saliency_fn <- function(model, x_batch) {
  x_tensor <- tf$convert_to_tensor(x_batch)

  with(tf$GradientTape() %as% tape, {
    tape$watch(x_tensor)
    preds <- model(x_tensor)[[1]]  # global output
  })

  grads <- tape$gradient(preds, x_tensor)
  as.array(grads)
}


targs <- list(
  c("C", "loop_C"),
  c("N", "loop_N")
)


met_it <- list(core = all_params3[1:4],
               motif1 = c("AA_Y", "AA_W", "AA_F"),
               motif2 = c("AA_Q", "AA_E", "AA_D"),
               motif3 = c("AA_G", "AA_H", "AA_V"),
               SS = c("SS_G", "SS_B", "SS_H"))

for(y in seq_along(met_it)) {

  mets <- met_it[[y]]
  p_final <- list()


  for(x in seq_along(targs)) {

    targ <- targs[[x]]

    shap_vals <- saliency_fn(models[[targ[1]]], nn_in_C$val$x)


    dat_toplot <- lapply(mets, \(x) shap_vals[nn_in_C$val$y_global == 1,,which(all_params3 == x)] %>% as_tibble %>% mutate(metric = x)) %>%
      bind_rows() %>%
      pivot_longer(cols = -metric) %>%
      mutate(index = setNames(1:36, paste0("V", 1:36))[name]) %>%
      select(-name)

    non_zero_mean <- function(x) {
      mean(x[!x == 0], na.rm = TRUE)
    }

    center_func <- if(names(met_it)[y] == "core") {center_func <- non_zero_mean} else {center_func <- mean}

    dat_toplot2 <- map(known_dat$data[known_dat$target %in% targ], \(x) {

      x[, mets] %>%
        mutate(index = row_number())

    }) %>%
      bind_rows %>%
      group_by(index) %>%
      summarise(across(everything(), center_func),
                index = first(index)) %>%
      pivot_longer(cols = -index) %>%
      mutate(type = "average of residues within or adjacent to known GPCR ligands")

    dat_toplot3 <- map(c_dat$data[c_dat$target %in% targ], \(x) {

      x[, mets] %>%
        mutate(index = row_number())

    }) %>%
      bind_rows %>%
      pivot_longer(cols = -index) %>%
      mutate(type = "average of residues from ~40K candidate regions")


    max_index <- max(dat_toplot$index, na.rm = TRUE)

    p <- ggplot2::ggplot(dat_toplot %>%
                           group_by(index, metric) %>%
                           summarise(across(everything(), mean, na.rm = T))) +
      #ggplot2::geom_violin(aes(x = name, y = value), trim = TRUE) +
      ggplot2::geom_point(aes(x = index, y = value)) +
      ylab("model sensitivity") +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5)) +
      ggplot2::facet_grid(rows = vars(metric), scales = "free_y") +
      theme_bw() +
      theme(
        strip.background = element_blank(),     # removes grey box
        strip.text = element_text(
          color = "black",
          size = 10
        )  )


    p2 <- ggplot2::ggplot(data = dat_toplot3, aes(x = index, y = value)) +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 1, mapping = aes(shape = type)
      )

    if(names(met_it)[y] == "core") {
      p2 <- p2 +
        stat_summary(
          fun.min = ~ quantile(.x, 0.25),
          fun.max = ~ quantile(.x, 0.75),
          geom = "errorbar",
          width = 0.2, linewidth = 0.2
        )
    }

    p2 <- p2 +
      ylab("input data") +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5)) +
      ggplot2::geom_point(data = dat_toplot2, aes(x = index, y = value, shape = type), size = 2) +
      ggplot2::facet_grid(rows = vars(name), scales = "free_y") +
      scale_shape_manual(
        values = c(
          "average of residues within or adjacent to known GPCR ligands" = 21,
          "average of residues from ~40K candidate regions" = 15
        )
      ) +
      theme_bw() +
      theme(
        strip.background = element_blank(),     # removes grey box
        strip.text = element_text(
          color = "black",
          size = 10
        )  )

    nn_input_comb[nn_input_comb$peps == "ANO8_w12-47", "pep_id"] <- "ANO8"

    seq_dat <- nn_input_comb %>%
      filter((target %in% targ & known == 1) | peps == "ANO8_w12-47") %>%
      mutate(meta_data = map2(meta_data, pep_id, \(x,y) x %>%
                                mutate(metric = y) %>%
                                mutate(index_og = index) %>%
                                mutate(index = row_number()))) %>%
      pull(meta_data) %>%
      bind_rows


    aa_mat <- seq_dat %>%
      select(AA, metric) %>%
      group_by(metric) %>%
      mutate(pos = row_number()) %>%
      ungroup() %>%
      pivot_wider(names_from = metric, values_from = AA, id_cols = pos) %>%
      mutate(across(-pos, ~replace_na(., "X"))) %>%
      select(-pos) %>%
      as.matrix() %>%
      t()

    D <- ape::dist.aa(aa_mat)
    tre <- ape::nj(D)
    tre <- ape::ladderize(tre)


    tree_plot <- ggtree::ggtree(tre, layout = "roundrect")

    msa_ord <- ggtree::get_taxa_name(tree_plot)

    tree_plot <- tree_plot +
      ggtree::geom_tree(layout = "roundrect") +
      theme(legend.position = "top")

    msa_ord <- rev(ggtree::get_taxa_name(tree_plot))




    seq_dat$metric <- factor(seq_dat$metric,
                             levels = msa_ord)

    p3 <- ggplot2::ggplot(seq_dat) +
      ggplot2::geom_tile(mapping = aes(x = index, y = metric, fill = known_idx)) +
      ggiraph::geom_point_interactive(data = seq_dat %>%
                                        mutate(tt_value = "test"),
                                      aes(x = index, y = metric, tooltip = tt_value, data_id = index), pch = 15, size = 2.5, color = "grey85") +

      scale_fill_discrete(palette = ggsci::pal_simpsons()) +
      ylab("") +
      xlab("") +

      scale_x_continuous(expand = expansion(mult = c(0, 0)),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5)) +
      ggplot2::geom_text(aes(x = index, y = metric, label = AA), size = 1.8, fontface = "bold", color = "black") +
      theme_bw()

    p3 <- patchwork::wrap_plots(tree_plot, p3)



    p_final[[targ[1]]] <- list(p3, p2, p)

    ggsave(filename = paste0("~/AF2_analysis/shap_1ddcnnddfccdddd_", names(met_it)[y], targs[[x]][1], ".svg"), p3, width = 16, height = 12)


  }

  p_final <- unlist(p_final, recursive = FALSE)

  p_final2 <- Reduce(`+`, p_final) +  patchwork::plot_layout(ncol = 2, nrow = 3, byrow = FALSE, heights = c(1.4,1,1),
                                                             guides = "collect") &
    theme(legend.position = "top", legend.title = element_blank())

  ggsave(filename = paste0("~/AF2_analysis/shap_1dcnndddfccdddd_", names(met_it)[y], ".svg"), p_final2, width = 16, height = 16)


}

















library(fastshap)

pred_fun <- function(object, newdata) {

  n <- nrow(newdata)

  reshaped <- array(
    newdata,
    dim = c(n, 36, 110)
  )

  preds <- predict(object, reshaped)

  as.numeric(preds$global[,1])
}


x_val_flat <- matrix(
  x_val,
  nrow = dim(x_val)[1],
  ncol = prod(dim(x_val)[-1])
)

baseline <- mean(pred_fun(models[["C"]], x_val_flat))

shap_values <- fastshap::explain(
  object = models[["C"]],
  X = x_val_flat,
  pred_wrapper = pred_fun,
  nsim = 100,
  baseline = baseline
)































