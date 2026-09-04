




## =============================================================================
## Vanilla-gradient saliency for the 1D-CNN global head, + the matching input
## profiles, sequence/known_idx track and NJ tree.
##
## Run AFTER 10_1dcnn_new6.R, in the same session. Needs: models, nn_input,
## nn_input_comb, known_dat, c_dat, all_params3, seq_len, n_channels.
## =============================================================================

library(keras3); library(tensorflow); library(tidyverse)

.need_obj <- c("models", "nn_input", "nn_input_comb", "known_dat", "c_dat",
               "all_params3", "seq_len", "n_channels")
.miss_obj <- .need_obj[!vapply(.need_obj, exists, logical(1))]
if (length(.miss_obj))
  stop("10_2_model_importance: run 10_1dcnn_new6.R first -- missing: ",
       paste(.miss_obj, collapse = ", "))

.need_pkg <- c("ape", "ggtree", "ggiraph", "ggsci", "patchwork")
.miss_pkg <- .need_pkg[!vapply(.need_pkg, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss_pkg))
  stop("10_2_model_importance: install.packages(c('",
       paste(.miss_pkg, collapse = "','"), "'))")

## output naming. Stamped so successive runs don't overwrite each other.
out_prefix <- "model_importance"
stamped <- function(..., ext = ".svg", dir = "~/AF2_analysis")
  file.path(dir, paste0(paste0(...), "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ext))

## extra (non-known) windows to show in the sequence track alongside the knowns;
## character(0) for knowns only. Ids exactly as they appear in nn_input_comb$peps.
extra_peps <- c("NUCB1_w45-80", "NUCB2_w48-83")

## rebuild the (n, seq_len, n_channels) input array from a dataset list-column.
## Replaces the per-model nn_in_* objects that 10_1dcnn_new6.R no longer creates.
build_x <- function(dataset) {
  n <- length(dataset[["data"]])
  aperm(array(unlist(dataset[["data"]]), dim = c(seq_len, n_channels, n)), c(3, 1, 2))
}

## Gradient saliency on the global head.
##
## `use_logit = TRUE` differentiates the PRE-ACTIVATION logit instead of the
## sigmoid probability. The global head ends in a sigmoid, whose slope p*(1-p)
## collapses toward 0 for confidently-scored windows -- and this script slices
## known positives, which are exactly the saturated ones. Differentiating the
## sigmoid therefore shrinks gradients for the windows the model is most sure
## about, independent of how important a channel actually is. The logit has no
## such ceiling, so magnitudes stay comparable across windows.
##
## The logit is rebuilt as h %*% W + b from the final Dense layer's own weights
## and the activation it feeds on, so nothing is recomputed by hand and the
## sigmoid is simply never applied. Signs are unchanged: sigmoid is monotonic,
## so d(logit)/dx and d(p)/dx always share a sign (they differ by the strictly
## positive factor 1/(p*(1-p))).
saliency_fn <- function(model, x_batch, use_logit = TRUE) {
  x_tensor <- tf$convert_to_tensor(x_batch, dtype = tf$float32)

  if (use_logit) {
    g         <- keras3::get_layer(model, "global")
    pen_model <- keras3::keras_model(model$input, g$input)   # activation feeding the head
    W <- tf$convert_to_tensor(g$kernel)
    b <- tf$convert_to_tensor(g$bias)
  }

  with(tf$GradientTape() %as% tape, {
    tape$watch(x_tensor)
    out <- if (use_logit) tf$matmul(pen_model(x_tensor), W) + b
           else           model(x_tensor)[[1]]
  })

  grads <- tape$gradient(out, x_tensor)
  as.array(grads)
}


targs <- list(
  c("C", "loop_C"),
  c("N", "loop_N")
)


## Channel groups to profile. Built by intersecting with the ACTIVE all_params3
## so a renamed/removed channel drops out loudly instead of silently yielding an
## empty slice (the old hard-coded one-hot AA_* groups did exactly that once the
## amino acids were re-encoded as properties).
met_it <- list(
  core   = c("cons_rs_n", "min_afm", "mean_afm", "relASA"),
  aa     = c("AA_hydro", "AA_charge", "AA_mw", "AA_pI"),
  angles = c("Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin"),
  energy = c("NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy"),
  SS     = c("SS_G", "SS_B", "SS_H")
)
met_it <- lapply(met_it, intersect, all_params3)
if (any(lengths(met_it) == 0))
  message("10_2_model_importance: dropping group(s) absent from all_params3: ",
          paste(names(met_it)[lengths(met_it) == 0], collapse = ", "))
met_it <- met_it[lengths(met_it) > 0]

for(y in seq_along(met_it)) {

  mets <- met_it[[y]]
  p_final <- list()


  for(x in seq_along(targs)) {

    targ <- targs[[x]]

    ## saliency for THIS terminus's own validation windows (previously always
    ## read the C-terminus object, so the N pass was scored on C data).
    val_x    <- build_x(nn_input[[targ[1]]]$val)
    val_ypos <- which(nn_input[[targ[1]]]$val$known == 1)
    if (length(val_ypos) == 0) { message("no val positives for ", targ[1], " -- skipping"); next }

    shap_vals <- saliency_fn(models[[targ[1]]], val_x)


    ## matrix(): keeps the slice 2-d (n x seq_len) even when the val set has a
    ## single positive, which the db-only filter makes entirely possible.
    dat_toplot <- lapply(mets, \(x)
        matrix(shap_vals[val_ypos, , which(all_params3 == x)], nrow = length(val_ypos)) %>%
          as_tibble(.name_repair = ~ paste0("V", seq_along(.))) %>% mutate(metric = x)) %>%
      bind_rows() %>%
      pivot_longer(cols = -metric) %>%
      mutate(index = setNames(seq_len(seq_len), paste0("V", seq_len(seq_len)))[name]) %>%
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
                           summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))) +
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

    ## label the extra windows by their window id (they have no pep_id)
    for (.ep in extra_peps)
      nn_input_comb[nn_input_comb$peps == .ep, "pep_id"] <- sub("_w.*$", "", .ep)


    seq_dat <- nn_input_comb %>%
      filter((target %in% targ & known == 1) | peps %in% extra_peps) %>%
      mutate(meta_data = map2(meta_data, pep_id, \(x,y) x %>%
                                mutate(metric = y) %>%
                                mutate(index_og = index) %>%
                                mutate(index = row_number()))) %>%
      pull(meta_data) %>%
      bind_rows


    aa_mat <- seq_dat %>%
      filter(index > 7) %>%
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



    ## wrap_elements(): p3 is itself a 2-plot patchwork (tree + sequence track).
    ## Left bare, `Reduce(`+`)` flattens the FIRST one into 2 top-level panels but
    ## nests later ones as 1, so the combined figure came to 7 panels against a
    ## 2x3 layout. Wrapping pins it to exactly one cell per terminus, giving the
    ## intended 2 columns (C, N) x 3 rows (tree+seq, input, saliency).
    p_final[[targ[1]]] <- list(patchwork::wrap_elements(p3), p2, p)

    ggsave(filename = stamped(out_prefix, "_", names(met_it)[y], "_", targ[1]), p3, width = 16, height = 12)


  }

  p_final <- unlist(p_final, recursive = FALSE)

  p_final2 <- Reduce(`+`, p_final) +  patchwork::plot_layout(ncol = 2, nrow = 3, byrow = FALSE, heights = c(1.4,1,1),
                                                             guides = "collect") &
    theme(legend.position = "top", legend.title = element_blank())

  ggsave(filename = stamped(out_prefix, "_", names(met_it)[y]), p_final2, width = 16, height = 16)


}

















## ---- optional: KernelSHAP over the flattened window ------------------------
## OFF by default -- nsim=100 over seq_len*n_channels features is very slow.
## Previously broken: `x_val` was never defined and the reshape hard-coded 110
## channels, which no longer matches n_channels.
run_shap  <- FALSE
shap_term <- "C"

if (run_shap) {
  library(fastshap)

  pred_fun <- function(object, newdata) {
    reshaped <- array(as.matrix(newdata), dim = c(nrow(newdata), seq_len, n_channels))
    as.numeric(predict(object, reshaped)$global[, 1])
  }

  x_val      <- build_x(nn_input[[shap_term]]$val)
  x_val_flat <- matrix(x_val, nrow = dim(x_val)[1], ncol = prod(dim(x_val)[-1]))

  baseline <- mean(pred_fun(models[[shap_term]], x_val_flat))

  shap_values <- fastshap::explain(
    object       = models[[shap_term]],
    X            = x_val_flat,
    pred_wrapper = pred_fun,
    nsim         = 100,
    baseline     = baseline
  )
}































