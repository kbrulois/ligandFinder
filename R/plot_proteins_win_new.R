# Forked variant of make_protein_plot that overlays NN-derived peptide windows
# from a `new_nn_input` tibble (e.g. `nn_input_comb` produced by 10_1dcnn_new5.R).
#
# Each window renders as a stacked horizontal segment (one segment per overlap
# layer, capped at 3); hover shows the overall score; click opens a per-window
# smoothed prediction-score plot whose x-axis is in protein-residue coordinates
# (mapped via meta_data$index_og).
#
# Helpers (smoother_func, assemble_anno_feats, expand_by_residue, ...) and
# package data (anno_feats, pep_nudges, species_dat, sim_mats, all_mets,
# all_desc_colors, desc_colors, tt_lut, DBC, desc_vars, id_map) are reused from
# plot_proteins.R.

nn_class_cols <- setNames(
  c("#FED439FF", "#370335FF", "#8A9197FF", "#D2AF81FF",
    "#D5E4A2FF", "#197EC0FF", "grey85", "#075149FF"),
  c("CT_cleavage_context", "DB", "gap", "NT_cleavage_context",
    "pep_other", "pep_pocket", "padding", "none")
)

nn_win_target_cols <- c(
  "db_N"       = "#197EC0FF",
  "db_C"       = "#C80813FF",
  "chym_N"     = "#46732EFF",
  "chym_C"     = "#FD7446FF",
  "pep_end_N"  = "#370335FF",
  "pep_end_C"  = "#8A9197FF"
)

## plot.margin = c(top, right, bottom, left), in lines.
## These are now IDENTICAL. The strip used to carry an extra 2 lines on the left
## to hand-compensate for the y-axis gutter the detail panels have and it did
## not -- fragile, since that gutter's true width depends on the tick labels.
## Instead the strip renders an invisible y title + invisible tick labels of the
## same width (see nn_axis_lab), so both reserve the same gutter and the equal
## margins line the two panels up exactly.
nn_det_margs  <- c(0, 1.6, 2, 16.6)
nn_det_margs2 <- c(0, 1.6, 2, 16.6)

## Hover text for the y-axis row labels in the main track. DRAFT -- inferred
## from how each metric is computed in the pipeline, not from documentation, so
## correct anything that is wrong. Keys must match the rendered label exactly.
nn_label_tt <- c(
  ## conservation (Aminode alignments)
  "cons_rs"       = "Aminode evolutionary rate score at this residue -- lower = more conserved.",
  "blos_wt_all"   = "BLOSUM62 similarity to the human residue, weighted, averaged over all aligned species.",
  "blos_wt_all_n" = "blos_wt_all scaled by species_limit (the number of species actually aligned here).",
  "blos_wt_mam"   = "BLOSUM62 similarity to the human residue, weighted, mammals only.",
  "gran_wt_all"   = "Grantham distance (physicochemical dissimilarity) to the human residue, weighted, all species.",
  ## AlphaMissense
  "mean_afm"      = "Mean AlphaMissense pathogenicity across all substitutions at this residue.",
  "min_afm"       = "Minimum AlphaMissense pathogenicity across all substitutions at this residue.",
  ## structure
  "relASA"        = "Relative solvent-accessible surface area from DSSP: 0 = buried, 1 = fully exposed.",
  "SS"            = "DSSP secondary structure assignment.",
  "topo"          = "Membrane topology; 'e' marks extracellular.",
  ## model scores
  ## The "c" in nn4c/xgb4c means context: 10_score_AA.R replicates every base
  ## feature at _lag1/_lag2/_lead1/_lead2, so both models see i-2 .. i+2 -- a
  ## 5-residue receptive field, against the window model's 36.
  "pep_nn4c_s6"   = "Residue-level MLP prediction, 5-residue receptive field (nn4c), smoothed over 6 residues.",
  "pep_xgb4c_s6"  = "Residue-level XGBoost prediction, 5-residue receptive field (xgb4c), smoothed over 6 residues.",
  ## UniProt annotation rows
  "modification"  = "UniProt post-translational modification at this residue.",
  "SV"            = "UniProt sequence variant (dbSNP identifier where one exists).",
  "domain"        = "UniProt domain annotation covering this residue."
)

## bottom_p's right margin, in points -- the `&` theme applied to the seq_p /
## main_p patchwork. The detail panels and the NN strip MUST use the same value:
## det_left_in already pins their left edges to the main track's, so any
## difference in the right margin makes their panels a different width, and the
## residue positions then drift apart toward the C-terminus (1.6 lines = 17.6pt
## against 10pt was ~0.74 residues of lag at p_width = num_res/7).
nn_marg_r_pt <- 10

## Top margin for seq_p, in points. Its residue-index labels are drawn above the
## panel (geom_text vjust = -2), so with t = 0 they fall outside the SVG viewport
## and get clipped by its top edge -- which reads as the detail panel above
## cutting them off. seq_p sets margin(1, 0, 0, 0) itself, but the `&` theme on
## the patchwork overrides it, so the room has to be made here.
nn_seq_top_pt <- 8

## half-height of a window rect, in `layer` units
nn_win_half_h <- 0.17

## Shared y-tick formatter. Fixed width (5 chars, e.g. " 0.25") so the reserved
## gutter is identical in the detail panel and the strip no matter what range
## the data happens to span.
nn_axis_lab <- function(b) formatC(b, width = 5, format = "f", digits = 2)

## Absolute width, in inches, to the LEFT of a plot's panel (margin + axis title
## + tick labels + ticks). The NN strip, the detail popups and the main residue
## track are separate girafe widgets rendered at the same width_svg, so they line
## up iff this value matches. Handles patchwork (the main track is seq_p / plot).
panel_left_in <- function(p) {
  g <- if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)
  lay <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  if (!nrow(lay)) return(NA_real_)
  l <- min(lay$l)
  if (l <= 1) return(0)
  sum(grid::convertWidth(g$widths[seq_len(l - 1)], "in", valueOnly = TRUE))
}


## Mirror of panel_left_in for the RIGHT side: absolute width, in inches, from
## the panel's right edge to the plot's. The main track's panel sits ~20pt from
## its SVG edge -- the 10pt subplot margin plus roughly as much again that the
## patchwork inserts -- so matching only the subplot margin still leaves the
## detail/strip panels wider than the main track's, and residue positions drift
## apart toward the C-terminus.
panel_right_in <- function(p) {
  g <- if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)
  lay <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  if (!nrow(lay)) return(NA_real_)
  r <- max(lay$r)
  if (r >= length(g$widths)) return(0)
  sum(grid::convertWidth(g$widths[seq(r + 1, length(g$widths))], "in", valueOnly = TRUE))
}


## Inline the shared JS/CSS that save_html() writes to libdir and references
## relatively -- without which a page is broken unless dependency_files/ travels
## beside it. Everything ggiraph pulls in is plain text (~288K total), so this is
## a straight substitution: no base64, and no pandoc, which the usual
## rmarkdown::pandoc_self_contained_html() route would require and which is not
## installed here. Replacements are fixed-string on the exact emitted tags.
nn_inline_deps <- function(html_file, lib_dir) {
  if (!file.exists(html_file) || !dir.exists(lib_dir)) return(invisible(FALSE))
  h <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  base <- fs::path_file(lib_dir)
  files <- fs::dir_ls(lib_dir, recurse = TRUE, type = "file",
                      glob = "*.js") |> c(fs::dir_ls(lib_dir, recurse = TRUE,
                                                     type = "file", glob = "*.css"))
  n <- 0L
  for (f in files) {
    rel <- paste0(base, "/", fs::path_rel(f, lib_dir))
    body <- paste(readLines(f, warn = FALSE), collapse = "\n")
    ## a literal </script> inside inlined JS would close the tag early
    body <- gsub("</script>", "<\\/script>", body, fixed = TRUE)
    if (grepl("\\.js$", f)) {
      tag <- paste0('<script src="', rel, '">')
      if (grepl(tag, h, fixed = TRUE)) {
        h <- sub(tag, paste0("<script>\n", body, "\n"), h, fixed = TRUE); n <- n + 1L
      }
    } else {
      tag <- paste0('<link href="', rel, '" rel="stylesheet" />')
      if (grepl(tag, h, fixed = TRUE)) {
        h <- sub(tag, paste0("<style>\n", body, "\n</style>"), h, fixed = TRUE); n <- n + 1L
      }
    }
  }
  writeLines(h, html_file)
  invisible(n)
}


assign_overlap_layers <- function(start, end, gap = 0) {
  out <- integer(length(start))
  layer_end <- numeric(0)
  for (i in order(start, end)) {
    k <- which(layer_end < start[i] - gap)   # rows already clear of this window
    if (length(k)) {
      k <- k[1]                              # lowest such row
      out[i] <- k
      layer_end[k] <- end[i]
    } else {
      layer_end <- c(layer_end, end[i])
      out[i] <- length(layer_end)
    }
  }
  out
}

make_detail_panel <- function(per_index, meta_data, title = NULL,
                              x_range = NULL, seq_offset = -0.45,
                              win_start = NULL, win_end = NULL,
                              left_in = NULL, right_in = NULL) {
  stopifnot(nrow(per_index) == nrow(meta_data))

  joined <- per_index %>%
    dplyr::select(-dplyr::any_of("index")) %>%
    dplyr::bind_cols(
      meta_data %>% dplyr::select(index_og = index, AA)
    ) %>%
    tidyr::pivot_longer(cols = -c(index_og, AA),
                        names_to = "class", values_to = "value") %>%
    dplyr::mutate(class = factor(class, levels = names(nn_class_cols)))

  if (is.null(x_range)) x_range <- range(joined$index_og, na.rm = TRUE)

  legend_pos <- "top"
  title_hjust <- 0.5
  if (!is.null(win_start) && !is.null(win_end) &&
      diff(x_range) > 0) {
    win_center <- (win_start + win_end) / 2
    rel <- (win_center - x_range[1]) / diff(x_range)
    rel <- max(0.02, min(0.98, rel))
    legend_pos <- c(rel, 1.02)
    title_hjust <- rel
  }
  if (!isTRUE(getOption("lf.nn_center_title", TRUE))) title_hjust <- 0

  p <- ggplot2::ggplot(joined, ggplot2::aes(x = index_og, y = value, color = class)) +
    ggplot2::geom_smooth(method = "loess", span = 0.6, se = FALSE,
                         linewidth = 0.6, na.rm = TRUE) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(tooltip = sprintf("res %d (%s)\n%s: %.3f",
                                     index_og, AA, class, value),
                   data_id = as.character(index_og)),
      pch = 21, stroke = 0.6, size = 1.6
    ) +
    ggplot2::scale_color_manual(values = nn_class_cols) +
    ggplot2::scale_x_continuous(
      limits = x_range,
      expand = ggplot2::expansion(add = c(seq_offset, -0.4)),
      breaks = seq(0, x_range[2], by = 10)
    ) +
    ggplot2::scale_y_continuous(labels = nn_axis_lab) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ## Softmax over the 8 classes: the values at each residue sum to 1, so
    ## "probability" is accurate here -- unlike the window score in the title,
    ## which is a raw uncalibrated sigmoid.
    ggplot2::labs(title = title, x = NULL, y = "per-residue class probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   ## Centred on the legend, which floats at x = rel -- the
                   ## window's centre -- just above the panel, so the title tracks
                   ## the window instead of the panel. Falls back to the panel
                   ## centre when no window is given.
                   ## options(lf.nn_center_title = FALSE) restores hjust = 0.
                   plot.title = ggplot2::element_text(hjust = title_hjust),
                   legend.position = legend_pos,
                   legend.justification = "center",
                   legend.direction = "horizontal",
                   legend.background = ggplot2::element_rect(fill = "white", color = NA),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.margin = if (is.null(left_in))
                       unit(nn_det_margs, "lines")
                     else
                       ggplot2::margin(t = nn_det_margs[1], b = nn_det_margs[3],
                                       unit = "lines") +
                       ggplot2::margin(r = nn_marg_r_pt, unit = "pt") +
                       ggplot2::margin(l = left_in, unit = "in") +
                       ggplot2::margin(r = if (is.null(right_in)) 0 else right_in,
                                       unit = "in"))
  p
}

make_protein_plot_win <- function(old_nn_input,
                                  new_nn_input,
                                  plot_dir,
                                  pep_tp) {

  tryCatch({

    ## ---- progress logging ------------------------------------------------
    ## Prints elapsed seconds at each stage so a stall can be located. Silence
    ## with options(lf.plot_verbose = FALSE).
    .t0 <- Sys.time()
    .step <- function(msg) {
      if (!isTRUE(getOption("lf.plot_verbose", TRUE))) return(invisible())
      message(sprintf("  [%-8s %6.1fs] %s",
                      as.character(old_nn_input[["gene"]])[1],
                      as.numeric(difftime(Sys.time(), .t0, units = "secs")), msg))
      utils::flush.console()
    }
    .step("start")

    # Normalize list-columns: bind_cols/pivot_longer downstream choke on
    # grouped/rowwise tibbles nested in per-residue list-cols (e.g. aa_scores,
    # cons, dssp, af_missense, alignment_AA). Force every nested data.frame to a
    # plain tibble so the copied make_protein_plot body behaves as it does there.
    old_nn_input <- old_nn_input %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(is.list),
        \(col) purrr::map(col, \(d) {
          if (is.data.frame(d)) tibble::as_tibble(dplyr::ungroup(d)) else d
        })
      ))

    # ----- prepare NN windows for this protein -----
    gene_tp <- old_nn_input[["gene"]]

    nn_anno <- NULL
    if (!is.null(new_nn_input) && nrow(new_nn_input) > 0) {
      nn_w <- new_nn_input %>%
        dplyr::filter(gene == gene_tp) %>%
        dplyr::filter(!is.na(pred))

      if (nrow(nn_w) > 0) {
        coords <- stringr::str_match(nn_w$peps, "_w(\\d+)-(\\d+)$")
        nn_w <- nn_w %>%
          dplyr::mutate(start = as.integer(coords[, 2]),
                        end   = as.integer(coords[, 3])) %>%
          dplyr::filter(!is.na(start), !is.na(end)) %>%
          dplyr::arrange(dplyr::desc(pred))

        if (nrow(nn_w) > 0) {
          nn_w$layer <- assign_overlap_layers(nn_w$start, nn_w$end)
          # tolerate older nn_input_comb that lack the newer columns
          for (col in c("end_type", "pred_raw", "rank_cat",
                        "nn_closest_peptide", "nn_closest_sim")) {
            if (!col %in% names(nn_w)) nn_w[[col]] <- NA
          }
          # segment colour is pred_raw; on older inputs that lack it every window
          # would otherwise render as na.value grey, so fall back to `pred`.
          if (all(is.na(nn_w$pred_raw))) {
            message("plot_proteins_win_new: no pred_raw for ", gene_tp,
                    " -- colouring windows by `pred` instead")
            nn_w$pred_raw <- nn_w$pred
          }
          nn_w <- nn_w %>%
            dplyr::mutate(
              target_short = stringr::str_replace(as.character(target), "^loop_", ""),
              wt_combo     = paste0(win_type, "_", target_short),
              panel_id     = sprintf("nnwin_%d", dplyr::row_number()),
              feature      = wt_combo,
              # single-line label above each bar: window range (no gene / no "w"
              # prefix), category rank, raw score. win_type and the nearest known
              # peptide stay in the tooltip.
              label_txt    = sprintf("%d-%d  rank: %s  score: %.2f",
                                     start, end, as.character(rank_cat), pred_raw),
              tooltip      = sprintf(paste0("%s\nrank_cat: %s\npred: %.3f  (raw %.3f)\n",
                                            "nn: %s  (sim %.2f)\ntype: %s  end_type: %s  terminus: %s"),
                                     peps, as.character(rank_cat), pred, pred_raw,
                                     dplyr::coalesce(as.character(nn_closest_peptide), "NA"),
                                     nn_closest_sim, win_type,
                                     dplyr::coalesce(as.character(end_type), "NA"),
                                     target_short)
            )
          nn_anno <- nn_w
        }
      }
    }

    ## sim_mats and species_dat ship in inst/extdata but are NOT package data and
    ## NOT in the plot bundle -- save_plot_bundle.R keeps species_dat and drops
    ## sim_mats, and plot_only.R loads neither. Without sim_mats the per-species
    ## alignment block below throws, its tryCatch returns NULL, and the species
    ## rows vanish with no message. Fall back to the shipped copies.
    if (!exists("sim_mats", inherits = TRUE))
      sim_mats <- readRDS(system.file("extdata", "sim_mats.rds",
                                      package = "ligandFinder"))
    if (!exists("species_dat", inherits = TRUE))
      species_dat <- readRDS(system.file("extdata", "species_dat.rds",
                                         package = "ligandFinder"))

    # ----- begin near-verbatim copy of make_protein_plot body -----
    .step(sprintf("nn windows prepared (%d)", if (is.null(nn_anno)) 0L else nrow(nn_anno)))
    to_plot <- expand_by_residue(old_nn_input)

    to_plot <- to_plot %>%
      group_by(gene) %>%
      mutate(max = max(index)) %>%
      mutate(across(starts_with(c("pep_nn", "pep_xgb", "chem_nn", "chem_xgb")),
                    .fns = ~smoother_func(x = ., append_name = "s"),
                    .unpack = TRUE))

    ## all_mets is a session object (bundled by save_plot_bundle.R) whose score
    ## names are hardcoded in 10_1dcnn_new6.R and carry a model suffix that
    ## 10_score_AA_xgboost.R derives from nn[["neural_net"]][8] -- so it drifts
    ## whenever the model list changes, and when it stops matching the score
    ## tracks disappear with no error. Union it with the score columns
    ## expand_by_residue actually produced, so a renamed model still plots.
    score_cols <- grep("^(pep|chem)_(nn|xgb)[^_]*(_s[0-9]+)?$",
                       names(to_plot), value = TRUE)
    all_mets_use <- if (exists("all_mets", inherits = TRUE))
      union(unname(all_mets), score_cols) else score_cols
    .step(sprintf("score tracks (%d): %s", length(score_cols),
                  if (length(score_cols)) paste(score_cols, collapse = ", ")
                  else "(none in to_plot)"))

    transfrom_aligments <- function(x, vals_to = "AA") {
      x %>%
        mutate(index = row_number()) %>%
        pivot_longer(cols = -index, names_to = "metric", values_to = vals_to)
    }

    .step("expand_by_residue done")
    to_plot_aln <- tryCatch({old_nn_input %>%
        mutate(sim_mat = map(alignment_AA, \(x) {
          ref_seq <- x[["Homo_sapiens"]]
          Map(\(sm) {
            x %>%
              mutate(across(everything(), \(y) {
                sm[cbind(ref_seq, y)]
              }))
          }, sim_mats)
        })) %>%
        mutate(alignment_AA = map(alignment_AA, ~transfrom_aligments(.))) %>%
        mutate(sim_mat = map(sim_mat, \(x) {

          tmp <- map(names(x), ~transfrom_aligments(x[[.]], vals_to = .))

          tmp <- purrr::reduce(tmp, function(x, y) {
            bind_cols(x, y %>% select(!any_of(names(x))))
          })

          bind_rows(
            species_dat %>%
              mutate(metric = as.character(aminode)) %>%
              mutate(index = 0) %>%
              mutate(blos = myo_sim,
                     gran = myo_sim) %>%
              select(index, metric, blos, gran) %>%
              filter(metric %in% unique(tmp[["metric"]])),
            tmp)

        })) %>%
        mutate(alignment_final = map2(alignment_AA, sim_mat, ~right_join(.x, .y, by = c("index", "metric")))) %>%
        select(gene, alignment_final) %>%
        unnest(alignment_final)
    }, error = function(e) {
      ## never swallow this silently: a missing sim_mats or a shape change in
      ## alignment_AA both land here and just delete the species rows.
      .step(paste0("species alignment SKIPPED: ", conditionMessage(e)))
      NULL
    })

    to_plot_c <- to_plot %>%
      pivot_longer(cols = any_of(all_mets_use),
                   names_to = "metric",
                   values_to = "value")

    ## DBC (the NTC / CTC cleavage rows) deliberately omitted -- the dibasic
    ## anchor is now shown on the window rects in the NN strip instead.
    desc_vars <- c(
                   "modification",
                   "SV",
                   "domain",
                   "topo",
                   "SS") %>% rev

    to_plot_desc <- to_plot %>%
      pivot_longer(cols = any_of(desc_vars),
                   names_to = "metric",
                   values_to = "value_desc")

    df <- bind_rows(to_plot_c,
                    to_plot_aln,
                    to_plot_desc)

    mets_in_plot <- unique(df[["metric"]])

    met_order <- c(all_mets_use, desc_vars, levels(species_dat[["aminode"]]))

    mets_in_plot <- mets_in_plot[match(met_order, mets_in_plot)] %>% .[!is.na(.)]

    df <- df %>%
      mutate(metric_type = case_when(metric %in% species_dat[["aminode"]] ~ "blosum62\n-------------\ngrantham",
                                     metric %in% desc_vars ~ "discrete",
                                     TRUE ~ "")) %>%
      mutate(metric_type = factor(metric_type, levels = c("", "discrete", "blosum62\n-------------\ngrantham"))) %>%
      mutate(metric = factor(metric,
                             levels = mets_in_plot))

    y_axis_colors <- tibble(color = c("#709AE1FF", "#FD7446FF", "#46732EFF", "#370335FF", "#8A9197FF"),
                            class = c("Mammalia", "Aves", "Lepidosauria", "Actinopteri", "Amphibia")
    )

    y_axis_colors <- left_join(species_dat, y_axis_colors, by = "class") %>%
      mutate(metric = as.character(aminode)) %>%
      select(color, metric)

    y_axis_colors <- bind_rows(tibble(metric = mets_in_plot,
                                      color = "black"),
                               y_axis_colors)

    y_axis_colors <- setNames(y_axis_colors[["color"]],
                              y_axis_colors[["metric"]])

    cons_dat <- df %>% filter(metric_type == "blosum62\n-------------\ngrantham")

    if(nrow(cons_dat) > 0) {seq_offset <- 0.6} else {seq_offset <- -0.45}

    ## make_cm_script_text() expects the OUTER pep_tp (it does pull(data)[[1]]
    ## itself), so keep a copy before this reassignment.
    pep_tp_outer <- pep_tp
    pep_tp <- pep_tp %>% pull(data) %>% `[[`(1)

    max_index <- max(df$index, na.rm = TRUE)
    ind_tp <- c(1, seq(from = 0, to = max_index, by = 10))

    features <- old_nn_input %>% pull(features) %>% `[[`(1)
    AA_sequence <- old_nn_input %>% pull(sequence_uni) %>% stringr::str_split(., "", simplify = TRUE) %>% `c`

    seq_dat <- bind_rows(tibble(AA = NA, index = 0, index_tp = NA),
                         tibble(AA = AA_sequence) %>%
                           mutate(index = row_number()) %>%
                           mutate(index_tp = if_else(index %in% ind_tp, index, NA))
    )

    anno_feat_dat <- assemble_anno_feats(features, peps = pep_tp)

    anno1 <- old_nn_input %>%
      pull(annotations) %>%
      `[[`(1) %>%
      filter(annotation_name == "comment") %>%
      filter(annotation_type %in% c("subcellular location", "tissue specificity", "disease")) %>%
      group_by(annotation_type) %>%
      summarise(annotation = list(paste0(annotation, collapse = "; "))) %>%
      pivot_wider(names_from = annotation_type, values_from = annotation)

    anno2 <- old_nn_input %>%
      pull(annotations) %>%
      `[[`(1) %>%
      filter(annotation_name == "dbReference" & name_1 == "disease") %>%
      group_by(annotation_type) %>%
      summarise(annotation = paste0(annotation, collapse = "; ")) %>%
      {paste0(.[["annotation_type"]], ": ", .[["annotation"]])}

    p_title <- old_nn_input %>%
      select(accession, gene, files) %>%
      mutate(uniprot_name = setNames(id_map$`Entry Name`, id_map$Entry)[accession]) %>%
      mutate(Aminode = paste0("http://www.aminode.org/?gene=", gene)) %>%
      mutate(UniProt = paste0("https://www.uniprot.org/uniprotkb/", accession, "/entry")) %>%
      mutate(chatGPT = paste0("https://chat.openai.com/?q=What+are+the+known+receptors+for+Gene:+", gene, ",+and+are+any+of+these+receptors+GPCRs?+If+so+list+them+first+and+discuss+their+functions.")) %>%
      mutate(GWAS = paste0("https://www.ebi.ac.uk/gwas/genes/", gene)) %>%
      mutate(PubMed = paste0("https://pubmed.ncbi.nlm.nih.gov/?term=", gene)) %>%
      mutate(Disease = paste0("https://chat.openai.com/?q=What+are+the+known+disease+associations+for+", gene, "+.+Give+PubMed+papers+published+in+the+last+three+years+to+support+your+claims.")) %>%
      mutate(AlphaFoldDB = paste0("https://alphafold.ebi.ac.uk/entry/AF-", accession, "-F1")) %>%
      mutate(`ChimeraX` = paste0("https://stacks.stanford.edu/file/tc396gg4330/v1/", gene, ".cxc"))

    title_text <- paste(paste0("Gene: ", p_title$gene),
                        paste0("Entry: ", p_title$accession),
                        paste0("Entry Name: ", p_title$uniprot_name),
                        sep = "    |    ")

    ## ChimeraX intentionally absent: the "Visualize in ChimeraX" button now
    ## carries the whole script, so a link to a hosted .cxc is redundant.
    all_links <- c("UniProt", "GWAS", "PubMed")

    meta_dat <- tibble(links = list(all_links))

    get_text_width <- function(l) {
      grid::convertWidth(
        grid::grobWidth(grid::textGrob(l, gp = grid::gpar(fontsize = 18))),
        "in",
        valueOnly = TRUE
      )
    }

    meta_dat <- meta_dat %>%
      pivot_longer(everything()) %>%
      unnest(value) %>%
      rowwise() %>%
      mutate(val_len = get_text_width(value)) %>%
      ungroup() %>%
      mutate(index = lag(cumsum(val_len))) %>%
      mutate(index = replace_na(index, 0)) %>%
      mutate(index = index * 10) %>%
      ## spacing of 10 per link, derived from the row count -- the old
      ## seq.default(0, 70, by = 10) hardcoded exactly 8 links and errors on any
      ## other number.
      mutate(index2 = (dplyr::row_number() - 1) * 10) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(tt_value = if(value %in% all_links) {p_title[[value]]} else {NA})

    ## peptides / selection deliberately dropped from the header; anno1's
    ## column names still supply the remaining subtitle labels below.
    subtitle_text <- anno1

    subtitle_text <- lapply(names(subtitle_text), \(x) {paste0(x, ": ", paste0(subtitle_text[[x]][[1]], collapse = "; "))})

    subtitle_text <- paste(subtitle_text, collapse = "\n")

    meta_p <- ggplot2::ggplot(meta_dat) +
      ggiraph::geom_text_interactive(aes(x = index2,
                                         y = name,
                                         label = value,
                                         data_id = value,
                                         onclick = paste0('window.open("', tt_value , '")')), hjust = 0, size = 4) +
      ggplot2::scale_x_continuous(limits = c(0, max(80, max_index)), expand = grid::unit(0, "lines")) +
      ggplot2::scale_y_discrete(expand = grid::unit(0, "lines")) +
      ggplot2::ggtitle(label = title_text, subtitle = subtitle_text) +
      theme_void()

    seq_p <- ggplot2::ggplot(seq_dat) +
      ggiraph::geom_text_interactive(aes(x = index, y = 1,
                                         label = AA,
                                         data_id = index),
                                     size = 3.5
      ) +
      ggplot2::geom_text(data = seq_dat,
                         aes(x = index, y = 1, label = index_tp), vjust = -2, size = 2.5, inherit.aes = FALSE) +
      ggplot2::scale_x_continuous(limits = c(0, (max_index + 1)), expand = expansion(add = c(seq_offset, -0.4))) +
      ggplot2::scale_y_discrete(expand = grid::unit(0, "lines")) +
      theme_void() + theme(
        plot.margin = margin(1, 0, 0, 0),
        panel.spacing = unit(0, "pt")
      )

    .step("data assembled; building main_p")
    main_p <- ggplot2::ggplot(data = df) +

      ggplot2::geom_tile(data = df %>% filter(metric_type == ""),
                         mapping = aes(x = index, y = metric, fill = value))

    if(nrow(cons_dat) > 0) {

      main_p <- main_p +

        ggplot2::geom_tile(data = cons_dat,
                           mapping = aes(x = index, y = metric, fill = gran),
                           width = 1,
                           height = 0.5,
                           position = position_nudge(y = -0.25)) +

        ggplot2::geom_tile(data = cons_dat,
                           mapping = aes(x = index, y = metric, fill = blos),
                           width = 1,
                           height = 0.5,
                           position = position_nudge(y = 0.25)) +

        ## Which half of each species row is which. blos is nudged +0.25 (upper),
        ## gran -0.25 (lower), and the facet strip that would have said so is
        ## blanked by strip.text.y.left = element_blank(). A caption rather than a
        ## positioned geom_text: the conservation facet is the bottom one, so
        ## anything nudged below its last row lands on main_p's x-axis numbers,
        ## and vjust has no value that clears both those and the rows above.
        ggplot2::labs(caption = paste("upper band: BLOSUM62 similarity",
                                      "lower band: Grantham distance",
                                      sep = "     "))

    }

    main_p <- main_p +

      ggiraph::geom_point_interactive(data = df %>%
                                        filter(metric_type != "discrete" & index != 0) %>%
                                        mutate(tt_value = case_when(metric_type == "" ~ paste0(metric, ": ", round(value, 2)),
                                                                    metric_type == "blosum62\n-------------\ngrantham" ~ tryCatch({paste0(metric, "\n", "blosum62: ", round(blos, 2), "\n", "grantham: ", round(gran, 2))}, error = function(e) NA),
                                                                    TRUE ~ NA)),
                                      aes(x = index, y = metric, tooltip = tt_value, data_id = index), pch = 15, size = 2.5, color = "grey85") +

      scale_fill_viridis_c(option = "H", name = "") +

      scale_y_discrete(labels = function(labs) {
        purrr::map_chr(labs, ~ glue::glue("<span style='color:{y_axis_colors[.x]}'>{.x}</span>"))
      }) +

      scale_x_continuous(expand = grid::unit(0, "lines"),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5))

    main_p <- main_p +
      ggplot2::geom_text(data = df %>% filter(metric_type != "discrete"),
                         aes(x = index, y = metric, label = AA), size = 1.8, fontface = "bold", color = "black")

    for(y in desc_vars[desc_vars %in% mets_in_plot]) {

      dat_toplot <- df %>% filter(metric == !!y)

      if(y == "domain") {
        desc_col_scale <- all_desc_colors[[y]]
        all_doms <- unique(dat_toplot[["value_desc"]])
        common_doms <- all_doms[all_doms %in% names(all_desc_colors[[y]])]
        uncommon_doms <- all_doms[!all_doms %in% names(all_desc_colors[[y]])]
        desc_col_scale <- all_desc_colors[[y]][all_desc_colors[[y]] %in% common_doms]
        uncommon_doms <- color_scalify(uncommon_doms, colors = desc_colors[!desc_colors %in% desc_col_scale])
        desc_col_scale <- c(desc_col_scale, uncommon_doms)
      } else if(y == "SV") {

        desc_vals <- unique(dat_toplot[["value_desc"]])
        desc_col_scale <- rep("#709AE1FF", length(desc_vals))
        names(desc_col_scale) <- desc_vals

      } else {
        desc_col_scale <- all_desc_colors[[y]]
      }

      main_p <- main_p +

        ggnewscale::new_scale_fill() +

        ggplot2::geom_tile(data = dat_toplot,
                           mapping = aes(x = index, y = metric, fill = value_desc), show.legend = FALSE) +
        scale_fill_manual(values = desc_col_scale, na.value = "white", name = "")

      if(y %in% c("modification", "topo", "SS")) {

        main_p <- main_p +
          ggiraph::geom_point_interactive(data = dat_toplot %>%
                                            filter(!is.na(value_desc)) %>%
                                            mutate(tt_value = tt_lut[[y]][value_desc]), aes(x = index, y = metric, tooltip = tt_value, data_id = index), pch = 15, size = 2.5, color = "grey85", show.legend = FALSE) +
          ggplot2::geom_text(data = dat_toplot, aes(x = index, y = metric, label = value_desc), size = 1.8, fontface = "bold", color = "black")

      }

      if(y == "domain") {

        main_p <- main_p +
          ggiraph::geom_point_interactive(data = dat_toplot %>%
                                            filter(!is.na(value_desc)) %>%
                                            mutate(tt_value = value_desc), aes(x = index, y = metric, tooltip = tt_value, data_id = index), pch = 15, size = 2.5, color = "grey85", show.legend = FALSE)
      }

      if(y == "SV") {

        dat_toplot <- dat_toplot %>%
          filter(!is.na(value_desc)) %>%
          mutate(SV = "sv") %>%
          mutate(href_text = stringr::str_replace(value_desc, "dbSNP:", "https://www.ncbi.nlm.nih.gov/snp/"))

        main_p <- main_p +
          ggiraph::geom_point_interactive(data = dat_toplot, aes(x = index, y = metric, tooltip = value_desc,
                                                                 onclick = paste0('window.open("', href_text , '")'), data_id = index), pch = 15, size = 2.5, color = "grey85") +
          ggplot2::geom_text(data = dat_toplot,
                             aes(x = index, y = metric, label = SV), size = 1.8, fontface = "bold", color = "black")

      }
    }

    main_p <- main_p +
      ggplot2::facet_grid(rows = vars(metric_type), scales = "free_y", switch = "y", space = "free_y") +

      ggplot2::geom_linerange(
        data = anno_feat_dat %>%
          filter(ggp == "v_rec") %>%
          mutate(start = start - 0.5,
                 end = end + 0.5) %>%
          pivot_longer(cols = c("start", "end")),
        mapping = aes(x = value, ymin = -Inf, ymax = Inf, color = feature), lineend = "round") +

      ggplot2::geom_curve(
        data = anno_feat_dat %>%
          filter(ggp == "v_rec") %>%
          mutate(metric_type = "") %>%
          mutate(metric_type = factor(metric_type, levels = c("", "discrete", "blosum62\n-------------\ngrantham"))),
        mapping = aes(x = start - 0.5, xend = end + 0.5, y = Inf, yend = Inf, color = feature),
        curvature = -1, lineend = "round") +

      ggplot2::geom_curve(
        data = anno_feat_dat %>%
          filter(ggp == "v_rec") %>%
          mutate(metric_type = if(nrow(cons_dat) > 0) {"blosum62\n-------------\ngrantham"} else {"discrete"}) %>%
          mutate(metric_type = factor(metric_type, levels = c("", "discrete", "blosum62\n-------------\ngrantham"))),
        mapping = aes(x = start - 0.5, xend = end + 0.5, y = -Inf, yend = -Inf, color = feature),
        curvature = 1, lineend = "round")

    h_rec_dat <- anno_feat_dat %>%
      filter(ggp == "h_rec") %>%
      mutate(metric = "AA_seq") %>%
      mutate(metric_type = "") %>%
      mutate(metric_type = factor(metric_type, levels = c("", "discrete", "blosum62\n-------------\ngrantham"))) %>%
      mutate(nudge = unname(pep_nudges[feature]))

    if (!"tooltip" %in% names(h_rec_dat)) h_rec_dat$tooltip <- NA_character_
    if (!"onclick" %in% names(h_rec_dat)) h_rec_dat$onclick <- NA_character_
    if (!"data_id" %in% names(h_rec_dat)) h_rec_dat$data_id <- NA_character_

    h_rec_dat <- h_rec_dat %>%
      mutate(row_id = dplyr::row_number()) %>%
      mutate(tt_value = dplyr::coalesce(tooltip, paste0(feature)),
             click_id = dplyr::coalesce(data_id,
                                        paste0("hrec_", feature, "_", row_id)))

    top_y <- df %>% filter(metric_type == "") %>% mutate(metric = droplevels(metric)) %>% pull(metric) %>% levels %>% length

    main_p <- main_p + ggiraph::geom_segment_interactive(
      data = h_rec_dat,
      aes(x = start - 0.3, xend = end + 0.3,
          y = top_y + nudge, yend = top_y + nudge,
          color = feature,
          tooltip = tt_value),
      inherit.aes = FALSE,
      linewidth = 1.2,
      lineend = "round") +

      ggplot2::geom_curve(data = anno_feat_dat %>%
                            filter(ggp == "arch") %>%
                            mutate(metric_type = "discrete") %>%
                            mutate(metric_type = factor(metric_type, levels = c("", "discrete", "blosum62\n-------------\ngrantham"))),
                          aes(x = start, xend = end, y = -Inf, yend = -Inf, color = feature),
                          curvature = -0.3,
                          lineend = "round",
                          inherit.aes = FALSE) +

      scale_color_manual(values = setNames(anno_feats$color, anno_feats$feature), name = "disulfide bond source: ") +

      coord_cartesian(clip = "off") +

      theme_bw() +

      theme(plot.margin = unit(c(0,1,2,1), "lines"),
            axis.title = element_blank(),
            axis.text.y = ggtext::element_markdown(size = 9),
            panel.spacing = grid::unit(0, "lines"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y.left = element_blank(),
            strip.placement = "outside",
            legend.position = "bottom",
            legend.justification = "left",
            legend.ticks = element_line(color = "black", linewidth = 0.1),
            legend.key.height = unit(0.2, "cm"),
            ## set here, not next to labs(): this theme() runs after theme_bw(),
            ## which would otherwise reset it
            plot.caption = element_text(hjust = 0, size = 7, colour = "grey25",
                                        margin = margin(t = 2)))

    p <- list(plot = main_p,
              meta_p = meta_p,
              seq_p = seq_p,
              title_p = title_text,
              subtitle_p = subtitle_text,
              links = p_title)

    num_res = df %>% count(metric, gene) %>% group_by(gene) %>% summarise(max = max(n, na.rm = FALSE)) %>% pull(max, gene)
    num_mets = df %>% count(index, gene) %>% group_by(gene) %>% summarise(max = max(n, na.rm = FALSE)) %>% pull(max, gene)

    .step("main_p built")
    p_width <- num_res/7
    p_height <- (num_mets  + 3)/4.3

    top_p <- p[["meta_p"]] +
      theme(panel.spacing = unit(0, "pt"),
            plot.margin = margin(t = 10, r = 10, b = 0, l = 80))

    ## plot.margin set per subplot rather than through `&`, so seq_p can carry a
    ## top margin for its index labels while main_p keeps t = 0. Left/right stay
    ## identical on both -- panel_left_in / panel_right_in measure the result, and
    ## the detail panels and strip match whatever it comes out to.
    bottom_p <- (p[["seq_p"]] +
                   theme(plot.margin = margin(t = nn_seq_top_pt, r = 10,
                                              b = 10, l = 80))) /
                (p[["plot"]] +
                   theme(plot.margin = margin(t = 0, r = 10, b = 10, l = 80)))
    bottom_p <- bottom_p + patchwork::plot_layout(heights = c(0.7, p_height)) &
      theme(panel.spacing = unit(0, "pt"))

    n_subtitle_lines <- length(stringr::str_split(subtitle_text, "\n",
                                                  simplify = TRUE))
    top_height_svg <- max(
      0.5 / (0.5 + 0.7 + p_height) * p_height,
      0.5 + n_subtitle_lines * 0.22
    )
    bot_height_svg  <- (0.7 + p_height) / (0.5 + 0.7 + p_height) * p_height

    .step("girafe: top_wgt")
    top_wgt <- ggiraph::girafe(ggobj = top_p,
                               width_svg = p_width,
                               height_svg = top_height_svg,
                               options = list(
                                 ggiraph::opts_sizing(rescale = FALSE),
                                 ggiraph::opts_selection(type = "none")
                               ))

    .step("girafe: bottom (main track)")
    wgt <- ggiraph::girafe(ggobj = bottom_p,
                           width_svg = p_width,
                           height_svg = bot_height_svg,
                           options = list(
                             ggiraph::opts_sizing(rescale = FALSE),
                             ggiraph::opts_selection(type = "none")
                           ))
    wgt$styles <- c(
      wgt$styles,
      "text { user-select: text; pointer-events: all; }"
    )

    wgt$styles <- c(
      wgt$styles,
      "
  text {
    user-select: none;
    pointer-events: all;
  }

  .start-residue {
    fill: orange !important;
    font-weight: bold;
  }

  .selection-rect {
    fill: none;
    stroke-width: 2px;
    rx: 6px;
    ry: 6px;
    pointer-events: none;
  }
  "
    )

    wgt <- htmlwidgets::onRender(
      wgt,
      paste0("
           const GENE   = '", old_nn_input[["gene"]], "';
           const COLORS = ", jsonlite::toJSON(rep(desc_colors, 4)), ";

           window.commentLines     = [];
           window.aliasLines       = [];
           window.aliasNames       = [];
           window.selectionRects   = [];
           window.selectedResidues = [];
           window.colorIndex       = 0;

           function setStatus(msg) {
             const el = document.getElementById('selection-status');
             if (el) el.textContent = msg || '';
           }

           function clearStartHighlight() {
             document.querySelectorAll('text.start-residue')
             .forEach(el => el.classList.remove('start-residue'));
           }

           function highlightStartResidue(idx) {
             clearStartHighlight();
             const el = document.querySelector('text[data-id=\"' + idx + '\"]');
             if (el) el.classList.add('start-residue');
           }

           function getTextEl(idx) {
             return document.querySelector('text[data-id=\"' + idx + '\"]');
           }

           /* Append into the SAME layer as the residue letters. The old
              document.querySelector('svg') returned the FIRST svg on the page --
              the meta/title widget -- so selection rectangles were drawn over the
              title area instead of over seq_p. getBBox() returns coordinates in
              the element's own user space, so sharing t1's parent is what makes
              the rect land on the residues. */
           function getTextLayer(el) {
             if (el && el.parentNode) return el.parentNode;
             const svgs = document.querySelectorAll('svg.ggiraph-svg');
             const svg = svgs.length ? svgs[svgs.length - 1]
                                     : document.querySelector('svg');
             if (!svg) return null;
             const g = svg.querySelector('g');
             return g || svg;
           }

           function drawSelectionRect(start, end, color) {
             const t1 = getTextEl(start);
             const t2 = getTextEl(end);
             if (!t1 || !t2) return null;

             const b1 = t1.getBBox();
             const b2 = t2.getBBox();
             const padY = 1;
             const trimX = 1;
             const x1 = Math.min(b1.x, b2.x);
             const x2 = Math.max(b1.x + b1.width, b2.x + b2.width);
             const top    = Math.min(b1.y, b2.y);
             const bottom = Math.max(b1.y + b1.height, b2.y + b2.height);
             const centerY = (top + bottom) / 2;
             const height  = (bottom - top) + padY * 2;

             const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
             rect.setAttribute('x', x1 + trimX / 2);
             rect.setAttribute('y', centerY - height / 2);
             rect.setAttribute('width', (x2 - x1) - trimX);
             rect.setAttribute('height', height);
             rect.setAttribute('rx', 6);
             rect.setAttribute('ry', 6);
             rect.setAttribute('class', 'selection-rect');
             rect.setAttribute('stroke', color);
             rect.setAttribute('fill', 'none');

             const layer = getTextLayer(t1);
             if (!layer) return null;
             layer.appendChild(rect);
             return rect;
           }

           function refreshBox() {
             const ta = document.getElementById('saved-result');
             if (!ta) return;

             let out = [];
             /* The panel holds the COMPLETE ChimeraX script, not just the
                selections: `close` so re-running never stacks a second copy of
                the model, then the embedded per-gene script (model + aliases +
                setattr values), then the selections below. */
             const base = window.NN_BASE_CXC ||
                          (window.NN_AF_URL ? ('open ' + window.NN_AF_URL) : '');
             if (base) { out.push('close'); out.push(base); out.push(''); }
             out = out.concat(window.commentLines);
             out.push('');
             out = out.concat(window.aliasLines);

             if (window.aliasNames.length) {
               out.push('alias manual_pep ' + window.aliasNames.join('; '));
               out.push('manual_pep');
             }

             /* The open window's ChimeraX colouring, appended rather than
                written straight into the textarea: refreshBox owns this box and
                rebuilds it on every manual selection, so a direct write would be
                wiped by the next one. Recomputed here, so manual selections and
                the window colouring coexist. */
             if (typeof window.nn_current_window_cxc === 'function') {
               const wcxc = window.nn_current_window_cxc();
               if (wcxc) { out.push(''); out.push(wcxc); }
             }

             ta.value = out.join('\\n');
           }

           window.nn_refresh_box = refreshBox;

           window.loadCXC = function(file) {
             const reader = new FileReader();
             reader.onload = function(e) {
               window.commentLines   = [];
               window.aliasLines     = [];
               window.aliasNames     = [];
               window.selectionRects.forEach(r => r?.remove());
               window.selectionRects = [];
               window.colorIndex     = 0;

               e.target.result.split(/\\r?\\n/).forEach(line => {
                 if (line.startsWith('#')) {
                   window.commentLines.push(line);

                   const m = line.match(/_(\\d+)-(\\d+)$/);
                   if (m) {
                     const start = +m[1];
                     const end   = +m[2];
                     const color = COLORS[window.colorIndex % COLORS.length];
                     const rect  = drawSelectionRect(start, end, color);
                     window.selectionRects.push(rect);
                     window.colorIndex++;
                   }
                 } else if (
                   line.startsWith('alias ') &&
                   !line.startsWith('alias manual_pep')
                 ) {
                   window.aliasLines.push(line);
                   const m = line.match(/^alias\\s+(\\S+)/);
                   if (m) window.aliasNames.push(m[1]);
                 }
               });

               refreshBox();
             };
             reader.readAsText(file);
           };

           window.undoLast = function() {
             if (!window.commentLines.length) return;

             window.commentLines.pop();
             window.aliasLines.pop();
             window.aliasNames.pop();
             window.colorIndex = Math.max(0, window.colorIndex - 1);

             const r = window.selectionRects.pop();
             if (r) r.remove();

             window.selectedResidues = [];
             clearStartHighlight();
             setStatus('');
             refreshBox();
           };

           document.addEventListener('click', function(e) {
             const t = e.target;
             if (!t || !t.hasAttribute('data-id')) return;

             const did = t.getAttribute('data-id');
             if (did && did.indexOf('nnwin_') === 0) {
               if (typeof window.nn_show_panel === 'function') window.nn_show_panel(did);
               return;
             }

             const idx = +did;
             if (isNaN(idx)) return;

             window.selectedResidues.push(idx);

             if (window.selectedResidues.length === 1) {
               highlightStartResidue(idx);
               setStatus('Selection started at residue ' + idx);
               return;
             }

             const uniq  = [...new Set(window.selectedResidues)].sort((a,b)=>a-b);
             const start = uniq[0];
             const end   = uniq[uniq.length - 1];

             const alias = 'p' + start + '-' + end;
             if (window.aliasNames.includes(alias)) {
               window.selectedResidues = [];
               clearStartHighlight();
               setStatus('');
               return;
             }

             const color = COLORS[window.colorIndex % COLORS.length];
             window.colorIndex++;

             window.commentLines.push('#' + GENE + '_' + start + '-' + end);

             const deleteLabels =
               window.aliasLines.length === 0 ? 'label delete; ' : '';

             window.aliasLines.push(
               'alias ' + alias + ' ' +
                 'select clear; ' +
                 deleteLabels +
                 'select :' + start + '-' + end + '; ' +
                 'color sel ' + color + '; ' +
                 'select :' + end + '; ' +
                 'label sel text pep_' + start + '-' + end + '; ' +
                 'label height 4; ' +
                 'select clear'
             );

             window.aliasNames.push(alias);

             const rect = drawSelectionRect(start, end, color);
             window.selectionRects.push(rect);

             window.selectedResidues = [];
             clearStartHighlight();
             setStatus('');
             refreshBox();
           });

           window.downloadCXC = function() {
             const ta = document.getElementById('saved-result');
             /* the panel already holds the complete script -- just save it */
             if (!ta?.value) return;
             const blob = new Blob([ta.value + '\\n'], { type: 'text/plain' });
             const a = document.createElement('a');
             a.href = URL.createObjectURL(blob);
             a.download = GENE + '.cxc';
             document.body.appendChild(a);
             a.click();
             document.body.removeChild(a);
           };

           if (!document.getElementById('annotation-box')) {
             const box = document.createElement('div');
             box.id = 'annotation-box';
             box.style.position = 'fixed';
             box.style.top = '10px';
             box.style.left = '10px';
             box.style.width = '260px';
             box.style.zIndex = 9999;
             box.style.fontFamily = 'monospace';

             /* Collapsible. The box is position:fixed over the left edge, which is
                exactly where the y-axis metric labels are -- on a page no wider
                than the viewport there is no horizontal scroll to move them out
                from under it, so the labels (and their hover text) are otherwise
                unreachable. */
             box.innerHTML =
               '<div id=\"cxc-head\" style=\"cursor:pointer;user-select:none;font-weight:bold;margin-bottom:4px;\"' +
               ' onclick=\"nn_toggle_box()\">&#9662; ChimeraX</div>' +
               '<div id=\"cxc-body\">' +
               '<label for=\"cxcFile\">Load existing CXC File:</label>' +
               '<input type=\"file\" id=\"cxcFile\" accept=\".cxc\" onchange=\"loadCXC(this.files[0])\" style=\"width:100%; margin-bottom:4px;\" />' +
               '<div id=\"selection-status\" style=\"color:#b30000;font-weight:bold;margin-bottom:4px;\"></div>' +
               '<textarea id=\"saved-result\" style=\"width:100%; height:240px;\"></textarea>' +
               '<button onclick=\"undoLast()\" style=\"width:100%; margin-top:4px;\">Undo</button>' +
               '<button onclick=\"downloadCXC()\" style=\"width:100%; margin-top:4px;\">Visualize in ChimeraX</button>' +
               '</div>';

             window.nn_toggle_box = function() {
               var body = document.getElementById('cxc-body');
               var head = document.getElementById('cxc-head');
               if (!body) return;
               var hidden = body.style.display === 'none';
               body.style.display = hidden ? 'block' : 'none';
               if (head) head.innerHTML = (hidden ? '&#9662;' : '&#9656;') + ' ChimeraX';
               box.style.width = hidden ? '260px' : 'auto';
               if (typeof nn_place_box === 'function') nn_place_box();
             };

             document.body.appendChild(box);

             /* Keep the CXC box clear of the fixed meta header. The header's
                translateX is computed (see metaExtraOffset / adjust) to line its
                content up with the main track's panel, so nudging the header
                would fight that -- move the box down instead. Purely a fixed
                overlay: touches no SVG, no plot.margin, no transform, so residue
                alignment is unaffected. Re-measured on resize because the header
                height varies with the number of subtitle lines. */
             var nn_place_box = function() {
               var hdr = document.getElementById('top-meta-fixed');
               var h = hdr ? hdr.getBoundingClientRect().height : 0;
               box.style.top = (h + 14) + 'px';
               /* keep the box inside the viewport on short screens */
               /* Indent only as far as actually needed. The plot already carries
                  ~290px of its own left gutter (margin + y-axis labels), so a
                  fixed body indent stacks on top of that and the gap balloons.
                  Zero it, measure where the leftmost axis label really lands,
                  and add back only the shortfall -- usually 0 once the gutter
                  alone clears the panel. setProperty(...,'important') is needed
                  to beat `padding: 0 !important` on html, body. */
               var pad = function(px) {
                 document.body.style.setProperty('padding-left', px + 'px', 'important');
               };
               pad(0);
               var svgs = document.querySelectorAll('svg.ggiraph-svg');
               if (svgs.length) {
                 var main = svgs[svgs.length - 1], minX = Infinity;
                 main.querySelectorAll('text').forEach(function(e) {
                   var r = e.getBoundingClientRect();
                   /* + scrollX: getBoundingClientRect is viewport-relative, and on
                      a deep link the page has already scrolled by the time this
                      runs -- measuring in viewport coords made minX hugely
                      negative and the indent ballooned by the scroll distance.
                      The box is position:fixed, so its viewport x IS its page x
                      at scrollX = 0, which is the frame we want to compare in. */
                   if (r.width > 0) minX = Math.min(minX, r.left + window.scrollX);
                 });
                 if (isFinite(minX)) {
                   var need = box.getBoundingClientRect().right + 12 - minX;
                   pad(Math.max(0, Math.round(need)));
                 }
               }

               /* Only shrink the textarea when the viewport is genuinely short.
                  Guarded because innerHeight can read 0 in headless/preview
                  contexts, which would collapse the box to its 90px floor. */
               var ta = document.getElementById('saved-result');
               if (ta && window.innerHeight > 300) {
                 var avail = window.innerHeight - (h + 14) - 140;
                 ta.style.height = Math.max(90, Math.min(240, avail)) + 'px';
               }
             };
             nn_place_box();
             setTimeout(nn_place_box, 300);
             setTimeout(nn_place_box, 1200);
             window.addEventListener('resize', nn_place_box);
             window.nn_place_box = nn_place_box;

             /* Populate the panel immediately: the embedded script is available
                as soon as the page loads, so the box should never start empty. */
             refreshBox();
           }
           ")
    )

    ## ---- embedded ChimeraX base script ---------------------------------
    ## Generated HERE, from the same data the plot uses, and inlined into the
    ## page: the HTML is then genuinely self-contained -- no .cxc to host, none
    ## read at view time, and the script cannot drift from what is plotted.
    ## make_cm_script_text() lives in inst/scripts/generate_cm_sct.R, so it is
    ## reached through the search path once that file has been sourced; without
    ## it we fall back to opening the AlphaFold model and nothing else.
    nn_base_cxc <- NULL
    if (exists("make_cm_script_text", mode = "function", inherits = TRUE)) {
      nn_base_cxc <- tryCatch(
        paste(make_cm_script_text(old_nn_input, pep_tp_outer,
                                  attr_mode = "setattr"), collapse = "\n"),
        error = function(e) { .step(paste0("chimerax script FAILED: ",
                                           conditionMessage(e))); NULL })
    }
    .step(sprintf("chimerax base script: %s",
                  if (is.null(nn_base_cxc)) "unavailable (source generate_cm_sct.R)"
                  else sprintf("%d lines, %.0f KB",
                               length(strsplit(nn_base_cxc, "\n")[[1]]),
                               nchar(nn_base_cxc) / 1024)))

    nn_base_cxc_js <- htmltools::tags$script(htmltools::HTML(paste0(
      "window.NN_BASE_CXC = ",
      if (is.null(nn_base_cxc)) "null" else jsonlite::toJSON(nn_base_cxc, auto_unbox = TRUE),
      ";\n",
      "window.NN_AF_URL = ",
      jsonlite::toJSON(paste0("https://alphafold.ebi.ac.uk/files/AF-",
                              p_title$accession, "-F1-model_v6.pdb"),
                       auto_unbox = TRUE), ";\n",
      "window.NN_LABEL_TT = ",
      jsonlite::toJSON(as.list(nn_label_tt), auto_unbox = TRUE), ";\n",
      "
      /* Hover text for axis row labels. Delegated off document and matched on
         textContent, so it needs no cooperation from ggiraph and leaves
         axis.text.y as ggtext::element_markdown -- which is what gives the
         species rows their colours. An interactive-guide approach would have to
         take that element over. */
      (function(){
        var tip = null;
        function ensureTip(){
          if (tip) return tip;
          tip = document.createElement('div');
          tip.id = 'nn-axis-tip';
          tip.style.cssText = 'position:fixed;z-index:10000;background:#111;color:#fff;' +
            'font:12px/1.35 sans-serif;padding:5px 8px;border-radius:4px;max-width:320px;' +
            'pointer-events:none;display:none;box-shadow:0 2px 6px rgba(0,0,0,.3)';
          document.body.appendChild(tip);
          return tip;
        }
        /* ggplot/ggiraph emit axis text with pointer-events:none, so a real mouse
           passes straight through and the handler below never fires -- a
           synthetic dispatchEvent does fire, which makes this easy to miss.
           Re-enable hit-testing on ONLY the labels that have descriptions;
           blanket-enabling it on every <text> would put non-interactive labels
           drawn over the tiles in front of their click targets. Re-run after
           render because the widgets populate asynchronously. */
        function armLabels(){
          if (!window.NN_LABEL_TT) return;
          document.querySelectorAll('svg.ggiraph-svg text').forEach(function(el){
            var k = (el.textContent || '').trim();
            if (Object.prototype.hasOwnProperty.call(window.NN_LABEL_TT, k))
              el.style.pointerEvents = 'auto';
          });
        }
        if (document.readyState === 'loading')
          document.addEventListener('DOMContentLoaded', armLabels);
        else armLabels();
        window.addEventListener('load', armLabels);
        setTimeout(armLabels, 400);
        setTimeout(armLabels, 1500);

        document.addEventListener('mouseover', function(e){
          var el = e.target;
          if (!el || el.tagName !== 'text') return;
          var txt = window.NN_LABEL_TT && window.NN_LABEL_TT[(el.textContent || '').trim()];
          if (!txt) return;
          var t = ensureTip();
          t.textContent = txt;
          t.style.display = 'block';
          var r = el.getBoundingClientRect();
          t.style.left = Math.min(window.innerWidth - 330, r.right + 8) + 'px';
          t.style.top  = Math.max(4, r.top - 4) + 'px';
        });
        document.addEventListener('mouseout', function(e){
          if (e.target && e.target.tagName === 'text' && tip) tip.style.display = 'none';
        });
      })();
      "
    )))

    # ----- detail panels (one girafe per NN window) -----
    if (!is.null(nn_anno) && nrow(nn_anno) > 0) {

      max_layer <- max(nn_anno$layer, na.rm = TRUE)

      ## ---- align the NN widgets with the main residue track --------------
      ## bottom_p (seq_p / plot) carries margin(l = 80pt) PLUS the width of the
      ## metric row labels, which varies with the labels themselves -- so the old
      ## fixed 16.6-line margin could not match it. Measure where bottom_p's panel
      ## actually starts, measure where these plots start with no left margin,
      ## and use the difference. Falls back to the fixed margins if the grob
      ## cannot be measured.
      main_left_in <- tryCatch(panel_left_in(bottom_p), error = function(e) NA_real_)

      det_left_in <- strip_left_in <- NULL
      probe_det <- NULL
      if (is.finite(main_left_in)) {
        probe_det <- make_detail_panel(
          per_index = nn_anno$per_index[[1]], meta_data = nn_anno$meta_data[[1]],
          title = "probe", x_range = c(0, max_index + 1), seq_offset = seq_offset,
          win_start = nn_anno$start[1], win_end = nn_anno$end[1],
          left_in = 0, right_in = 0)
        d0 <- tryCatch(panel_left_in(probe_det), error = function(e) NA_real_)
        if (is.finite(d0)) det_left_in <- max(0, main_left_in - d0)
      }

      ## Same treatment for the right edge. Left alone, the panels share a left
      ## edge but not a width, which is a scale mismatch: zero error at the
      ## N-terminus growing toward the C-terminus.
      main_right_in <- tryCatch(panel_right_in(bottom_p), error = function(e) NA_real_)
      det_right_in <- NULL
      if (is.finite(main_right_in) && !is.null(probe_det)) {
        d1 <- tryCatch(panel_right_in(probe_det), error = function(e) NA_real_)
        if (is.finite(d1)) det_right_in <- max(0, main_right_in - d1)
      }

      ## The dibasic site the window was ANCHORED on -- not every dibasic the
      ## window happens to span. It sits at a fixed offset inside the window:
      ## positions 6-7 for an N-terminal window, 30-31 for a C-terminal one.
      ## Derived from the window geometry, so it does not depend on which dibasic
      ## annotation the input happens to carry. Terminus is the N/C suffix of
      ## target_short (db_N, db_C, chym_N, chym_C, pep_end_N, pep_end_C).
      ## Suppress with options(lf.nn_show_db_mark = FALSE).
      db_marks <- NULL
      if (isTRUE(getOption("lf.nn_show_db_mark", TRUE))) {
        term <- stringr::str_extract(as.character(nn_anno$target_short), "[NC]$")
        ## Anchor is at positions 6-7 (N) or 30-31 (C) of the *padded* 36-row
        ## window, but `peps` -- and so start/end -- records only the real residue
        ## range: 9.2_add_contact_data.R pads windows that run off either end of
        ## the protein (ANO8_w1-25 is 25 residues + 11 pad rows at the front).
        ## Measuring each anchor from its own edge is padding-invariant, and for a
        ## full 36-mer still lands on 6-7 and 30-31.
        anch <- ifelse(term == "N", nn_anno$start + 5L,
                ifelse(term == "C", nn_anno$end   - 6L, NA_integer_))
        db_marks <- tibble::tibble(
          layer = nn_anno$layer,
          start = nn_anno$start,
          end   = nn_anno$end,
          anch  = anch,
          ## anchor kind: db, chym, pep_end. target_short is the terminus (N/C);
          ## win_type is the kind -- together they form wt_combo (db_N, chym_C...).
          lab   = as.character(nn_anno$win_type)
        ) %>%
          dplyr::filter(!is.na(anch), anch >= start, anch + 1L <= end)

        .step(sprintf("dibasic anchors: %d N / %d C window(s), %d mark(s)",
                      sum(term == "N", na.rm = TRUE),
                      sum(term == "C", na.rm = TRUE),
                      nrow(db_marks)))
      }

      .step("building nn strip")
      nn_strip_p <- ggplot2::ggplot(nn_anno) +
        ## drawn as rects (not thick segments) so each window can carry a thin
        ## black border; fill is the raw score.
        ggiraph::geom_rect_interactive(
          ggplot2::aes(xmin = start - 0.3, xmax = end + 0.3,
                       ymin = layer - nn_win_half_h, ymax = layer + nn_win_half_h,
                       fill = pred_raw,
                       tooltip = tooltip,
                       data_id = panel_id,
                       onclick = sprintf("nn_show_panel(&quot;%s&quot;)", panel_id)),
          color = "black", linewidth = 0.25) +
        ggplot2::geom_label(
          ggplot2::aes(x = (start + end) / 2,
                       y = layer - 0.45,
                       label = label_txt),
          size = 2.4, color = "black", fill = "white",
          label.size = 0.3, label.r = grid::unit(0.15, "lines"),
          label.padding = grid::unit(0.12, "lines"),
          vjust = 0.5, lineheight = 0.9
        ) +
        ## Same limits as seq_p (:454) and the detail panels (:118). The strip
        ## used c(1, max_index), which with the identical expansion put its data
        ## range one residue to the right of everything else -- every window drew
        ## ~1 residue off. oob_keep stays: a window can overrun the limits.
        ggplot2::scale_x_continuous(limits = c(0, max_index + 1),
                                    expand = ggplot2::expansion(add = c(seq_offset, -0.4)),
                                    oob = scales::oob_keep
        ) +
        ggplot2::scale_y_reverse(limits = c(max_layer + 0.5, -0.1),
                                 breaks = seq_len(max_layer),
                                 labels = nn_axis_lab) +
        ## Segment colour is the raw (uncalibrated) global score. Limits are
        ## pinned to [0,1] -- the score is a sigmoid output -- so the colours mean
        ## the same thing across proteins instead of rescaling per gene; anything
        ## out of range is squished to the end rather than dropped to NA. The
        ## win_type/terminus that used to drive colour now lives in the label.
        ggplot2::scale_fill_viridis_c(
          option = "H", limits = c(0, 1), oob = scales::squish,
          na.value = "grey40", name = "window-level 1D-CNN prediction",
          guide = ggplot2::guide_colourbar(barwidth = grid::unit(6, "lines"),
                                           barheight = grid::unit(0.4, "lines"),
                                           title.position = "left",
                                           title.vjust = 1)) +
        ## Reserve the SAME left gutter as the detail panels. These are separate
        ## girafe widgets stacked in HTML, so nothing aligns them automatically:
        ## the detail panel's y title + tick labels push its panel right, and the
        ## strip (no y axis) started further left. Rendering an invisible y title
        ## and invisible tick labels of identical width -- via the shared
        ## nn_axis_lab() formatter -- makes both panels begin at the same x, so
        ## the residue positions line up and the axis text sits outside it.
        ggplot2::labs(y = "score") +
        ggplot2::theme_void() +
        ggplot2::theme(
          ## same measured offset as the detail panels: the strip reserves an
          ## identical invisible y-axis gutter, so one value aligns both.
          plot.margin = if (is.null(det_left_in))
              unit(nn_det_margs2, "lines")
            else
              ggplot2::margin(t = nn_det_margs2[1], b = nn_det_margs2[3],
                              unit = "lines") +
              ggplot2::margin(r = nn_marg_r_pt, unit = "pt") +
              ggplot2::margin(l = det_left_in, unit = "in") +
              ggplot2::margin(r = if (is.null(det_right_in)) 0 else det_right_in,
                              unit = "in"),
          axis.title.y = ggplot2::element_text(size = 11, angle = 90, color = NA,
                                               margin = ggplot2::margin(r = 2.75)),
          axis.text.y  = ggplot2::element_text(size = 8.8, hjust = 1, color = NA,
                                               margin = ggplot2::margin(r = 2.2)),
          ## theme_bw() reserves 2.75pt for y ticks, theme_void() reserves none;
          ## without this the strip sits 0.038in left of the detail panel.
          axis.ticks.y = ggplot2::element_line(color = NA),
          axis.ticks.length = unit(2.75, "pt"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, "lines"),
          legend.title = ggplot2::element_text(size = 8),
          legend.text = ggplot2::element_text(size = 7)
        )

      ## Drawn after the rects so it sits on top: a white band spanning both
      ## anchor residues, tagged with the anchor kind. White because the rect fill
      ## runs the whole of viridis-C underneath it.
      if (!is.null(db_marks) && nrow(db_marks) > 0) {
        nn_strip_p <- nn_strip_p +
          ggplot2::geom_rect(
            data = db_marks,
            ggplot2::aes(xmin = anch - 0.5, xmax = anch + 1.5,
                         ymin = layer - nn_win_half_h, ymax = layer + nn_win_half_h),
            inherit.aes = FALSE, fill = "white", colour = "black", linewidth = 0.2) +
          ggplot2::geom_label(
            data = db_marks,
            ggplot2::aes(x = anch + 0.5, y = layer, label = lab),
            inherit.aes = FALSE, size = 2.6, colour = "black", fill = "white",
            label.size = 0.15, label.r = grid::unit(0.05, "lines"),
            label.padding = grid::unit(0.10, "lines"))
      }

      .step("girafe: nn strip")
      nn_wgt <- ggiraph::girafe(
        ggobj = nn_strip_p,
        width_svg = p_width,
        height_svg = max(2.5, 0.42 * max_layer + 1.0),   # taller: the label is now 2 lines
        options = list(
          ggiraph::opts_sizing(rescale = FALSE),
          ## ggiraph's default selection css is fill:red;stroke:red, which hides
          ## the pred_raw fill on exactly the window you are looking at. Style the
          ## stroke only -- the fill stays, and a heavier border than the hover
          ## state (1.5) is the sole selection cue.
          ggiraph::opts_selection(type = "single", only_shiny = FALSE,
                                  css = "stroke:black;stroke-width:3;"),
          ggiraph::opts_hover(css = "stroke:black;stroke-width:1.5;cursor:pointer;")
        )
      )

      .step(sprintf("building %d detail panels ...", nrow(nn_anno)))
      detail_widgets <- purrr::map(seq_len(nrow(nn_anno)), function(i) {
        if (i %% 10 == 1) .step(sprintf("  detail panel %d/%d", i, nrow(nn_anno)))
        detail_p <- make_detail_panel(
          per_index = nn_anno$per_index[[i]],
          meta_data = nn_anno$meta_data[[i]],
          # raw (uncalibrated) score, matching the strip label and the rect fill
          ## Non-breaking spaces, not plain ones: svglite emits no xml:space
          ## attribute, so SVG's default whitespace handling collapses a run of
          ## spaces to a single one (measured: "A     B" renders exactly as wide
          ## as "A B"). U+00A0 is not collapsed.
          title = sprintf("peptide window ID: %s\u00a0\u00a0\u00a0\u00a0\u00a0window-level 1D-CNN prediction: %.3f",
                          nn_anno$peps[i], nn_anno$pred_raw[i]),
          x_range = c(0, max_index + 1),
          seq_offset = seq_offset,
          win_start = nn_anno$start[i],
          win_end   = nn_anno$end[i],
          left_in   = det_left_in,
          right_in  = det_right_in
        )
        ggiraph::girafe(
          ggobj = detail_p,
          width_svg = p_width,
          height_svg = 2.5,
          options = list(
            ggiraph::opts_sizing(rescale = FALSE),
            ggiraph::opts_selection(type = "none")
          )
        )
      })

      detail_panel_tags <- purrr::pmap(
        list(detail_widgets, nn_anno$panel_id, nn_anno$start, nn_anno$end,
             nn_anno$peps),
        function(w, id, s, e, pp) {
          htmltools::tags$div(
            ## data-peps is the deep-link key: panel_id is positional (windows are
            ## ordered by desc(pred)) so it shifts whenever scores change, while
            ## peps -- GENE_w<start>-<end> -- is stable and readable in a URL.
            class = "nn-detail-panel", id = id,
            `data-peps`      = as.character(pp),
            `data-win-start` = as.integer(s),
            `data-win-end`   = as.integer(e),
            w
          )
        }
      )

      ## ---- ChimeraX colouring for the open window -----------------------
      ## Per residue, which class has the highest prediction. Emitted as data so
      ## the browser can build a .cxc on demand -- the alternative, scraping the
      ## rendered tooltips, breaks the moment the tooltip format changes.
      nn_cxc_dat <- purrr::map(seq_len(nrow(nn_anno)), function(i) {
        pi <- nn_anno$per_index[[i]]
        md <- nn_anno$meta_data[[i]]
        cls <- setdiff(names(pi), "index")
        m <- as.matrix(pi[, cls, drop = FALSE])
        storage.mode(m) <- "double"
        wm <- max.col(replace(m, is.na(m), -Inf), ties.method = "first")
        ## padded rows have no residue index, and a row that is entirely NA has
        ## no winner -- drop both rather than colouring them arbitrarily
        keep <- !is.na(md$index) & rowSums(!is.na(m)) > 0
        list(peps  = nn_anno$peps[i],
             start = as.integer(nn_anno$start[i]),
             end   = as.integer(nn_anno$end[i]),
             index = as.integer(md$index[keep]),
             class = cls[wm[keep]],
             value = round(unname(m[cbind(seq_len(nrow(m)), wm)])[keep], 4))
      })
      names(nn_cxc_dat) <- nn_anno$panel_id

      nn_cxc_js <- htmltools::tags$script(htmltools::HTML(paste0(
        "window.NN_CXC = ", jsonlite::toJSON(nn_cxc_dat, auto_unbox = TRUE), ";\n",
        "window.NN_CLASS_COLS = ",
        jsonlite::toJSON(as.list(nn_class_cols), auto_unbox = TRUE), ";\n",
        "window.NN_GENE = ", jsonlite::toJSON(as.character(gene_tp), auto_unbox = TRUE), ";\n",
        "
        /* Build a ChimeraX command file colouring each residue of one window by
           its highest-scoring class. Same palette as the detail panel, and the
           same alias idiom generate_cm_sct.R writes: reset to the base colour,
           then one select/color pair per run of consecutive residues sharing a
           class. One fixed alias, so each window replaces the last. */
        window.nn_window_cxc = function(panelId) {
          var d = window.NN_CXC && window.NN_CXC[panelId];
          if (!d || !d.index || !d.index.length) return null;
          var cols = window.NN_CLASS_COLS || {};
          var idx = [].concat(d.index), cl = [].concat(d.class);
          var parts = [], i = 0;
          while (i < idx.length) {
            var j = i;
            while (j + 1 < idx.length && cl[j + 1] === cl[i] &&
                   idx[j + 1] === idx[j] + 1) j++;
            var col = cols[cl[i]];
            if (col) {
              var rng = (idx[i] === idx[j]) ? (':' + idx[i])
                                            : (':' + idx[i] + '-' + idx[j]);
              parts.push('select ' + rng + '; color sel ' + col);
            }
            i = j + 1;
          }
          if (!parts.length) return null;
          /* Fixed alias name, deliberately not per-peptide: loading a second
             window then REPLACES the first in ChimeraX instead of leaving a
             pile of nn_GENE_wA-B aliases behind. Combined with the reset at the
             head of the body, running it clears the previous window's colouring
             and label before applying this one. Trade-off: you cannot keep two
             windows aliased at once and toggle between them. */
          var alias = 'nn_window';
          /* No `select all; color #D2AF81FF` reset here any more. The composed
             file starts with `close`, so nothing from a previous window can
             survive to be reset -- and the reset actively wiped the base
             script's own colouring, which now runs immediately above this. */
          var body = ['label delete']
                       .concat(parts)
                       .concat(['select :' + d.end,
                                'label sel text ' + String(d.peps).replace(/[^A-Za-z0-9_-]/g, '_'),
                                'label height 4', 'select clear']).join('; ');
          return ['# ' + d.peps + '  max-prediction colouring  (' +
                    idx.length + ' residues, ' + parts.length + ' runs)',
                  'alias ' + alias + ' ' + body,
                  alias].join('\\n');
        };

        /* Only while a panel is actually shown -- the .active class survives
           nn_hide_panel, so gate on .shown or the block would linger after close. */
        window.nn_current_window_cxc = function() {
          var c = document.getElementById('nn-detail-container');
          if (!c || !c.classList.contains('shown')) return null;
          var a = c.querySelector('.nn-detail-panel.active');
          return a ? window.nn_window_cxc(a.id) : null;
        };

        window.nn_download_cxc = function() {
          var c = document.getElementById('nn-detail-container');
          var a = c && c.querySelector('.nn-detail-panel.active');
          if (!a) return;
          var txt = window.nn_current_window_cxc();
          if (!txt) return;
          var blob = new Blob([txt + '\\n'], { type: 'text/plain' });
          var el = document.createElement('a');
          el.href = URL.createObjectURL(blob);
          el.download = (a.getAttribute('data-peps') || a.id) + '.cxc';
          document.body.appendChild(el); el.click(); document.body.removeChild(el);
          setTimeout(function(){ URL.revokeObjectURL(el.href); }, 1000);
        };
        "
      )))

      detail_container <- htmltools::tags$div(
        id = "nn-detail-container",
        `data-max-res` = as.integer(max_index),
        ## Absolutely positioned, NOT float:right. girafe_container_std centres
        ## the SVG in its container, and a float narrows the line box it centres
        ## within -- so a floated button shifted every detail SVG left by half its
        ## own width (~23.6px, ~2.3 residues) relative to the strip and main
        ## track, but only once the window was wider than the SVG and there was
        ## slack to centre in. Out of flow, all three centre identically.
        htmltools::tags$button(
          "close",
          onclick = "nn_hide_panel()",
          style = "position:absolute; top:4px; right:4px; z-index:2;"
        ),
        ## Also absolutely positioned -- see the note above; anything in normal
        ## flow here shifts the centred detail SVG out of alignment.
        htmltools::tags$button(
          "CXC",
          onclick = "nn_download_cxc()",
          title = "ChimeraX colouring for the window currently open",
          style = "position:absolute; top:4px; right:52px; z-index:2;"
        ),
        detail_panel_tags
      )

      nn_panel_css <- htmltools::tags$style(htmltools::HTML(
        "html, body {
         margin: 0 !important;
         padding: 0 !important;
         overflow: visible !important;
         min-height: 100vh;
       }
       /* Clear the fixed CXC panel (left:10, 260 wide). Indenting the BODY
          shifts every widget -- strip, detail panels, main track -- by the same
          amount, so their relative alignment is untouched; the top meta bar is
          position:fixed and recomputes its own offset from the main panel.
          Later rule of equal specificity, so it beats the padding:0 above. */
       /* Fallback only -- nn_place_box() measures the real shortfall and
          overrides this. ~70px is what the measurement lands on in practice:
          the y-axis labels start around x=216 and the panel ends at 270. */
       body { padding-left: 70px !important; }
       #top-meta-fixed {
         position: fixed !important;
         top: 0 !important;
         left: 0 !important;
         z-index: 1000 !important;
         background: white !important;
         border-bottom: 1px solid #ddd;
         padding-bottom: 2px;
         transform: translateX(260px);
         box-shadow: 0 2px 4px rgba(0,0,0,0.05);
         width: max-content;
       }
       #nn-detail-container {
         display: block;
         position: relative;   /* anchor for the absolutely positioned close button */
         width: 100%;
         /* horizontal padding would shift the detail SVG right of the main
            track's -- a flat ~0.36 residue offset at p_width = num_res/7 */
         padding: 4px 0;
         background: white;
         border-top: 1px solid #ccc;
         margin-top: 6px;
       }
       /* After widgets init, hide container until a panel is shown */
       #nn-detail-container.initialized:not(.shown) {
         display: none;
       }
       /* When shown, only display the active panel */
       #nn-detail-container.shown .nn-detail-panel:not(.active) {
         display: none;
       }
       .nn-detail-panel {
         display: block;
         background: white;
         min-height: 240px;
       }
       /* Cross-widget residue hover highlight */
       .nn-res-hover {
         stroke: #FF6600 !important;
         stroke-width: 1.5 !important;
         fill: #FF6600 !important;
       }"
      ))

      nn_panel_js <- htmltools::tags$script(htmltools::HTML(
        "(function(){
         function init(){
           nn_status('JS loaded ' + new Date().toLocaleTimeString());
         }
         function postInit(){
           var c = document.getElementById('nn-detail-container');
           if (!c) return;
           var hw = c.querySelectorAll('div.html-widget');
           var scripts = c.querySelectorAll('script[type=\"application/json\"]');
           var svgs0 = c.querySelectorAll('svg.ggiraph-svg').length;
           var firstInner0 = hw.length > 0 ? hw[0].innerHTML.length : -1;
           nn_status('pre: hw=' + hw.length + ' scripts=' + scripts.length +
                     ' svgs=' + svgs0 + ' inner0=' + firstInner0);
           if (window.HTMLWidgets && window.HTMLWidgets.staticRender) {
             try { window.HTMLWidgets.staticRender(); }
             catch(e) { nn_status('staticRender ERR: ' + e.message); }
           } else {
             nn_status('no HTMLWidgets.staticRender');
           }
           setTimeout(function(){
             var svgs1 = c.querySelectorAll('svg.ggiraph-svg').length;
             var firstInner1 = hw.length > 0 ? hw[0].innerHTML.length : -1;
             nn_status('post: svgs=' + svgs1 + ' inner0=' + firstInner1);
             c.classList.add('initialized');
             /* Deep link. The strip may still be rendering, so retry briefly. */
             if (location.hash) {
               var tries = 0;
               (function tryJump(){
                 if (window.nn_jump_to && window.nn_jump_to(location.hash, false)) return;
                 if (++tries < 8) setTimeout(tryJump, 250);
               })();
             }
           }, 50);
         }
         if (document.readyState === 'loading') {
           document.addEventListener('DOMContentLoaded', init);
         } else {
           init();
         }
         if (document.readyState === 'complete') {
           setTimeout(postInit, 50);
         } else {
           window.addEventListener('load', function(){ setTimeout(postInit, 50); });
         }
         window.addEventListener('hashchange', function(){
           if (window.nn_jump_to) window.nn_jump_to(location.hash, true);
         });
       })();
       function nn_status(msg) {}

       /* ---- deep linking: page.html#ANO8_w37-72 -------------------------
          Opens that window's detail panel and centres it horizontally. The
          fragment matches, in order: data-peps (GENE_w37-72), the panel id
          (nnwin_3), the bare window (w37-72) or its range (37-72). Case
          insensitive. Returns false if nothing matched, so the caller can
          retry while the widgets are still rendering.                     */
       function nn_center_window(id, smooth) {
         var svgs = document.querySelectorAll('svg.ggiraph-svg');
         var sel = '[data-id=' + JSON.stringify(id) + ']';
         var el = null, i;
         for (i = 0; i < svgs.length; i++) {
           el = svgs[i].querySelector(sel);
           if (el) break;
         }
         var cx = null, r;
         if (el) {
           r = el.getBoundingClientRect();
           cx = r.left + r.width / 2 + window.scrollX;
         } else {
           /* strip rect not found -- fall back to the main track's residues */
           var pnl = document.getElementById(id);
           if (!pnl || !svgs.length) return;
           var main = svgs[svgs.length - 1];
           var a = main.querySelector('[data-id=' + JSON.stringify(pnl.getAttribute('data-win-start')) + ']');
           var b = main.querySelector('[data-id=' + JSON.stringify(pnl.getAttribute('data-win-end')) + ']');
           if (!a || !b) return;
           var ra = a.getBoundingClientRect(), rb = b.getBoundingClientRect();
           cx = (ra.left + ra.width / 2 + rb.left + rb.width / 2) / 2 + window.scrollX;
         }
         var left = Math.max(0, cx - window.innerWidth / 2);
         try { window.scrollTo({ left: left, behavior: smooth ? 'smooth' : 'auto' }); }
         catch (err) { window.scrollTo(left, window.scrollY); }
       }

       window.nn_jump_to = function(token, smooth) {
         var c = document.getElementById('nn-detail-container');
         if (!c || !token) return false;
         var want = String(token).replace(/^#/, '');
         try { want = decodeURIComponent(want); } catch (e) {}
         want = want.trim();
         if (!want) return false;
         var lc = want.toLowerCase();
         var panels = c.querySelectorAll('.nn-detail-panel'), hit = null;
         for (var i = 0; i < panels.length; i++) {
           var pn = panels[i];
           var peps = (pn.getAttribute('data-peps') || '').toLowerCase();
           var s = pn.getAttribute('data-win-start'), e = pn.getAttribute('data-win-end');
           if (pn.id.toLowerCase() === lc || peps === lc ||
               ('w' + s + '-' + e) === lc || (s + '-' + e) === lc ||
               peps.replace(/^.*_/, '') === lc) { hit = pn; break; }
         }
         if (!hit) { nn_status('jump: no match for ' + want); return false; }
         window.nn_show_panel(hit.id);
         /* Opening the panel changes layout, and the body indent and the fixed
            top bar are both computed from measurements -- settle those FIRST, or
            we scroll to a position that is stale by the time they run. */
         if (typeof window.nn_place_box === 'function') {
           try { window.nn_place_box(); } catch (err) {}
         }
         try { window.dispatchEvent(new Event('resize')); } catch (err) {}
         setTimeout(function(){ nn_center_window(hit.id, smooth !== false); }, 60);
         return true;
       };
       window.nn_show_panel = function(id) {
         nn_status('SHOW called id=' + id);
         try {
           var c = document.getElementById('nn-detail-container');
           if (!c) { nn_status('  no nn-detail-container'); return; }
           c.classList.add('shown');
           var panels = c.querySelectorAll('.nn-detail-panel');
           var matched = 0, activePanel = null;
           for (var i = 0; i < panels.length; i++) {
             var isMatch = panels[i].id === id;
             panels[i].classList.toggle('active', isMatch);
             if (isMatch) { matched++; activePanel = panels[i]; }
           }
           /* keep the URL in sync so the address bar is always a shareable
              deep link to whatever window is currently open */
           if (activePanel) {
             var pp = activePanel.getAttribute('data-peps');
             if (pp && window.history && history.replaceState) {
               try { history.replaceState(null, '', '#' + encodeURIComponent(pp)); }
               catch (err) {}
             }
           }
           if (window.nn_refresh_box) { try { window.nn_refresh_box(); } catch (err) {} }
           var rect = c.getBoundingClientRect();
           nn_status('  matched ' + matched + '/' + panels.length +
                     ' box ' + Math.round(rect.width) + 'x' + Math.round(rect.height));
           // Panel alignment is handled entirely in R: det_left_in / det_right_in
           // measure bottom_p's panel edges and match them via plot.margin.
         } catch (err) {
           nn_status('  ERROR ' + err.message);
         }
       };
       window.nn_hide_panel = function() {
         var c = document.getElementById('nn-detail-container');
         if (c) c.classList.remove('shown');
         if (window.nn_refresh_box) { try { window.nn_refresh_box(); } catch (err) {} }
       };
       // Cross-widget residue hover: highlight all elements sharing a numeric data-id
       (function(){
         var residueIdRe = /^\\d+$/;
         function findResidueDataId(el) {
           var n = el, hops = 0;
           while (n && n.getAttribute && hops < 4) {
             var d = n.getAttribute('data-id');
             if (d && residueIdRe.test(d)) return d;
             n = n.parentNode;
             hops++;
           }
           return null;
         }
         function setHover(did, on) {
           var els = document.querySelectorAll('[data-id=\"' + did + '\"]');
           for (var i = 0; i < els.length; i++) {
             els[i].classList.toggle('nn-res-hover', on);
           }
         }
         document.addEventListener('mouseover', function(e) {
           var did = findResidueDataId(e.target);
           if (did) setHover(did, true);
         });
         document.addEventListener('mouseout', function(e) {
           var did = findResidueDataId(e.target);
           if (did) setHover(did, false);
         });
       })();
       document.addEventListener('click', function(e) {
         var did_seen = null;
         var via = '';
         var n = e.target;
         var hops = 0;
         while (n && n.getAttribute && hops < 6) {
           var did = n.getAttribute('data-id');
           if (did && did.indexOf('nnwin_') === 0) { did_seen = did; via = 'parent'; break; }
           n = n.parentNode;
           hops++;
         }
         if (!did_seen && document.elementsFromPoint) {
           var stack = document.elementsFromPoint(e.clientX, e.clientY);
           for (var i = 0; i < stack.length; i++) {
             var d = stack[i].getAttribute && stack[i].getAttribute('data-id');
             if (d && d.indexOf('nnwin_') === 0) { did_seen = d; via = 'point'; break; }
           }
         }
         nn_status('click <' + (e.target && e.target.tagName) + '> via=' + (via || 'none') + ' did=' + (did_seen || '(none)'));
         if (did_seen) window.nn_show_panel(did_seen);
       }, true);"
      ))

      top_wgt_fixed <- htmltools::tags$div(
        id = "top-meta-fixed",
        style = paste0(
          "position: fixed; top: 0; left: 0; z-index: 100; ",
          "background: white; border-bottom: 1px solid #ddd; ",
          "padding-bottom: 2px; width: max-content; ",
          "box-shadow: 0 2px 4px rgba(0,0,0,0.05);"
        ),
        top_wgt
      )

      top_meta_spacer <- htmltools::tags$div(
        id = "top-meta-spacer",
        style = "height: 60px;"
      )

      top_meta_pin_js <- htmltools::tags$script(htmltools::HTML(
        "(function(){
         var horizOffset = 0;
         var metaExtraOffset = 190;
         function svgScale(svg) {
           var rect = svg.getBoundingClientRect();
           var vb = svg.viewBox && svg.viewBox.baseVal;
           return (vb && vb.width > 0) ? rect.width / vb.width : 1;
         }
         function findPanelRect(svg) {
           // ggplot SVGs have multiple clipPath rects: one paper-level (x≈0,
           // full width) and one per panel (x>0 because of left margin + y-axis).
           // We want the panel: largest area among rects with x > 0.
           var rects = svg.querySelectorAll('defs clipPath rect');
           var best = null, bestArea = 0;
           for (var k = 0; k < rects.length; k++) {
             var x = parseFloat(rects[k].getAttribute('x'));
             var w = parseFloat(rects[k].getAttribute('width'));
             var h = parseFloat(rects[k].getAttribute('height'));
             if (!(w > 0 && h > 0)) continue;
             if (!(x > 0)) continue;  // skip paper-level rect at x=0
             var area = w * h;
             if (area > bestArea) { bestArea = area; best = {x:x, w:w, h:h}; }
           }
           return best;
         }
         function adjust() {
           var f = document.getElementById('top-meta-fixed');
           if (!f) return;
           var s = document.getElementById('top-meta-spacer');
           if (s) {
             var h = f.getBoundingClientRect().height;
             if (h > 0) s.style.height = (h + 4) + 'px';
           }
           // Bottom protein plot = LAST html-widget on page, excluding the
           // fixed bar's own widget and any inside the detail container.
           var all = document.querySelectorAll('div.html-widget');
           var detail = document.getElementById('nn-detail-container');
           var bottom = null;
           for (var i = all.length - 1; i >= 0; i--) {
             if (f.contains(all[i])) continue;
             if (detail && detail.contains(all[i])) continue;
             bottom = all[i]; break;
           }
           if (!bottom) return;
           var bSvg = bottom.querySelector('svg.ggiraph-svg');
           if (!bSvg) return;
           var bp = findPanelRect(bSvg);
           if (!bp) return;
           var bScale = svgScale(bSvg);
           var bRect = bSvg.getBoundingClientRect();
           // Convert panel left to PAGE coords so measurement works whether
           // the user has scrolled or not. At scroll=0 this equals the panel's
           // viewport-x — and that is the constant viewport-x we want the bar
           // pinned at.
           var panelPageLeft = bRect.left + window.pageXOffset + bp.x * bScale;

           var mSvg = f.querySelector('svg.ggiraph-svg');
           if (!mSvg) return;
           var mScale = svgScale(mSvg);
           var mRect = mSvg.getBoundingClientRect();
           // The bar is position:fixed, so mRect.left is the meta SVG's
           // viewport-x; the SVG's offset within the bar is constant:
           //   svgOffsetInBar = mRect.left - currentHorizOffset
           // After transform translateX(H), content viewport-x =
           //   H + svgOffsetInBar + 80 * mScale
           // We want this to equal panelPageLeft (held constant across scroll).
           var svgOffsetInBar = mRect.left - horizOffset;
           horizOffset = panelPageLeft - svgOffsetInBar - 80 * mScale + metaExtraOffset;
           f.style.transform = 'translateX(' + horizOffset + 'px)';
         }
         if (document.readyState === 'complete') setTimeout(adjust, 100);
         else window.addEventListener('load', function(){ setTimeout(adjust, 100); });
         setTimeout(adjust, 600);
         setTimeout(adjust, 1500);
         window.addEventListener('resize', adjust);
       })();"
      ))

      .step("detail panels done; assembling page")
      page <- htmltools::tagList(
        nn_panel_css,
        nn_base_cxc_js,
        nn_cxc_js,
        nn_panel_js,
        top_meta_pin_js,
        top_wgt_fixed,
        top_meta_spacer,
        nn_wgt,
        detail_container,
        wgt
      )

      out_file <- fs::path(plot_dir, old_nn_input[["gene"]], ext = "html")
      .step("writing html")
      htmltools::save_html(
        page,
        file = out_file,
        libdir = fs::path_file(fs::path(plot_dir, "dependency_files"))
      )
      ## options(lf.selfcontained = FALSE) to keep the external references
      if (isTRUE(getOption("lf.selfcontained", TRUE))) {
        n <- nn_inline_deps(out_file, fs::path(plot_dir, "dependency_files"))
        .step(sprintf("inlined %s dependencies -> %.1f MB",
                      n, file.size(out_file) / 2^20))
      }
    } else {
      page <- htmltools::tagList(nn_base_cxc_js, top_wgt, wgt)
      out_file <- fs::path(plot_dir, old_nn_input[["gene"]], ext = "html")
      .step("writing html")
      htmltools::save_html(
        page,
        file = out_file,
        libdir = fs::path_file(fs::path(plot_dir, "dependency_files"))
      )
      ## options(lf.selfcontained = FALSE) to keep the external references
      if (isTRUE(getOption("lf.selfcontained", TRUE))) {
        n <- nn_inline_deps(out_file, fs::path(plot_dir, "dependency_files"))
        .step(sprintf("inlined %s dependencies -> %.1f MB",
                      n, file.size(out_file) / 2^20))
      }
    }

  }, error = function(e) {
    message("make_protein_plot_win failed for gene ",
            old_nn_input[["gene"]], ": ", conditionMessage(e))
    message(paste(utils::limitedLabels(sys.calls()), collapse = "\n"))
  })

}
