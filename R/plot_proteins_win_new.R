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

nn_det_margs <- c(0, 1.6,2,16.6)
nn_det_margs2 <- c(0, 1.6,2,18.6)


assign_overlap_layers <- function(start, end) {
  ord <- order(-(end - start))
  layer_end <- numeric(0)
  out <- integer(length(start))
  for (i in ord) {
    placed <- FALSE
    for (k in seq_along(layer_end)) {
      if (start[i] > layer_end[k]) {
        out[i] <- k
        layer_end[k] <- end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      layer_end <- c(layer_end, end[i])
      out[i] <- length(layer_end)
    }
  }
  out
}

make_detail_panel <- function(per_index, meta_data, title = NULL,
                              x_range = NULL, seq_offset = -0.45,
                              win_start = NULL, win_end = NULL) {
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
  if (!is.null(win_start) && !is.null(win_end) &&
      diff(x_range) > 0) {
    win_center <- (win_start + win_end) / 2
    rel <- (win_center - x_range[1]) / diff(x_range)
    rel <- max(0.02, min(0.98, rel))
    legend_pos <- c(rel, 1.02)
  }

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
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::labs(title = title, x = NULL, y = "score") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = legend_pos,
                   legend.justification = "center",
                   legend.direction = "horizontal",
                   legend.background = ggplot2::element_rect(fill = "white", color = NA),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.margin = unit(nn_det_margs, "lines"))
  p
}

make_protein_plot_win <- function(old_nn_input,
                                  new_nn_input,
                                  plot_dir,
                                  pep_tp) {

  tryCatch({

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
          nn_w <- nn_w %>%
            dplyr::mutate(
              target_short = stringr::str_replace(as.character(target), "^loop_", ""),
              wt_combo     = paste0(win_type, "_", target_short),
              panel_id     = sprintf("nnwin_%d", dplyr::row_number()),
              feature      = wt_combo,
              # label above each bar: category rank | win type_terminus | raw score
              #                       nearest known peptide window
              label_txt    = sprintf("#%s  %s  raw %.2f\nnn: %s",
                                     as.character(rank_cat), wt_combo, pred_raw,
                                     dplyr::coalesce(as.character(nn_closest_peptide), "NA")),
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

    # ----- begin near-verbatim copy of make_protein_plot body -----
    to_plot <- expand_by_residue(old_nn_input)

    to_plot <- to_plot %>%
      group_by(gene) %>%
      mutate(max = max(index)) %>%
      mutate(across(starts_with(c("pep_nn", "pep_xgb", "chem_nn", "chem_xgb")),
                    .fns = ~smoother_func(x = ., append_name = "s"),
                    .unpack = TRUE))

    transfrom_aligments <- function(x, vals_to = "AA") {
      x %>%
        mutate(index = row_number()) %>%
        pivot_longer(cols = -index, names_to = "metric", values_to = vals_to)
    }

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
    }, error = function(e) return(NULL))

    to_plot_c <- to_plot %>%
      pivot_longer(cols = any_of(all_mets %>% unname),
                   names_to = "metric",
                   values_to = "value")

    desc_vars <- c(DBC,
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

    met_order <- c(all_mets, desc_vars, levels(species_dat[["aminode"]]))

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

    all_links <- c("UniProt",
                   "Aminode",
                   "chatGPT",
                   "GWAS",
                   "PubMed",
                   "Disease",
                   "AlphaFoldDB",
                   "ChimeraX")

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
      mutate(index2 = seq.default(0, 70, by = 10)) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(tt_value = if(value %in% all_links) {p_title[[value]]} else {NA})

    subtitle_text <- bind_cols(tibble(peptides = list(pep_tp$roi_name)),
                               tibble(selection = list(pep_tp$selection)),
                               anno1)

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
                           position = position_nudge(y = 0.25))

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

      if(y %in% c("NTC", "CTC")) {

        main_p <- main_p +
          ggiraph::geom_point_interactive(data = dat_toplot %>%
                                            filter(!is.na(value_desc)) %>%
                                            mutate(tt_value = value_desc),
                                          aes(x = index, y = metric, tooltip = tt_value, data_id = index), pch = 15, size = 2.5, color = "grey85", show.legend = FALSE) +
          ggplot2::geom_text(data = dat_toplot %>%
                               filter(!is.na(value_desc)),
                             aes(x = index, y = metric, label = AA), size = 1.8, fontface = "bold", color = "black")

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
            legend.key.height = unit(0.2, "cm"))

    p <- list(plot = main_p,
              meta_p = meta_p,
              seq_p = seq_p,
              title_p = title_text,
              subtitle_p = subtitle_text,
              links = p_title)

    num_res = df %>% count(metric, gene) %>% group_by(gene) %>% summarise(max = max(n, na.rm = FALSE)) %>% pull(max, gene)
    num_mets = df %>% count(index, gene) %>% group_by(gene) %>% summarise(max = max(n, na.rm = FALSE)) %>% pull(max, gene)

    p_width <- num_res/7
    p_height <- (num_mets  + 3)/4.3

    top_p <- p[["meta_p"]] +
      theme(panel.spacing = unit(0, "pt"),
            plot.margin = margin(t = 10, r = 10, b = 0, l = 80))

    bottom_p <- p[["seq_p"]] / p[["plot"]]
    bottom_p <- bottom_p + patchwork::plot_layout(heights = c(0.7, p_height)) &
      theme(panel.spacing = unit(0, "pt"),
            plot.margin = margin(t = 0, r = 10, b = 10, l = 80))

    n_subtitle_lines <- length(stringr::str_split(subtitle_text, "\n",
                                                  simplify = TRUE))
    top_height_svg <- max(
      0.5 / (0.5 + 0.7 + p_height) * p_height,
      0.5 + n_subtitle_lines * 0.22
    )
    bot_height_svg  <- (0.7 + p_height) / (0.5 + 0.7 + p_height) * p_height

    top_wgt <- ggiraph::girafe(ggobj = top_p,
                               width_svg = p_width,
                               height_svg = top_height_svg,
                               options = list(
                                 ggiraph::opts_sizing(rescale = FALSE),
                                 ggiraph::opts_selection(type = "none")
                               ))

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

           function getTextLayer() {
             const svg = document.querySelector('svg');
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

             getTextLayer().appendChild(rect);
             return rect;
           }

           function refreshBox() {
             const ta = document.getElementById('saved-result');
             if (!ta) return;

             let out = [];
             out = out.concat(window.commentLines);
             out.push('');
             out = out.concat(window.aliasLines);

             if (window.aliasNames.length) {
               out.push('alias manual_pep ' + window.aliasNames.join('; '));
               out.push('manual_pep');
             }

             ta.value = out.join('\\n');
           }

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

             box.innerHTML =
               '<label for=\"cxcFile\">Load existing CXC File:</label>' +
               '<input type=\"file\" id=\"cxcFile\" accept=\".cxc\" onchange=\"loadCXC(this.files[0])\" style=\"width:100%; margin-bottom:4px;\" />' +
               '<div id=\"selection-status\" style=\"color:#b30000;font-weight:bold;margin-bottom:4px;\"></div>' +
               '<textarea id=\"saved-result\" style=\"width:100%; height:240px;\"></textarea>' +
               '<button onclick=\"undoLast()\" style=\"width:100%; margin-top:4px;\">Undo</button>' +
               '<button onclick=\"downloadCXC()\" style=\"width:100%; margin-top:4px;\">View or Save Selections</button>';

             document.body.appendChild(box);
           }
           ")
    )

    # ----- detail panels (one girafe per NN window) -----
    if (!is.null(nn_anno) && nrow(nn_anno) > 0) {

      max_layer <- max(nn_anno$layer, na.rm = TRUE)

      nn_strip_p <- ggplot2::ggplot(nn_anno) +
        ggiraph::geom_segment_interactive(
          ggplot2::aes(x = start - 0.3, xend = end + 0.3,
                       y = layer, yend = layer,
                       color = wt_combo,
                       tooltip = tooltip,
                       data_id = panel_id,
                       onclick = sprintf("nn_show_panel(&quot;%s&quot;)", panel_id)),
          linewidth = 3, lineend = "round") +
        ggplot2::geom_label(
          ggplot2::aes(x = (start + end) / 2,
                       y = layer - 0.45,
                       label = label_txt),
          size = 2.4, color = "black", fill = "white",
          label.size = 0.3, label.r = grid::unit(0.15, "lines"),
          label.padding = grid::unit(0.12, "lines"),
          vjust = 0.5, lineheight = 0.9
        ) +
        ggplot2::scale_x_continuous(limits = c(1, max_index),
                                    expand = ggplot2::expansion(add = c(seq_offset, -0.4)),
                                    oob = scales::oob_keep
        ) +
        ggplot2::scale_y_reverse(limits = c(max_layer + 0.5, -0.1)) +
        ggplot2::scale_color_manual(values = nn_win_target_cols,
                                    na.value = "grey40",
                                    name = NULL, drop = FALSE, guide = "none") +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.margin = unit(nn_det_margs2, "lines"),
          legend.position = "none",
          legend.key.size = unit(0.4, "lines"),
          legend.text = ggplot2::element_text(size = 8)
        )

      nn_wgt <- ggiraph::girafe(
        ggobj = nn_strip_p,
        width_svg = p_width,
        height_svg = max(2.5, 0.42 * max_layer + 1.0),   # taller: the label is now 2 lines
        options = list(
          ggiraph::opts_sizing(rescale = FALSE),
          ggiraph::opts_selection(type = "single", only_shiny = FALSE),
          ggiraph::opts_hover(css = "stroke-width:5;cursor:pointer;")
        )
      )

      detail_widgets <- purrr::map(seq_len(nrow(nn_anno)), function(i) {
        detail_p <- make_detail_panel(
          per_index = nn_anno$per_index[[i]],
          meta_data = nn_anno$meta_data[[i]],
          title = sprintf("%s  (score %.3f)",
                          nn_anno$peps[i], nn_anno$pred[i]),
          x_range = c(0, max_index + 1),
          seq_offset = seq_offset,
          win_start = nn_anno$start[i],
          win_end   = nn_anno$end[i]
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
        list(detail_widgets, nn_anno$panel_id, nn_anno$start, nn_anno$end),
        function(w, id, s, e) {
          htmltools::tags$div(
            class = "nn-detail-panel", id = id,
            `data-win-start` = as.integer(s),
            `data-win-end`   = as.integer(e),
            w
          )
        }
      )

      detail_container <- htmltools::tags$div(
        id = "nn-detail-container",
        `data-max-res` = as.integer(max_index),
        htmltools::tags$button(
          "close",
          onclick = "nn_hide_panel()",
          style = "float:right;"
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
         width: 100%;
         padding: 4px;
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
       })();
       function nn_status(msg) {}
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
           var rect = c.getBoundingClientRect();
           nn_status('  matched ' + matched + '/' + panels.length +
                     ' box ' + Math.round(rect.width) + 'x' + Math.round(rect.height));
           // Align detail panel's left edge with main_p panel via scaleX anchored at right edge.
           // Right edges naturally align (same SVG width, same right margin), so scaling down
           // toward the right shrinks the wider detail panel until its left edge meets main_p's.
           if (activePanel) {
             var detailSvg = activePanel.querySelector('svg.ggiraph-svg');
             var mainEl = c.nextElementSibling;
             var mainSvg = mainEl ? mainEl.querySelector('svg.ggiraph-svg') : null;
             if (detailSvg && mainSvg) {
               var dRects = detailSvg.querySelectorAll('defs clipPath rect');
               var dBest = 0, dlm = NaN, ddw = NaN;
               for (var k = 0; k < dRects.length; k++) {
                 var dw = parseFloat(dRects[k].getAttribute('width'));
                 var dh = parseFloat(dRects[k].getAttribute('height'));
                 if (!(dw > 0 && dh > 0)) continue;
                 if (dw * dh > dBest) {
                   dBest = dw * dh;
                   dlm = parseFloat(dRects[k].getAttribute('x'));
                   ddw = dw;
                 }
               }
               var mRects = mainSvg.querySelectorAll('defs clipPath rect');
               var mBestH = 0, mlm = NaN, mdw = NaN;
               for (var k = 0; k < mRects.length; k++) {
                 var mw = parseFloat(mRects[k].getAttribute('width'));
                 var mh = parseFloat(mRects[k].getAttribute('height'));
                 if (!(mw > 0 && mh > 0)) continue;
                 if (mh > mBestH) {
                   mBestH = mh;
                   mlm = parseFloat(mRects[k].getAttribute('x'));
                   mdw = mw;
                 }
               }
               if (isFinite(dlm) && isFinite(mlm) && ddw > 0 && mdw > 0) {
                 var rightX = dlm + ddw;
                 var scale = mdw / ddw;
                 detailSvg.style.transformOrigin = rightX + 'px 0';
                 detailSvg.style.transform = 'scaleX(' + scale + ')';
                 nn_status('  align d[' + Math.round(dlm) + '+' + Math.round(ddw) + ']' +
                           ' m[' + Math.round(mlm) + '+' + Math.round(mdw) + ']' +
                           ' s=' + scale.toFixed(3));
               } else {
                 nn_status('  align skipped: dlm=' + dlm + ' mlm=' + mlm +
                           ' ddw=' + ddw + ' mdw=' + mdw);
               }
             }
           }
         } catch (err) {
           nn_status('  ERROR ' + err.message);
         }
       };
       window.nn_hide_panel = function() {
         var c = document.getElementById('nn-detail-container');
         if (c) c.classList.remove('shown');
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

      page <- htmltools::tagList(
        nn_panel_css,
        nn_panel_js,
        top_meta_pin_js,
        top_wgt_fixed,
        top_meta_spacer,
        nn_wgt,
        detail_container,
        wgt
      )

      out_file <- fs::path(plot_dir, old_nn_input[["gene"]], ext = "html")
      htmltools::save_html(
        page,
        file = out_file,
        libdir = fs::path_file(fs::path(plot_dir, "dependency_files"))
      )
    } else {
      page <- htmltools::tagList(top_wgt, wgt)
      out_file <- fs::path(plot_dir, old_nn_input[["gene"]], ext = "html")
      htmltools::save_html(
        page,
        file = out_file,
        libdir = fs::path_file(fs::path(plot_dir, "dependency_files"))
      )
    }

  }, error = function(e) {
    message("make_protein_plot_win failed for gene ",
            old_nn_input[["gene"]], ": ", conditionMessage(e))
    message(paste(utils::limitedLabels(sys.calls()), collapse = "\n"))
  })

}
