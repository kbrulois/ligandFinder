




smoother_func <- function(x, window_size = 6, append_name = "a") {
  bind_cols(
    lapply(window_size, \(y) {
      tibble(!!paste0(append_name, y) := slider::slide_dbl(x, ~mean(., na.rm = TRUE), .before = y %/% 2, .after = y %/% 2))
    }
    ))
}



desc_colors <- c("#FED439FF", "#709AE1FF", "#FD7446FF",
                 "#D5E4A2FF", "#46732EFF", "#71D0F5FF",
                 "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF", "#8A9197FF",
                 "#197EC0FF", "#F05C3BFF")

color_scalify <- function(x,
                          colors = desc_colors
) {
  setNames(colors[1:length(x)], x)
}

ss_features <- c(setNames("Alpha helix (4-12)", "H"),
                 setNames("Isolated beta-bridge residue", "B"),
                 setNames("Strand", "E"),
                 setNames("3-10 helix", "G"),
                 setNames("Pi helix", "I"),
                 setNames("Turn", "T"),
                 setNames("Bend", "S"),
                 setNames("Kappa helix", "P"))

res_mod_lut <- c(setNames("P", "phosphorylation"),
                 setNames("S", "sulfation"),
                 setNames("A", "amidation"),
                 setNames("D", "deamidated"),
                 setNames("H", "hydroxylation"),
                 setNames("M", "methylation"),
                 setNames("N", "acetylation"),
                 setNames("G", "glycosylation"),
                 setNames("O", "other"))

res_mod_lut_inv <- names(res_mod_lut)
names(res_mod_lut_inv) <- unname(res_mod_lut)

topo_lut <-  c(setNames("signal peptide", "s"),
               setNames("intracellular", "i"),
               setNames("TM", "t"),
               setNames("extracellular", "e"),
               setNames("luminal", "l"),
               setNames("vesicular", "v"),
               setNames("data not available", "-"))

DBC <- c("NTC", "CTC")

all_desc_colors <- list(topo = color_scalify(c("s", "i", "e", "l", "-", "t", "", "v")),
                        SS = color_scalify(c("-", "H", "T", "P", "S", "E", "G", "B", "I")),
                        modification = color_scalify(unname(res_mod_lut)),
                        domain = color_scalify(c("Ig-like C2-type", "Cadherin", "EGF-like", "Fibronectin type-III",
                                                 "EGF-like; calcium-binding", "Sushi", "TSP type-1", "Ig-like",
                                                 "LDL-receptor class A", "Collagen-like", "Laminin EGF-like",
                                                 "Ig-like V-type", "CUB", "Laminin G-like", "Peptidase S1")),
                        NTC = color_scalify(DBC),
                        CTC = color_scalify(DBC)
)



tt_lut <- list(SS = ss_features,
               topo = topo_lut,
               modification = res_mod_lut_inv)

anno_feats <- tibble(feature = c("Dibasic",
                                 "Cysteine",
                                 "both",
                                 "afdsb",
                                 "unidsb",
                                 "phs",
                                 "gpcrdb_gtp",
                                 "sven",
                                 "top200NC",
                                 "disha",
                                 "uniprot_peptide",
                                 "new_pep"),
                     ggp = c("v_rec",
                             "v_rec",
                             "arch",
                             "arch",
                             "arch",
                             "h_rec",
                             "h_rec",
                             "h_rec",
                             "h_rec",
                             "h_rec",
                             "h_rec",
                             "h_rec"),
                     color = c("black",
                               "blue",
                               "#F05C3BFF",
                               "#197EC0FF",
                               "#C80813FF",
                               "#D5E4A2FF",
                               "#709AE1FF",
                               "#46732EFF",
                               "#FD7446FF",
                               "#C80813FF",
                               "#1A9993FF",
                               "#FED439FF")
)

c("#FED439FF", "#709AE1FF", "#FD7446FF",
  "#D5E4A2FF", "#46732EFF", "#71D0F5FF",
  "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF", "#8A9197FF",
  "#197EC0FF", "#F05C3BFF")

anno_feat_colors <- setNames(anno_feats[["color"]], anno_feats[["feature"]])

pep_nudges <- c(`new_pep` = 2.25, `disha` = 2.0, `top200NC` = 1.75, `uniprot_peptide` = 1.5, `phs` = 1.25, `gpcrdb_gtp` = 0.75, `sven` = 1)

assemble_anno_feats <- function(features, peps) {

  tmp <- features %>%
    as_tibble %>%
    mutate(type = if_else(type == "disulfide bond", "unidsb", type)) %>%
    filter(case_when(type %in% c("afdsb", "unidsb", "Dibasic", "Cysteine") ~ type %in% c("afdsb", "unidsb", "Dibasic", "Cysteine"),
                     source %in% c("phs", "gpcrdb_gtp", "sven", "top200NC", "disha", "uni_pep") ~ source %in% c("phs", "gpcrdb_gtp", "sven", "top200NC", "disha", "uni_pep"))) %>%
    mutate(feature = case_when(type %in% c("afdsb", "unidsb", "Dibasic", "Cysteine") ~ type,
                               source %in% c("phs", "gpcrdb_gtp", "sven", "top200NC", "disha", "uni_pep") ~ source), .before = everything()) %>%
    mutate(end = case_when(feature == "unidsb" & is.na(end) & stringr::str_detect(description, "Interchain") ~ stringr::str_extract(description, "C-\\d+") %>% stringr::str_remove(., "C-") %>% as.integer(.),
                           TRUE ~ end)) %>%
    filter(!is.na(end)) %>%
    select(feature, start, end)


  dsb <- tmp %>%
    filter(feature %in% c("afdsb", "unidsb")) %>%
    mutate(group = paste0(start, "_", end)) %>%
    summarise(feature = if(n() == 2) {"both"} else {first(feature)},
              start = first(start),
              end = first(end), .by = group) %>%
    select(-group)

  tmp <- bind_rows(tmp %>% filter(!feature %in% c("afdsb", "unidsb")),
                   dsb)

  if("start" %in% colnames(peps)) {
    tmp <- bind_rows(tmp, peps)
  }

  tmp %>%
    {left_join(., anno_feats, by = "feature")}

}


assemble_res_mod <- function(features, sequence_uni) {

  res_mod <- features %>%
    filter(type %in% c("glycosylation site", "modified residue")) %>%
    mutate(description = case_when(grepl("^Phospho", description) ~ "phosphorylation",
                                   grepl("^Sulfo", description) ~ "sulfation",
                                   grepl(" amide$", description) ~ "amidation",
                                   grepl("^Deamidated ", description) ~ "deamidated",
                                   grepl("hydroxy", description, ignore.case = TRUE) ~ "hydroxylation",
                                   grepl("methyl", description, ignore.case = TRUE) ~ "methylation",
                                   grepl("N-acetyl", description) ~ "acetylation",
                                   type == "glycosylation site" ~ "glycosylation",
                                   TRUE ~ "other"), .before = everything()) %>%
    mutate(type = res_mod_lut[description])

  tmp <- expand_features(res_mod,
                         sequence_uni = sequence_uni,
                         to_expand = tibble(type = c("anything"),
                                            source = c("uniprot"),
                                            filter_by = c("source")))
  if(ncol(tmp) == 0) {
    tmp <- tibble(uniprot_type = rep(NA, nchar(sequence_uni)))
  }
  tmp %>%
    select(!any_of("uniprot_desc")) %>%
    dplyr::rename(modification = uniprot_type)

}

assemble_dom <- function(features, sequence_uni) {

  dom <- features %>%
    filter(type == "domain") %>%
    mutate(type = stringr::str_remove(description, " \\d+"))

  tmp <- expand_features(dom,
                         sequence_uni = sequence_uni,
                         to_expand = tibble(type = c("anything"),
                                            source = c("uniprot"),
                                            filter_by = c("source")))

  if(ncol(tmp) == 0) {
    tmp <- tibble(uniprot_type = rep(NA, nchar(sequence_uni)))
  }
  tmp %>%
    select(!any_of("uniprot_desc")) %>%
    dplyr::rename(domain = uniprot_type)

}

assemble_SV <- function(features, sequence_uni) {

  SV <- features %>%
    filter(type == "sequence variant")

  tmp <- expand_features(SV,
                         sequence_uni = sequence_uni,
                         to_expand = tibble(type = c("anything"),
                                            source = c("uniprot"),
                                            filter_by = c("source")))

  if(ncol(tmp) == 0) {
    return(tibble(SV = rep(NA, nchar(sequence_uni))))
  } else {
    return(tmp %>%
             mutate(SV = if_else(!is.na(uniprot_desc) & stringr::str_detect(uniprot_desc, "dbSNP:rs\\d+"), stringr::str_extract(uniprot_desc, "dbSNP:rs\\d+"), uniprot_type))
    )
  }
}

assemble_DBC <- function(features, sequence_uni) {

  dbc <- features %>%
    filter(source == "dbc") %>%
    group_by(type, start) %>%
    summarise(start = start[!is.na(start)][1]) %>%
    mutate(type = stringr::str_extract(type, "^NTC|^CTC")) %>%
    mutate(source = type) %>%
    mutate(end = NA) %>%
    mutate(description = NA) %>%
    ungroup()

  tmp <- expand_features(dbc,
                         sequence_uni = sequence_uni,
                         to_expand = tibble(type = c("anything"),
                                            source = unique(dbc$source),
                                            filter_by = c("source")))
  if(ncol(tmp) == 0) {
    tmp <- tibble(NTC = rep(NA, nchar(sequence_uni)),
                  CTC = rep(NA, nchar(sequence_uni)))
  }
  tmp %>%
    dplyr::rename_with(.fn = ~sub("_type", "", .))


}

expand_by_residue <- function(x,
                              dat_to_expand = c("ref_index",
                                                "topo",
                                                "features_expanded",
                                                "modification",
                                                "domain",
                                                "SV",
                                                "DBC",
                                                "dssp",
                                                "af_missense",
                                                "cons",
                                                "alignment_AA",
                                                "aa_scores")
) {

  x <- x %>%
    mutate(features_expanded = map2(features, sequence_uni, expand_features)) %>%
    mutate(modification = map2(features, sequence_uni, assemble_res_mod)) %>%
    mutate(domain = map2(features, sequence_uni, assemble_dom)) %>%
    mutate(SV = map2(features, sequence_uni, assemble_SV)) %>%
    mutate(DBC = map2(features, sequence_uni, assemble_DBC)) %>%
    mutate(ref_index = map(sequence_uni, \(x) {tibble(index = 1:nchar(x),
                                                      AA = stringr::str_split(x, "", simplify = TRUE) %>% `c`)}))

  x <- x %>%
    mutate(topo = map2(sequence_uni, topo, \(x, y) tibble(AA = stringr::str_split(x, "", simplify = TRUE) %>% c,
                                                          topo = stringr::str_split(y, "", simplify = TRUE) %>% c))) %>%
    mutate(to_expand = pmap(pick(any_of(dat_to_expand)),
                            bind_cols, .name_repair = "minimal")) %>%
    mutate(to_expand = map(to_expand, \(x) x[, !duplicated(colnames(x))])) %>%
    mutate(to_expand = map(to_expand, \(x) x[, colnames(x) != ""]))

  to_return <- x %>%
    select(!where(is.list), to_expand) %>%
    unnest(to_expand)

  to_return <- to_return %>%
    mutate(topo2 = if_else(has_topo & ("e" %in% topo), topo, "e"))

  return(to_return)

}


expand_features <- function(features,
                            sequence_uni,
                            fill_type = c("type", "desc"),
                            to_expand = tibble(type = c("gpcr_pep", "sven_pep", "peptide", "sequence variant", "signal peptide", "Dibasic"),
                                               source = c("gpcrdb_gtp", "sven", "uniprot", "uniprot", "uniprot", "sites"),
                                               filter_by = c("source", "source", "type_source", "type_source", "type_source", "type_source"))
) {

  dat_size <- nchar(sequence_uni)

  feat_clean <- features %>%
    filter(!is.na(start)) %>%
    filter(start <= dat_size) %>%
    filter(is.na(end) | end <= dat_size)

  dat_ts <- dat_s <- tibble(.rows = dat_size)

  to_expand_s <- to_expand %>%
    filter(filter_by == "source")

  if(nrow(to_expand_s) > 0) {

    cols <- to_expand_s %>% pull(source) %>% unique()

    col_names <- expand.grid(cols = cols,
                             fill_type = fill_type) %>%
      mutate(cols = paste0(cols, "_", fill_type)) %>%
      pull(cols)

    dat_s[col_names] <- NA

    for(x in cols) {

      feat_sub <- feat_clean %>%
        filter(source == x)

      for(i in seq_len(nrow(feat_sub))) {

        start_val <- feat_sub[["start"]][i]
        end_val <- feat_sub[["end"]][i]

        if(is.na(end_val)) {
          dat_s[[paste0(x, "_type")]][start_val] <- feat_sub[["type"]][i]
          dat_s[[paste0(x, "_desc")]][start_val] <- feat_sub[["description"]][i]
        } else {
          dat_s[[paste0(x, "_type")]][start_val:end_val] <- feat_sub[["type"]][i]
          dat_s[[paste0(x, "_desc")]][start_val:end_val] <- feat_sub[["description"]][i]
        }
      }

    }

    dat_s <- dat_s %>%
      select(where(~ !all(is.na(.x))))
  }

  #####different type of expansion
  to_expand_ts <- to_expand %>%
    filter(filter_by == "type_source")

  if(nrow(to_expand_ts) > 0) {

    to_expand_ts <- to_expand_ts %>%
      distinct(source, type)

    cols <- to_expand_ts %>%
      mutate(type2 = paste0(source, "_", type)) %>%
      pull(type2)

    col_names <- expand.grid(cols = cols,
                             fill_type = fill_type) %>%
      mutate(cols = paste0(cols, "_", fill_type)) %>%
      pull(cols)

    dat_ts[col_names] <- NA

    for(x in seq_along(cols)) {

      type_val <- to_expand_ts %>% dplyr::slice(x) %>% pull(type)
      source_val <- to_expand_ts %>% dplyr::slice(x) %>% pull(source)
      x <- cols[x]

      feat_sub <- feat_clean %>%
        filter(type == type_val & source == source_val)

      for(i in seq_len(nrow(feat_sub))) {

        start_val <- feat_sub[["start"]][i]
        end_val <- feat_sub[["end"]][i]

        if(is.na(end_val)) {
          dat_ts[[paste0(x, "_type")]][start_val] <- feat_sub[["type"]][i]
          dat_ts[[paste0(x, "_desc")]][start_val] <- feat_sub[["description"]][i]
        } else {
          dat_ts[[paste0(x, "_type")]][start_val:end_val] <- feat_sub[["type"]][i]
          dat_ts[[paste0(x, "_desc")]][start_val:end_val] <- feat_sub[["description"]][i]
        }
      }

    }

    dat_ts <- dat_ts %>%
      select(where(~ !all(is.na(.x))))
  }

  return(bind_cols(dat_s, dat_ts))

}




make_protein_plot <- function(input,
                              pep_tp
) {

  tryCatch({
    to_plot <- expand_by_residue(input)


    to_plot <- to_plot %>%
      group_by(gene) %>%
      mutate(max = max(index)) %>%
      mutate(across(starts_with(c("pep_nn", "pep_xgb", "chem_nn", "chem_xgb")), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE))

    transfrom_aligments <- function(x, vals_to = "AA") {
      x %>%
        mutate(index = row_number()) %>%
        pivot_longer(cols = -index, names_to = "metric", values_to = vals_to)
    }

    to_plot_aln <- tryCatch({input %>%
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

    features <- input %>% pull(features) %>% `[[`(1)
    AA_sequence <- input %>% pull(sequence_uni) %>% stringr::str_split(., "", simplify = TRUE) %>% `c`

    seq_dat <- bind_rows(tibble(AA = NA, index = 0, index_tp = NA),
                         tibble(AA = AA_sequence) %>%
                           mutate(index = row_number()) %>%
                           mutate(index_tp = if_else(index %in% ind_tp, index, NA))
    )

    anno_feat_dat <- assemble_anno_feats(features, peps = pep_tp)

    anno1 <- input %>%
      pull(annotations) %>%
      `[[`(1) %>%
      filter(annotation_name == "comment") %>%
      filter(annotation_type %in% c("subcellular location", "tissue specificity", "disease")) %>%
      group_by(annotation_type) %>%
      summarise(annotation = list(paste0(annotation, collapse = "; "))) %>%
      pivot_wider(names_from = annotation_type, values_from = annotation)

    anno2 <- input %>%
      pull(annotations) %>%
      `[[`(1) %>%
      filter(annotation_name == "dbReference" & name_1 == "disease") %>%
      group_by(annotation_type) %>%
      summarise(annotation = paste0(annotation, collapse = "; ")) %>%
      {paste0(.[["annotation_type"]], ": ", .[["annotation"]])}



    p_title <- input %>%
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
                           #color = "black",
                           #lwd = 0.1,
                           width = 1,
                           height = 0.5,
                           position = position_nudge(y = -0.25)) +

        ggplot2::geom_tile(data = cons_dat,
                           mapping = aes(x = index, y = metric, fill = blos),
                           #color = "black",
                           #lwd = 0.1,
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
      mutate(nudge = unname(pep_nudges[feature])) %>%
      mutate(tt_value = paste0(feature))

    top_y <- df %>% filter(metric_type == "") %>% mutate(metric = droplevels(metric)) %>% pull(metric) %>% levels %>% length


    main_p <- main_p + ggiraph::geom_segment_interactive(
      data = h_rec_dat,
      aes(x = start - 0.3, xend = end + 0.3, y = top_y + nudge, yend = top_y + nudge, color = feature, tooltip = tt_value),
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
              seq_p_index = seq_p_index,
              title_p = title_text,
              subtitle_p = subtitle_text,
              links = p_title)

    num_res = df %>% count(metric, gene) %>% group_by(gene) %>% summarise(max = max(n, na.rm = FALSE)) %>% pull(max, gene)
    num_mets = df %>% count(index, gene) %>% group_by(gene) %>% summarise(max = max(n, na.rm = FALSE)) %>% pull(max, gene)


    p_width <- num_res/7
    p_height <- (num_mets  + 3)/4.3

    final_p <- p[["meta_p"]] /  p[["seq_p"]]  / p[["plot"]]

    final_p <- final_p + patchwork::plot_layout(heights = c(0.5, 0.7, p_height)) &
      theme(panel.spacing = unit(0, "pt"),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 80)  # l = 50 increases left margin
      )
    wgt <- ggiraph::girafe(ggobj = final_p,
                           width_svg = p_width,
                           height_svg = p_height,
                           options = list(
                             ggiraph::opts_sizing(rescale = FALSE),
                             #ggiraph::opts_hover(css = "stroke-width: 3px; transition: all 0.3s ease;"),
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

    # -----------------------------
    # JAVASCRIPT
    # -----------------------------

    wgt <- htmlwidgets::onRender(
      wgt,
      paste0("
           const GENE   = '", input[["gene"]], "';
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

           // ---- TIGHT RECT DRAWING ----
             function drawSelectionRect(start, end, color) {
  const t1 = getTextEl(start);
  const t2 = getTextEl(end);
  if (!t1 || !t2) return null;

  const b1 = t1.getBBox();
  const b2 = t2.getBBox();

  const padY  = 1;   // vertical padding
  const trimX = 1;   // horizontal tightening

  // Horizontal bounds
  const x1 = Math.min(b1.x, b2.x);
  const x2 = Math.max(b1.x + b1.width, b2.x + b2.width);

  // Vertical visual center (baseline-safe)
  const top    = Math.min(b1.y, b2.y);
  const bottom = Math.max(b1.y + b1.height, b2.y + b2.height);
  const centerY = (top + bottom) / 2;
  const height  = (bottom - top) + padY * 2;

  const rect = document.createElementNS(
    'http://www.w3.org/2000/svg',
    'rect'
  );

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

           // -----------------------------
             // LOAD CXC + RESTORE RECTANGLES
           // -----------------------------
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

           // -----------------------------
             // UNDO (REPEATABLE)
           // -----------------------------
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

           // -----------------------------
             // CLICK SELECTION
           // -----------------------------
             document.addEventListener('click', function(e) {
               const t = e.target;
               if (!t || !t.hasAttribute('data-id')) return;

               const idx = +t.getAttribute('data-id');
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

           // -----------------------------
             // SAVE
           // -----------------------------
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

           // -----------------------------
             // UI
           // -----------------------------
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


    htmlwidgets::saveWidget(widget = wgt,
                            file = fs::path(plot_dir, input[["gene"]], ext = "html"),
                            selfcontained = TRUE,
                            libdir = fs::path(plot_dir, "dependency_files"),
                            title = input[["gene"]])

  }, error = function(e) message(e[["message"]]))


}













plot_dir <- "~/AF2_analysis/all_peps_v4/"


mets <- list(cons = c("blos_wt_all_n", "cons_rs", "blos_wt_mam", "blos_wt_all", "gran_wt_all"),
             af_missense = c("mean_afm", "min_afm"),
             dssp = c("relASA"),
             aa_scores = c("pep_xgb4c", "chem_xgb3c", "pep_nn4c", "chem_nn4c")
)

mets$aa_scores <- c(mets$aa_scores, paste0(mets$aa_scores, "_s6"))

all_mets <- do.call(`c`, mets) %>% unname

names(all_mets) <- rep(names(mets), sapply(mets, length))


secretome <- readRDS("~/AF2_analysis/secretome_latest.rds.rds")

peps_tp <- readRDS("~/AF2_analysis/peps_tp.rds")


interv <- 1:2


dir.create(plot_dir)
unlink(plot_dir)

genes <- c("ANO8", "GDF5")


the_input <- secretome %>%
  filter(!accession %in% c("A0AAG2TCD0", "A0AAG2UXZ5")) %>%
  filter(gene %in% !!genes) %>%
  mutate(aa_scores = map(aa_scores, \(x) x[, colnames(x) %in% mets[["aa_scores"]]])) %>%
  group_split(gene)

pep_input <- peps_tp %>%
  filter(gene %in% genes) %>%
  group_split(gene)

input_gene <- map_chr(the_input, \(x) x$gene)

pep_input_gene <- map_chr(pep_input, \(x) x$gene)

identical(input_gene, pep_input_gene)





to_iterate <- c(seq.default(1, length(input_gene), by = 100), length(input_gene) + 1)

chunks <- length(to_iterate) - 1

for(x in 56:chunks) {

  interv <- to_iterate[x]:(to_iterate[x + 1] - 1)

  message("computing ", x, " of ", chunks)

  start <- Sys.time()
  #future::plan(strategy = future::multisession(workers = 6))
  future::plan(strategy = future::sequential())


  purrr::walk2(.x = the_input[interv],
               .y = pep_input[interv],
               make_protein_plot)

  end <- Sys.time()

  end - start

  gc()


}


input <- the_input[[1]]
pep_tp <- pep_input[[1]]












