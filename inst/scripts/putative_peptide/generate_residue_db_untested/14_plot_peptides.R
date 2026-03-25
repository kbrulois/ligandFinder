



secretome_roi <- readRDS("~/AF2_analysis/secretome_roi.rds")

new_peps <- list.files("~/AF2_analysis/new_peps_html") %>% stringr::str_remove(., ".html$")


pepoi <- tibble(gene = c("SCG3BA2", "ANO8"),
                start = c(70,
                end = c(82)))

secretome_roi <- secretome_roi %>%
                    mutate(link = paste0("https://stacks.stanford.edu/file/mh260xs0996/all_peps_v2/all_peps_v2/", gene, ".html")) %>%
                    mutate(selection = case_when(pep_xgb4c_max > 0.6 & pep_nn4c_max > 0.6 ~ "both",
                                                 pep_xgb4c_max > 0.6 & !pep_nn4c_max > 0.6 ~ "xgb",
                                                 !pep_xgb4c_max > 0.6 & pep_nn4c_max > 0.6 ~ "nn",
                                                 TRUE ~ "low_score")) %>%
                    mutate(roi_type2 = if_else(roi_type == "phs" & selection == "low_score", "phs_low_score", roi_type)) %>%
                    mutate(roi_type3 = case_when(roi_type2 == "phs" & `percent_ol_phs:gpcrdb_gtp` > 0 ~ "phs_known",
                                                 roi_type2 == "phs" & `percent_ol_phs:uni_pep` > 50 ~ "phs_ol_uniprot_peptide",
                                                 roi_type2 == "phs" & `percent_ol_phs:disha` > 20 ~ "phs_ol_disha",
                                                 roi_type2 == "phs" & `percent_ol_phs:disha` == 0 ~ "phs_no_ol",
                                                 roi_type2 == "phs" ~ "phs_no_ol",
                                                 TRUE ~ roi_type2)) %>%
                    mutate(location2 = ifelse(stringr::str_detect(location, "l$"), "Sec", "TM")) %>%
                    mutate(db_status = case_when(has_db3_start & has_db3_end ~ "flankingDB",
                                                 has_db3_start | has_db3_end ~ "singleDB",
                                                 !has_db3_start & !has_db3_end ~ "DBfree")) %>%
                    mutate(db_status = case_when(db_status != "singleDB" ~ db_status,
                                                 has_db3_start & hsterm == "N" ~ "singleDBmatched",
                                                 has_db3_end & hsterm == "C" ~ "singleDBmatched",
                                                 TRUE ~ "singleDBunmatched"))

data.table::fwrite(secretome_roi %>% select(!where(is.list)), "~/AF2_analysis/roi_initial8.csv")






secretome_roi <- data.table::fread("~/AF2_analysis/roi_initial8.csv") %>% as_tibble



secretome_roi <- secretome_roi %>%
                    mutate(id = paste0(gene,"_", start, "x", end)) %>%
                    mutate(link = paste0("https://stacks.stanford.edu/file/rn350qj0909/all_peps_v5/", gene, ".html")) %>%
                    mutate(cxc_link = paste0("https://stacks.stanford.edu/file/tc396gg4330/v1/", gene, ".cxc")) %>%
                    mutate(roi_type3 = factor(roi_type3, levels = c("phs_no_ol", "phs_ol_disha", "disha", "phs_ol_uniprot_peptide", "gpcrdb_gtp", "sven", "top200NC", "phs_known", "phs_low_score"))) %>%
                    mutate(db_status = factor(db_status, levels = c("flankingDB", "singleDBmatched", "singleDBunmatched", "DBfree"))) %>%
                    mutate(db_score = db_status == "flankingDB") %>%
                    mutate(dbc_score = max_dbc > 0.5) %>%
                    mutate(ol_with_uni = `percent_ol_uni_pep:phs` > 0) %>%
                    mutate(sel_score = selection == "both") %>%
                    mutate(no_cys = !has_cys) %>%
                    mutate(final_score = no_cys + sel_score + ol_with_uni + dbc_score + db_score + has_cterm_amid) %>%
                    mutate(tt_value = paste(paste0(gene),
                                            paste0(start, "-", end),
                                            paste0("has_TM: ", has_TM),
                                            paste0("has_cys: ", has_cys),
                                            paste0("has_cterm_amid: ", has_cterm_amid),
                                            paste0("Localization: ", location2, "; ", location),
                                            paste0("NTC_tot_start: ", round(NTC_tot_start, 2)),
                                            paste0("CTC_tot_end: ", round(CTC_tot_end, 2)),
                                            paste0("nn_max: ", round(pep_nn4c_max, 2)),
                                            paste0("nn_c: ", round(pep_nn4c_cterm, 2)),
                                            paste0("nn_n: ", round(pep_nn4c_nterm, 2)),
                                            paste0("nn: ", round(pep_nn4c_entire, 2)),
                                            paste0("xgb_max: ", round(pep_xgb4c_max, 2)),
                                            paste0("xgb_c: ", round(pep_xgb4c_cterm, 2)),
                                            paste0("xgb_n: ", round(pep_xgb4c_nterm, 2)),
                                            paste0("xgb: ", round(pep_xgb4c_entire, 2)),
                                            paste0("PepRank: ", round(PeptideRanker_score, 2)),
                                            paste0("db_per_AA: ", round(db_per_AA, 2)),
                                            paste0("db_status: ", db_status), sep = "\n"))






plot_dir <- "~/Desktop/tmp2"

dir.create(plot_dir)

params <-   c("rank", colnames(secretome_roi)[c(12, 15, 132, 133, 134, 139, 208:211, 240:243, 249, 250, 312, 315, 361, 325, 329, 374, 379, 378, 393, 396)])


for(param in params) {

  df <- secretome_roi %>%
    filter(!is.na(pep_xgb4c_max) & !is.na(db_status)) %>%
    filter(db_status %in% c("flankingDB", "singleDBmatched")) %>%
    filter(roi_length > 4)

  df2 <- df %>%
      filter(roi_type == "phs") %>%
      filter(selection != "low_score") %>%
      filter(db_status %in% c("flankingDB", "singleDBmatched")) %>%
      filter(`percent_ol_phs:gpcrdb_gtp` == 0 & `percent_ol_phs:top200NC` < 10) %>%
      arrange(desc(final_score), desc(best_model_entire)) %>%
      mutate(rank = row_number()) %>%
      slice(1:200) %>%
      arrange(id) %>%
      mutate(plot = cut(row_number(), 4, labels = 1:4))

  if(param == "rank") {
    df <- left_join(df, df2 %>% select(id, rank), by = "id")
  }

  df <- df %>%
    arrange(!!rlang::sym(param))


p <- ggplot2::ggplot(df[order(df[[param]], na.last = FALSE), ]) +
  ggiraph::geom_point_interactive(aes(x = UMAP1,
                                      y = UMAP2,
                                      color = !!rlang::sym(param),
                                      data_id = id,
                                      tooltip = tt_value,
                                      onclick = paste0('window.open("', link , '");', 'window.open("', cxc_link , '")'))) +
  ggh4x::facet_nested(cols = vars(db_status), rows = vars(roi_type3), switch = "y")


if(param %in% c("hsterm", "has_TM", "has_cterm_amid", "has_cys", "ol_with_uni")) {
  p <- p +
  ggplot2::scale_color_discrete(palette = ggsci::pal_simpsons())
} else {
  p <- p +
    ggplot2::scale_color_viridis_c(option = "H")
}

  p <- p +
  theme_bw()


df2 <- df2  %>%
    group_split(plot)

p2 <- map(df2, \(to_plot) {
            all_id <- to_plot[["id"]]
            p <- ggplot2::ggplot(data = to_plot %>%
                                   mutate(index = 1) %>%
                                   mutate(id = factor(id, levels = rev(id)))) +
                  ggiraph::geom_tile_interactive(aes(x = index,
                                                     y = id,
                                                     fill = !!rlang::sym(param),
                                                     data_id = id,
                                                     tooltip = tt_value,
                                                     onclick = paste0('window.open("', link , '");', 'window.open("', cxc_link , '")')),
                                                 width = 0.6,
                                                 height = 0.6,
                                                 color = "black",
                                                 stroke = 0.2)

  if(param %in% c("hsterm", "has_TM")) {
    p <- p +
      ggplot2::scale_fill_discrete(palette = ggsci::pal_simpsons(), guide = "none") +
      ggplot2::scale_y_discrete(position = "right")
  } else {
    p <- p +
      ggplot2::scale_fill_viridis_c(option = "H", guide = "none") +
      ggplot2::scale_y_discrete(position = "right")
  }

  p <- p +
    theme_void() +
    theme(axis.text.y.right = element_text(size = 9, hjust = 0.1))
  p
})



final_p <- Reduce(`+`, p2) + p

design <-
"
12345
####5"


final_p <- final_p + patchwork::plot_layout(design = design, widths = c(1,1,1,1, 30), guides = "collect") &
  theme(panel.spacing = unit(0, "pt"),
        legend.position = "top")

wgt <- ggiraph::girafe(ggobj = final_p,
                       width_svg = 11,
                       height_svg = 15,
                       options = list(
                         ggiraph::opts_sizing(rescale = FALSE),
                         ggiraph::opts_hover(css = "r:8px;stroke:black;stroke-width:2px;"), # Increase radius to 8px and add a border
                         ggiraph::opts_tooltip(css = "font-family:sans-serif;color:black;font-size:big;background-color:white;") # Optional: customize tooltip style

                       ))

htmlwidgets::saveWidget(widget = wgt,
                        file = fs::path(plot_dir, param, ext = "html"),
                        selfcontained = TRUE)


}








































plot_dir <- "~/Desktop/tmp2"

dir.create(plot_dir)

params <-   c("rank", colnames(secretome_roi)[c(12, 15, 132, 133, 134, 139, 208:211, 240:243, 249, 250, 312, 315, 361, 325, 329, 374, 379, 378, 393, 396)])


param <- "pep_nn4c_max"

  df <- secretome_roi %>%
    filter(!is.na(pep_xgb4c_max) & !is.na(db_status)) %>%
    filter(db_status %in% c("flankingDB", "singleDBmatched")) %>%
    filter(roi_length > 4)

  df2 <- df %>%
    filter(roi_type == "phs") %>%
    filter(selection != "low_score") %>%
    filter(db_status %in% c("flankingDB", "singleDBmatched")) %>%
    filter(`percent_ol_phs:gpcrdb_gtp` == 0 & `percent_ol_phs:top200NC` < 10) %>%
    arrange(desc(final_score), desc(best_model_entire)) %>%
    mutate(rank = row_number()) %>%
    slice(1:400) %>%
    arrange(id) %>%
    mutate(plot = cut(row_number(), 8, labels = 1:8))

  if(param == "rank") {
    df <- left_join(df, df2 %>% select(id, rank), by = "id")
  }



  df2 <- df2  %>%
    group_split(plot)

  p2 <- map(df2, \(to_plot) {
    all_id <- to_plot[["id"]]
    p <- ggplot2::ggplot(data = to_plot %>%
                           mutate(index = 1) %>%
                           mutate(id = factor(id, levels = rev(id)))) +
      ggiraph::geom_tile_interactive(aes(x = index,
                                         y = id,
                                         fill = !!rlang::sym(param),
                                         data_id = id,
                                         tooltip = tt_value,
                                         onclick = paste0('window.open("', link , '");', 'window.open("', cxc_link , '")')),
                                     width = 0.6,
                                     height = 0.6,
                                     color = "black",
                                     stroke = 0.2)

    if(param %in% c("hsterm", "has_TM")) {
      p <- p +
        ggplot2::scale_fill_discrete(palette = ggsci::pal_simpsons(), guide = "none") +
        ggplot2::scale_y_discrete(position = "right")
    } else {
      p <- p +
        ggplot2::scale_fill_viridis_c(option = "H", guide = "none") +
        ggplot2::scale_y_discrete(position = "right")
    }

    p <- p +
      theme_void() +
      theme(axis.text.y.right = element_text(size = 9, hjust = 0.1))
    p
  })



  final_p <- Reduce(`+`, p2)

  design <-
    "
12345678
"


  final_p <- final_p + patchwork::plot_layout(design = design, widths = c(1,1,1,1, 1,1), guides = "collect") &
    theme(panel.spacing = unit(0, "pt"),
          legend.position = "top")

  wgt <- ggiraph::girafe(ggobj = final_p,
                         width_svg = 11,
                         height_svg = 8,
                         options = list(
                           ggiraph::opts_sizing(rescale = FALSE),
                           ggiraph::opts_hover(css = "r:8px;stroke:black;stroke-width:2px;"), # Increase radius to 8px and add a border
                           ggiraph::opts_tooltip(css = "font-family:sans-serif;color:black;font-size:big;background-color:white;") # Optional: customize tooltip style

                         ))

  htmlwidgets::saveWidget(widget = wgt,
                          file = paste0("~/Desktop/top400_", param, ".html"),
                          selfcontained = TRUE)



























p <- ggplot2::ggplot(secretome_roi %>% filter(!is.na(pep_xgb4c_max) & !is.na(db_status)) %>% filter(roi_length > 4)) +
  ggiraph::geom_point_interactive(aes(x = UMAP1,
                                      y = UMAP2,
                                      color = selection,
                                      data_id = id,
                                      tooltip = tt_value,
                                      onclick = paste0('window.open("', link , '");', 'window.open("', cxc_link , '")'))) +
  ggh4x::facet_nested(rows = vars(db_status), cols = vars(roi_type3), switch = "y") +
  #ggplot2::scale_color_viridis_c(option = "H") +
  ggplot2::scale_color_manual(values = c(`both` = "#FD7446FF", `low_score` = "#8A9197FF", `nn` = "#709AE1FF", `xgb` = "#FED439FF")) +
  theme_bw()

wgt <- ggiraph::girafe(ggobj = p,
                       width_svg = 24,
                       height_svg = 12,
                       options = list(
                         ggiraph::opts_sizing(rescale = FALSE)
                       ))

htmlwidgets::saveWidget(widget = wgt,
                        file = fs::path("~/Desktop/test2", ext = "html"),
                        selfcontained = TRUE)







