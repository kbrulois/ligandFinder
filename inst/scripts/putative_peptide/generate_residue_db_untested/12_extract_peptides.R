



all_feat_cols <- unique(do.call(c, lapply(secretome$features, \(x) colnames(x))))



extract_roi <- \(x, y) {
  rois <- x %>%
    filter(source %in% c("gpcrdb_gtp", "sven", "phs")) %>%
    dplyr::select(-start, -end, -evidence) %>%
    dplyr::rename("roi_name" = "type", "roi_type" = "source")

  dummy_table <- tibble(accession = y, roi_name = NA, roi_type = NA)
  if(nrow(rois) > 0) {
    return(bind_cols(accession = rep(y, nrow(rois)), rois))
  } else {
    return(dummy_table)
  }
}

left_side <- secretome %>%
  rowwise() %>%
  reframe(extract_roi(features, accession))

secretome_roi <- left_join(left_side, secretome %>% select(!where(is.list)), by = "accession")

secretome_roi <- secretome_roi %>%
  mutate(known_ligand = droplevels(factor(ifelse(roi_type == "gpcrdb_gtp", "known", "putative")))) %>%
  mutate(SV = if_else(`SV_sequence variant_entire` > 0, "yes", "no")) %>%
  mutate(across(starts_with("percent_ol_"), ~replace_na(., replace = 0)))








subsetter <- secretome_roi$roi_length < 3 & !is.na(secretome_roi$roi_length)

to_umap <- secretome_roi[!subsetter, ] %>%
  dplyr::select(c(17:249, 315:318)) %>%
  dplyr::select(!where(is.character)) %>%
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

secretome_roi[["UMAP1"]] <- NA
secretome_roi[["UMAP2"]] <- NA

secretome_roi[["UMAP1"]][!subsetter] <- umap_res[["layout"]][,1]
secretome_roi[["UMAP2"]][!subsetter] <- umap_res[["layout"]][,2]



data.table::fwrite(secretome_roi %>% select(!where(is.list)), "~/AF2_analysis/roi_initial5.csv")










umap_gate <- readRDS("~/AF2_analysis/peptide_umap_gate.rds")

goi <- c("ANGPTL2", "ANO8", "COCH", "CRTAP", "DHH", "ENAM", "GRM4",
         "ITIH2", "KCNV2", "MCEMP1", "NUCB1", "PLVAP", "SHISA9", "SLC39A6",
         "SMOC2", "TMC7", "TRABD2A", "WNT1")

secretome_roi[["umap_gate"]] <- umap_gate

tm_focused <- secretome_roi %>%
  #filter(umap_gate != "no") %>%
  filter(!location %in% c("4l", "3l", "2l")) %>%
  filter(roi_length > 10) %>%
  filter(db_per_AA > 0.01) %>%
  filter(roi_type == "phs") %>%
  filter(pep_xgb4c_max > 0.6 | pep_nn4c_max > 0.6) %>%
  filter(`percent_ol_phs:gpcrdb_gtp` == 0) %>%
  filter(`percent_ol_phs:uni_pep` == 0) %>%
  filter(`percent_ol_phs:top200NC` == 0) %>%
  filter(has_db3_start & has_db3_end) %>%
  filter(CTC_tot_end > 0.5 | NTC_tot_start > 0.5) %>%
  #filter(dbc_term_start == "N" & dbc_term_end == "C") %>%
  select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end) %>%
  mutate(selection = "tm_focused")


dibasic_context <- secretome_roi %>%
  #filter(umap_gate != "no") %>%
  #filter(location %in% c("4l", "3l", "2l")) %>%
  filter(roi_length > 4) %>%
  filter(db_per_AA > 0.01) %>%
  filter(roi_type == "phs") %>%
  filter(pep_xgb4c_max > 0.6 | pep_nn4c_max > 0.6) %>%
  filter(`percent_ol_phs:gpcrdb_gtp` == 0) %>%
  filter(`percent_ol_phs:uni_pep` == 0) %>%
  filter(`percent_ol_phs:top200NC` == 0) %>%
  filter(has_db3_start | has_db3_end) %>%
  #filter(CTC_tot_end > 0.5 | NTC_tot_start > 0.5) %>%
  #filter(dbc_term_start == "N" & dbc_term_end == "C") %>%
  select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end) %>%
  mutate(selection = "dibasic_context")


secretoglobins <- secretome_roi %>%
  #filter(umap_gate != "no") %>%
  filter(grepl("^SCG", gene)) %>%
  filter(roi_length > 10) %>%
  #filter(db_per_AA > 0.01) %>%
  filter(roi_type %in% c("phs")) %>%
  filter(pep_xgb4c_max > 0.6 | pep_nn4c_max > 0.6 | CTC_tot_end > 0.5 | NTC_tot_start > 0.5) %>%
  filter(`percent_ol_phs:gpcrdb_gtp` == 0) %>%
  filter(`percent_ol_phs:uni_pep` == 0) %>%
  filter(`percent_ol_phs:top200NC` == 0) %>%
  filter(has_db3_start | has_db3_end) %>%
  #filter(dbc_term_start == "N" & dbc_term_end == "C") %>%
  select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end) %>%
  arrange(desc(pep_xgb4c_max)) %>%
  mutate(selection = "secretoglobin")



umap_gated <- secretome_roi %>%
  filter(umap_gate != "no") %>%
  #filter(!location %in% c("4l", "3l", "2l")) %>%
  filter(roi_length > 10) %>%
  filter(db_per_AA > 0.01) %>%
  filter(roi_type == "phs") %>%
  filter(pep_xgb4c_max > 0.6 | pep_nn4c_max > 0.6) %>%
  filter(`percent_ol_phs:gpcrdb_gtp` == 0) %>%
  filter(`percent_ol_phs:uni_pep` == 0) %>%
  filter(`percent_ol_phs:top200NC` == 0) %>%
  filter(has_db3_start & has_db3_end) %>%
  #filter(CTC_tot_end > 0.5 | NTC_tot_start > 0.5) %>%
  #filter(dbc_term_start == "N" & dbc_term_end == "C") %>%
  select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end) %>%
  mutate(selection = "umap_gated")


c_term_amid <- secretome_roi %>%
  #filter(umap_gate != "no") %>%
  #filter(!location %in% c("4l", "3l", "2l")) %>%
  filter(roi_length > 10) %>%
  filter(db_per_AA > 0.01) %>%
  filter(roi_type == "phs") %>%
  filter(pep_xgb4c_max > 0.6 | pep_nn4c_max > 0.6) %>%
  filter(`percent_ol_phs:gpcrdb_gtp` == 0) %>%
  filter(`percent_ol_phs:uni_pep` == 0) %>%
  filter(`percent_ol_phs:top200NC` == 0) %>%
  filter(has_db3_start & has_db3_end) %>%
  filter(grepl("G(?:KK|KR|RK|RR)", seq_dbcenter_end)) %>%
  #filter(CTC_tot_end > 0.5 | NTC_tot_start > 0.5) %>%
  #filter(dbc_term_start == "N" & dbc_term_end == "C") %>%
  select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end, location) %>%
  mutate(selection = "c_term_amid")


single_db <- secretome_roi %>%
  #filter(umap_gate != "no") %>%
  #filter(location %in% c("4l", "3l", "2l")) %>%
  filter(roi_length > 4) %>%
  filter(db_per_AA > 0.01) %>%
  filter(roi_type == "phs") %>%
  filter(pep_xgb4c_max > 0.6 | pep_nn4c_max > 0.6) %>%
  filter(`percent_ol_phs:gpcrdb_gtp` == 0) %>%
  filter(`percent_ol_phs:uni_pep` == 0) %>%
  filter(`percent_ol_phs:top200NC` == 0) %>%
  filter(has_db3_start | has_db3_end) %>%
  #filter(CTC_tot_end > 0.5 | NTC_tot_start > 0.5) %>%
  #filter(dbc_term_start == "N" & dbc_term_end == "C") %>%
  select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end) %>%
  mutate(selection = "single_db")





peps <- bind_rows(secretoglobins,
                  tm_focused,
                  dibasic_context,
                  umap_gated,
                  c_term_amid,
                  single_db)

peps <- peps %>%
        mutate(roi_name %>%
                 stringr::str_remove(., "^phs_") %>%
                 tibble_split(., "-", names = c('start', 'end'))) %>%
        mutate(across(all_of(c("start", "end")),  as.numeric)) %>%
        distinct(gene, roi_name, .keep_all = T)

all_peps <- secretome_roi %>%
                mutate(final_max = max(pep_xgb4c_max, pep_nn4c_max, na.rm = TRUE)) %>%
                group_by(gene) %>%
                slice_max(final_max, n = 1, with_ties = FALSE) %>%
                ungroup() %>%
                select(gene, roi_name, pep_xgb4c_max, pep_nn4c_max, NTC_tot_start, CTC_tot_end) %>%
                mutate(selection = "best")




peps_tp <- peps %>%
  mutate(feature = "new_pep") %>%
  nest(.by = "gene")

peps_tp2 <- all_peps %>%
              mutate(feature = "new_pep") %>%
              nest(.by = "gene") %>%
              filter(!gene %in% peps_tp$gene)

peps_tp <- bind_rows(peps_tp, peps_tp2)




