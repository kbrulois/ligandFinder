



out_file_name <- "bm_update2"
run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"
reindex <- FALSE

num_of_grps <- 16

future::plan(strategy = future::multicore(workers = num_of_grps))

run_dirs = c("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking",
             "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking_APACE")

res <- bind_rows(map(run_dirs, ~get_metrics(.)))



known_paris <- get_known_pairs()

res <- res %>%
  rowwise %>%
  mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_name, p2_name) %in% x) == 2)) ~ "known",
                                TRUE ~ "unknown"), .after = "iptm") %>%
  ungroup %>%
  mutate(known_pair2 = if_else(known_pair == "known", paste0(p1_name, ";", p2_name), NA))


col_types <- list(all_bw = bw_align[["name"]],
                  cp_bw = grep("_CP$", bw_align[["name"]], value = TRUE),
                  orient = c("EC_lig\\d_mid", "IC_lig\\d_mid", "mid_lig\\d_CT", "mid_lig\\d_NT"),
                  af_qc = c("^pLDDT", "^pae", "iptm", "iptm+ptm"),
                  sum_contacts = c("^ligContacts", "totalCP"),
                  sumcon_cp_bw = c(grep("_CP$", bw_align[["name"]], value = TRUE), c("^ligContacts", "totalCP")))

pca_res <- list()
nmf_res <- list()
umap_res <- list()

umap_config <- umap::umap.defaults
umap_config$min_dist <- 0.5
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200



for(col_type in names(col_types)[4:6]) {

  subsetter <- res %>%
    select(lig1_location,
           matches(col_types[[col_type]])) %>%
    select(!where(is.list)) %>%
    rowwise() %>%
    mutate(subsetter = anyNA(c_across(!any_of(c("lig1_location", "pae_files"))))) %>%
    mutate(subsetter2 = lig1_location == "I") %>%
    mutate(subsetter = subsetter | subsetter2) %>%
    pull(subsetter)

  dim_red_input <- res %>%
    select(-pae_files) %>%
    select(matches(col_types[[col_type]])) %>%
    select(!where(is.list)) %>%
    ungroup %>%
    dplyr::filter(!subsetter) %>%
    as.matrix

  message("dim red columns ", colnames(dim_red_input))

  message("doing PCA")

  pca_res[[col_type]] <- prcomp(t(dim_red_input), rank. = 10)

  for(i in 1:ncol(pca_res[[col_type]][["rotation"]])) {

    dim_red_name <- paste0("PC", i, "_", col_type)

    res[[dim_red_name]] <- NA

    res[[dim_red_name]][!subsetter] <- pca_res[[col_type]][["rotation"]][,i]

  }

  message("doing NMF")


  nmf_res[[col_type]] <- NMFN::nnmf(x = dim_red_input, k = 6)

  for(i in 1:ncol(nmf_res[[col_type]][["W"]])) {

    dim_red_name <- paste0("NMF", i, "_", col_type)

    res[[dim_red_name]] <- NA

    res[[dim_red_name]][!subsetter] <- nmf_res[[col_type]][["W"]][,i]

  }

  message("doing UMAP")

  umap_res[[col_type]] <- umap::umap(d = nmf_res[[col_type]][["W"]], config = umap_config)

  for(i in 1:ncol(umap_res[[col_type]][["layout"]])) {

    dim_red_name <- paste0("nmfUMAP", i, "_", col_type)

    res[[dim_red_name]] <- NA

    res[[dim_red_name]][!subsetter] <- umap_res[[col_type]][["layout"]][,i]

  }



}




gpcr_cols <- c("p1_name",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily")

res <- res %>%
  relocate(starts_with("PC"), .after = "iptm") %>%
  relocate(starts_with("UMAP"), .after = "iptm") %>%
  relocate(starts_with("NMF"), .after = "iptm") %>%
  {left_join(., gpcr_list %>% rename(p1_name = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_name")} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "iptm")

res <- res %>%
  select(-all_of(gpcr_cols[-1])) %>%
  {left_join(., gpcr_list %>% rename(p1_name = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_name")} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "iptm")



saveRDS(pca_res, paste0(out_file_name, "_pca.rds"))
saveRDS(umap_res, paste0(out_file_name, "_umap.rds"))
saveRDS(nmf_res, paste0(out_file_name, "_nmf.rds"))

if(reindex) {
res <- res %>%
  mutate(p2_index_wo_sp = as.character(as.numeric(stringr::str_remove(p2_range, "35-")) - 34), .after = p2_range)
}

res <- res %>%
  dplyr::rename(model_full = model_e) %>%
  mutate(model = stringr::str_extract(model_full, "model_\\d"), .after = "model_full") %>%
  mutate(pred = stringr::str_extract(model_full, "pred_\\d"), .after = "model")


spoc_input <- list.files(pattern = "^spoc_.*.csv$")

spoc <- bind_rows(
          map(spoc_input, ~data.table::fread(.) %>% as_tibble)
)

if(sum(spoc$complex_name %in% res$afpd_dir_name) == nrow(spoc)) {

  res <- left_join(res, spoc %>% rename(afpd_dir_name = complex_name), by = "afpd_dir_name")

}

res <- res %>%
  relocate(any_of(colnames(spoc)), .after = "run_name")

res <- res %>%
  mutate(lig1_NorC = stringr::str_extract(lig1_end, "[NC]$")) %>%
  mutate(lig1_NorC_position = stringr::str_remove(lig1_end, "[NC]$"))

to_join <- to_run %>%
  select(model_trunc, trunc_term, trunc_dir, trunc_size) %>%
  rename(model_name = model_trunc) %>%
  mutate(model_name = stringr::str_replace(model_name, ",\\d+-\\d+;", ";"))

sum(to_join$model_name %in% res$model_name)


res <- left_join(res, to_join, "model_name")




res <- res %>%
  mutate(known_type = case_when(run_name == "bm" & !is.na(known_pair2) ~ "WT",
                                run_name == "bm_truncs" & !is.na(known_pair2) ~ paste0(trunc_term, "_", trunc_dir),
                                TRUE ~ NA))

data.table::fwrite(res %>%
                     select(!where(is.list)),
                   paste0(out_file_name, "_w_spoc_truncs.csv"))

message("scp kbrulois@dtn.sherlock.stanford.edu:", getwd(), "/", out_file_name, "_w_spoc_truncs.csv", " ~/Desktop/", out_file_name, ".csv")

saveRDS(res, paste0(out_file_name, "_res.rds"))



