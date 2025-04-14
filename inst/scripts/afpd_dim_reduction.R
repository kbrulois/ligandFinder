
run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"

out_file_name <- "CXCL14vGPCRs"
reindex <- TRUE

input_files <- list.files(run_analysis_dir)

input_files <- paste0(c("jh_w", "jh_wo", "mm_w", "mm_wo"), ".rds")

input_files <- paste0("CXCL14vGPCRs", ".rds")

res <- bind_rows(
  map(input_files,
      ~readRDS(paste0(run_analysis_dir, "/", .)))
)


known_pairs <- list(c("GPR25", "CXL17"),
                    c("CCR9", "CCL25"),
                    c("GPR15", "GP15L"),
                    c("CML1", "RARR2"),
                    c("CML2", "RARR2"),
                    c("CCRL2", "RARR2"))

ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))

known_pairs2 <- ligand_list %>%
                    rowwise %>%
                    mutate(known = map2(.x = uniprot_name,
                         .y = receptor,
                         .f = \(x, y) {c(x, y)})) %>%
                    ungroup %>%
                    mutate(known_lgl = map_lgl(known, ~any(is.na(.)))) %>%
                    filter(!known_lgl) %>%
                    pull(known) %>%
                    unname

known_pairs <- c(known_pairs, known_pairs2)




res <- res %>%
  rowwise %>%
  mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_id, p2_id) %in% x) == 2)) ~ "known",
                                TRUE ~ "unknown"), .after = "model") %>%
  ungroup


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



for(col_type in names(col_types)) {

  subsetter <- res %>%
    select(ligand_location,
           matches(col_types[[col_type]])) %>%
    select(!where(is.list)) %>%
    rowwise() %>%
    mutate(subsetter = anyNA(c_across(-ligand_location))) %>%
    mutate(subsetter2 = ligand_location == "IC") %>%
    mutate(subsetter = subsetter | subsetter2) %>%
    pull(subsetter)

  dim_red_input <- res %>%
    select(matches(col_types[[col_type]])) %>%
    select(!where(is.list)) %>%
    ungroup %>%
    dplyr::filter(!subsetter) %>%
    as.matrix


  pca_res[[col_type]] <- prcomp(t(dim_red_input), rank. = 10)

  for(i in 1:ncol(pca_res[[col_type]][["rotation"]])) {

    dim_red_name <- paste0("PC", i, "_", col_type)

    res[[dim_red_name]] <- NA

    res[[dim_red_name]][!subsetter] <- pca_res[[col_type]][["rotation"]][,i]

  }



  nmf_res[[col_type]] <- NMFN::nnmf(x = dim_red_input, k = 6)

  for(i in 1:ncol(nmf_res[[col_type]][["W"]])) {

    dim_red_name <- paste0("NMF", i, "_", col_type)

    res[[dim_red_name]] <- NA

    res[[dim_red_name]][!subsetter] <- nmf_res[[col_type]][["W"]][,i]

  }


  umap_res[[col_type]] <- umap::umap(d = nmf_res[[col_type]][["W"]], config = umap_config)

  for(i in 1:ncol(umap_res[[col_type]][["layout"]])) {

    dim_red_name <- paste0("nmfUMAP", i, "_", col_type)

    res[[dim_red_name]] <- NA

    res[[dim_red_name]][!subsetter] <- umap_res[[col_type]][["layout"]][,i]

  }



}




gpcr_cols <- c("p1_id",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily",
               "bw_avail")

res <- res %>%
  relocate(starts_with("PC"), .after = "model") %>%
  relocate(starts_with("UMAP"), .after = "model") %>%
  relocate(starts_with("NMF"), .after = "model") %>%
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_id")} %>%
  relocate(all_of(gpcr_cols[-1]), .after = "model")






saveRDS(pca_res, paste0(out_file_name, "_pca.rds"))
saveRDS(umap_res, paste0(out_file_name, "_umap.rds"))
saveRDS(nmf_res, paste0(out_file_name, "_nmf.rds"))

if(reindex) {
res <- res %>%
  mutate(p2_index_wo_sp = as.character(as.numeric(stringr::str_remove(p2_range, "35-")) - 34), .after = p2_range)
}

res <- res %>%
  dplyr::rename(model_full = model) %>%
  mutate(model = stringr::str_extract(model_full, "model_\\d"), .after = "model_full") %>%
  mutate(pred = stringr::str_extract(model_full, "pred_\\d"), .after = "model")

data.table::fwrite(res %>%
                     select(!where(is.list)),
                   paste0(out_file_name, ".csv"))

message("scp kbrulois@dtn.sherlock.stanford.edu:", getwd(), "/", out_file_name, ".csv", " ~/Desktop/", out_file_name, ".csv")

saveRDS(res, paste0(out_file_name, "_res.rds"))



