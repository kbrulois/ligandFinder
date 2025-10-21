
runs <- readRDS('runs.rds')
runs <- runs %>%
  mutate(v2_computed = file.exists(paste0(input_path_models, "/", afpd_dir_name, "/metrics_v2.csv"))) %>%
  filter(v2_computed)

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))


out_file_name <- "bm_update5"
run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"
spoc_path <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/bm_final_spoc"


reindex <- FALSE


num_of_grps <- 20

future::plan(strategy = future::multicore(workers = num_of_grps))


run_dirs <- input_path_models <- c("/scratch/groups/ebutcher/deorphan/models/bm_sep28", )

names(run_dirs) <- sapply(run_dirs, basename) %>% unname

get_data <- function(input_path_models) {
tibble(afpd_dir_name = fs::dir_ls(input_path_models) %>% stringr::str_remove(., ".tar$") %>% basename(),
              file_parts = map(afpd_dir_name, ~stringr::str_split(., "_", simplify = TRUE))) %>%
  mutate(file_part_len = map_int(file_parts, length)) %>%
  mutate(parse_proteins(afpd_dir_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  mutate(data_files = furrr::future_map(afpd_dir_name, ~fs::dir_ls(paste0(input_path_models, "/", .)) %>% basename())) %>%
  mutate(num_files = furrr::future_map_int(afpd_dir_name, ~length(list.files(paste0(input_path_models, "/", .))))) %>%
  mutate(has_ranking_debug = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "ranking_debug.json", sep = "/"))
  })) %>%
  mutate(has_metrics_csv = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "metrics_v2.csv", sep = "/"))
  })) %>%
  mutate(has_contact_rds = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "metrics_v2c.rds", sep = "/"))
  }))

}

res <- map(run_dirs, get_data)

res <- map2(names(run_dirs), res, \(x, y) {
  y %>% mutate(run_name = x)
})

res <- bind_rows(res)

res <- res %>%
  filter(has_metrics_csv)

res <- res %>%
        mutate(metrics = furrr::future_map2(afpd_dir_name, run_name, \(x, y) {

          file <- paste(run_dirs[y], x, "metrics_v2.csv", sep = "/")
          if(file.exists(file)) {
            return(data.table::fread(file) %>% as_tibble)
          } else {
            return("none")
          }
        })) %>%
  mutate(contacts = furrr::future_map2(afpd_dir_name, run_name, \(x, y) {

    file <- paste(run_dirs[y], x, "metrics_v2c.rds", sep = "/")
    if(file.exists(file)) {
      return(readRDS(file))
    } else {
      return("none")
    }
  }))

res_cols <- colnames(res)

table(sapply(res$metrics, nrow))

res3 <- res %>%
  filter(has_metrics_csv) %>%
  filter(!map_lgl(metrics, is.character)) %>%
  filter(map_lgl(metrics, ~nrow(.) != 0)) %>%
  #filter(p1_id %in% c("NPY2R", "NPY5R")) %>%
  mutate(metrics = map(metrics, ~ select(.x, !any_of(res_cols)))) %>%
  unnest(metrics)


known_pairs <- get_known_pairs()

res3 <- res3 %>%
  rowwise %>%
  mutate(known_pair = case_when(any(map_lgl(known_pairs, \(x) sum(c(p1_name, p2_name) %in% x) == 2)) ~ "known",
                                TRUE ~ "unknown"), .after = "p1_name") %>%
  ungroup %>%
  mutate(known_pair2 = if_else(known_pair == "known", paste0(p1_name, ";", p2_name), NA))



gpcr_cols <- c("p1_id",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily")

res3 <- res3 %>%
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_id")} %>%
  relocate(all_of(gpcr_cols[-1]), .before = "p1_id")

res3 <- res3 %>%
  select(-all_of(gpcr_cols[-1])) %>%
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_id")} %>%
  relocate(all_of(gpcr_cols[-1]), .before = "p1_id")

res3 <- res3 %>%
  mutate(brinp = if_else(run_name == "bm", known_pair, paste0(p2_id, "_", p2_range))) %>%
  mutate(gpcr = if_else(p1_id %in% c("NPY2R", "NPY5R"), p1_id, "other"))

out_dir <- "/oak/stanford/groups/ebutcher/kevin"
local_dir <- "~/AF2_analysis"
file_name <- "brinp_v7.csv"

data.table::fwrite(res3 %>%
                     select(!where(is.list)),
                   paste0(out_dir, "/", file_name))

message("scp kbrulois@dtn.sherlock.stanford.edu:", out_dir, "/", file_name, " ", local_dir, "/", file_name)


res2 <- res %>%
  mutate(sum_contacts = furrr::future_map(contacts, \(x) {

    logi <- x != "none"
    x[logi] <- lapply(x[logi], summarize_contacts)
    x[!logi] <- NULL
    return(x)})) %>%
  unnest(sum_contacts)

left_join(res3, res2, by = , suffix = c("", "_get_rid_of_it"))













res %>%
  unnest()

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



for(col_type in names(col_types)[c(2,4,5,6)]) {

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

yo()




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


spoc_input <- list.files(path = spoc_path,
                         pattern = "^spoc_.*.csv$")

spoc <- bind_rows(
          map(paste0(spoc_path, "/", spoc_input), ~data.table::fread(.) %>% as_tibble)
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
                   paste0(out_file_name, "_w_spoc.csv"))

message("scp kbrulois@dtn.sherlock.stanford.edu:", getwd(), "/", out_file_name, "_w_spoc.csv", " ~/Desktop/", out_file_name, ".csv")

saveRDS(res, paste0(out_file_name, "_res.rds"))



