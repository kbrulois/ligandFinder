





.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

set_db_path("/home/groups/ebutcher/kevin/ligandFinder")
#demo("afpd_rename_files", package = "ligandFinder")


run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"

run_id <- "bm"
alg <- "AF3"

input_path_models <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models", "/", "benchmarking_APACE")
input_path_models <- "~/peptide_alg/AF3_test"
input_path_models <- "/scratch/groups/ebutcher/deorphan/models/benchmarking"

pq_path <- "~/ligandFinder_data/residue_db"
pq_path <- "/home/groups/ebutcher/kevin/ligandFinder/residue_db"
voronota_path <- "/usr/local/bin/voronota-contacts"
voronota_path <- "/home/groups/ebutcher/programs/voronota/bin/voronota-contacts"


res_db <- arrow::open_dataset(source = pq_path)




gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))



num_of_grps <- 16

future::plan(strategy = future::multicore(workers = num_of_grps))

run_dir <- input_path_models






tmp <- tibble(files = fs::dir_ls(run_dir) %>% basename(),
              file_parts = map(files, ~stringr::str_split(., "_", simplify = TRUE))) %>%
  mutate(file_part_len = map_int(file_parts, length)) %>%
  mutate(parse_proteins(files, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  mutate(data_files = furrr::future_map(files, ~fs::dir_ls(.) %>% basename()))


tmp <- tmp %>%
  mutate(num_files = map_int(data_files, length))


tmp <- tmp %>%
        mutate(log = furrr::future_map(files, \(x) {

          to_import <- paste0(c(input_path_models, x, "file_name_log.csv"), collapse = "/")
          if(file.exists(to_import)) {
          return(data.table::fread(to_import))
          } else {
            return(NULL)
          }

            }))


tmp <- tmp %>%
  mutate(to_filter = map_lgl(log, is.null)) %>%
  filter(!to_filter)

tmp <- tmp %>%
  mutate(mod_file_name_col_exists = map_lgl(log, \(x) {"mod_file_name" %in% colnames(x)}))


tmp <- tmp %>%
  mutate(furrr::future_map2_dfr(log, data_files, \(x, y) {

    if("mod_file_name" %in% colnames(x)) {
      mod_col <- "mod_file_name"
    } else {
      mod_col <- "new_file_name"
    }

    tibble(mod_file_ol = sum(x[["mod_file_name"]] %in% y),
           new_file_ol = sum(x[["new_file_name"]] %in% y),
           directory_files_in_log = list(y[y %in% x[[mod_col]]]),
           directory_files_not_in_log = list(y[!y %in% x[[mod_col]]]),
           log_files_in_directory = list(x[[mod_col]][x[[mod_col]] %in% y]),
           log_files_not_in_directory = list(x[[mod_col]][!x[[mod_col]] %in% y]))
  }))



