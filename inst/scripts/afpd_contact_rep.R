



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
alg <- "AF2v3"

input_path_models <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models", "/", "benchmarking_APACE")
input_path_models <- "~/peptide_alg/testing_set"

pq_path <- "~/ligandFinder_data/residue_db"
voronota_path <- "/usr/local/bin/voronota-contacts"
voronota_path <- "/usr/local/bin/voronota-contacts"


res_db <- arrow::open_dataset(source = pq_path)




gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

dirs <- fs::dir_ls(input_path_models) %>% stringr::str_subset(., ".tar$")

if(length(dirs) > 0) {
  lapply(dirs, \(x) untar(tarfile = x, exdir = path.expand("~/peptide_alg/testing_set/")))
}


runs <- parse_dirname(run_dir = input_path_models,
                      delim_proteins = "_",
                      delim_ranges = "x",
                      delim_start_end = "x") %>%
              mutate(parsed_pair = map(parsed_pair, ~pivot_wider(., names_from=c("protein", "annotation"), values_from=value))) %>%
              unnest(parsed_pair) %>%
              mutate(num_files = map_int(afpd_dir_name, ~length(list.files(paste0(input_path_models, "/", .))))) %>%
              mutate(complete = map_lgl(afpd_dir_name, \(x) {
                  file.exists(paste(input_path_models, x, "ranking_debug.json", sep = "/"))
                    })) %>%
              mutate(complete2 = furrr::future_map_lgl(afpd_dir_name, \(x) {
                file.exists(paste(input_path_models, x, "metrics_v1.csv", sep = "/"))
                      }))

jobs <- runs %>%
          filter(complete)


future::plan(strategy = future::sequential())

num_of_grps <- 16

jobs <- jobs %>%
  mutate(group = paste0("job", ntile(n = num_of_grps)))


future::plan(strategy = future::multicore(workers = num_of_grps))

start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(unique(jobs[["group"]])[1:2], \(job) {

  to_do <- jobs %>%
              filter(group == job)

  dirs2 <- to_do %>% pull(afpd_dir_name)


  proteins <- c(to_do[["p1_name"]], to_do[["p2_name"]]) %>% unique

  residue_data <- res_db %>%
    filter(uni_gene %in% proteins) %>%
    collect()


  lapply(dirs2, \(x) do_metrics(directory = x, job = job, residue_data = residue_data))

  to_do <- to_do %>%
    mutate(metrics = file.exists(paste0(input_path_models, "/", afpd_dir_name, "/metrics_v2.csv")))


  res <- bind_rows(
    map(to_do[["afpd_dir_name"]][to_do[["metrics"]]],
        ~data.table::fread(paste0(input_path_models, "/", ., "/metrics_v2.csv")) %>% as_tibble)
  )

  saveRDS(res, paste0(job, ".rds"))

  message("completed ", job)

})

end <- Sys.time()

end - start

yo()






