



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



gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

dirs <- fs::dir_ls(input_path_models) %>% stringr::str_subset(., ".tar$")

if(length(dirs) > 0) {
  lapply(dirs, \(x) untar(tarfile = x, exdir = path.expand("~/peptide_alg/testing_set/")))
}


comp_jobs <- parse_dirname(run_dir = input_path_models,
                           delim_proteins = "_",
                           delim_ranges = "x",
                           delim_start_end = "x") %>%
              mutate(parsed_pair = map(parsed_pair, ~pivot_wider(., names_from=c("protein", "annotation"), values_from=value))) %>%
              unnest(parsed_pair)

comp_jobs_sub <- comp_jobs %>%
  mutate(num_files = map_int(afpd_dir_name, ~length(list.files(paste0(input_path_models, "/", .)))))



jobs <- tibble(dir_name = list.files(input_path_models)) %>%
  mutate(complete = map_lgl(dir_name, \(x) {
    file.exists(paste(input_path_models, x, "ranking_debug.json", sep = "/"))
  }))

future::plan(strategy = future::multicore(workers = 16))

jobs <- jobs %>%
          mutate(complete2 = furrr::future_map_lgl(dir_name, \(x) {
            file.exists(paste(input_path_models, x, "metrics_v1.csv", sep = "/"))
          }))

jobs <- jobs %>%
          filter(complete & !complete2)

num_of_grps <- 16

jobs <- jobs %>%
  mutate(group = paste0("job", ntile(n = num_of_grps)))


future::plan(strategy = future::multicore(workers = num_of_grps))

start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(unique(jobs[["group"]]), \(job) {

  to_do <- jobs %>%
              filter(group == job)

  dirs <- to_do %>% pull(dir_name)

  lapply(dirs, \(directory) {

    tryCatch({

    metrics <- import_raw_metrics(dir_name = directory,
                                  run_name = run_id,
                                  algorithm = alg)

    metrics <- left_join(metrics,
                         gpcr_list %>%
                               rename(p1_name = uniprot_name) %>%
                               select(p1_name, `bw: full_table`),
                         by = "p1_name")
    metrics <- process_metrics(input_data = metrics)

    metrics <- bind_rows(
                        compute_RLdists(input_data = metrics %>%
                                    filter(seq_match == "match")),
                        metrics %>%
                         filter(seq_match2 == "different"))

    if("lig1_location" %in% colnames(metrics)) {

    metrics <- bind_rows(
                  metrics %>%
                    filter(seq_match2 == "match" & lig1_location != "I") %>%
                    mutate(pw_dist = map(pdb, bio3d::dm.pdb)) %>%
                    mutate(contacts = pmap(list(pw_dist = pw_dist,
                                  bw = `bw: full_table`,
                                  pdb.xyz = pdb.xyz,
                                  pae = pae), get_contacts)) %>%
                    unnest(contacts),

                    metrics %>%
                           filter(seq_match2 == "different" | lig1_location == "I")
                  )

    modify_file_names(input_path_models = input_path_models,
                      dir_name = directory,
                      run_name = run_id,
                      algorithm = alg,
                      metrics = metrics)

    }

    metrics <- metrics %>% select(!where(is.list))

    data.table::fwrite(metrics,
                       file = paste(input_path_models, directory, "metrics_v1.csv", sep = "/"))

    message("saved to ", paste(input_path_models, directory, "metrics_v1.csv", sep = "/"))


  }, error = function(e) message("problem with ", job, " ", directory))

})

  to_do <- to_do %>%
    mutate(metrics = file.exists(paste0(input_path_models, "/", dir_name, "/metrics_v1.csv")))


  res <- bind_rows(
    map(to_do[["dir_name"]][to_do[["metrics"]]],
        ~data.table::fread(paste0(input_path_models, "/", ., "/metrics_v1.csv")) %>% as_tibble)
  )

  saveRDS(res, paste0(job, ".rds"))

  message("completed ", job)

})

end <- Sys.time()

end - start

yo()






