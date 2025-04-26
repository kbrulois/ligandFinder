



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

set_db_path("/home/groups/ebutcher/kevin/ligandFinder")
#demo("afpd_rename_files", package = "ligandFinder")


run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"

run_id <- "deepX14_2"
alg <- "AF2v3"

input_path_models <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models", "/", run_id)

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))


comp_jobs <- parse_dirname(run_dir = input_path_models,
                           delim_proteins = "_",
                           delim_ranges = "x",
                           delim_start_end = "x") %>%
              mutate(parsed_pair = map(parsed_pair, ~pivot_wider(., names_from=c("protein", "annotation"), values_from=value))) %>%
              unnest(parsed_pair)

comp_jobs_sub <- comp_jobs %>%
  mutate(num_files = map_int(afpd_dir_name, ~length(list.files(paste0(input_path_models, "/", .)))))



jobs <- tibble(dir_name = list.files(input_path_models))

num_of_grps <- 16

jobs <- jobs %>%
  mutate(group = paste0("job", ntile(n = num_of_grps)))


future::plan(strategy = future::multicore(workers = 16))

start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(unique(jobs[["group"]]), \(job) {

  to_do <- jobs %>%
              filter(group == job)

  dirs <- to_do %>% pull(dir_name)

  lapply(dirs[3], \(directory) {

    tryCatch({

    metrics <- import_raw_metrics(dir_name = directory,
                                  run_name = run_id,
                                  algorithm = alg)

    message("number rows after import ", nrow(metrics))

    message("run id ", run_id)

    message("alg ", alg)

    metrics <- left_join(metrics,
                         gpcr_list %>%
                               rename(p1_name = uniprot_name) %>%
                               select(p1_name, `bw: full_table`),
                         by = "p1_name")

    message("number rows of gpcr_list ", nrow(gpcr_list))

    metrics <- process_metrics(input_data = metrics)

    metrics <- bind_rows(
                        compute_RLdists(input_data = metrics %>%
                                    filter(seq_match == "match")),
                        metrics %>%
                         filter(seq_match == "different"))

    message("lig location computed ", "lig1_location" %in% colnames(metrics))

    if("lig1_location" %in% colnames(metrics)) {

    metrics <- bind_rows(
                  metrics %>%
                    filter(seq_match == "match" & lig1_location != "I") %>%
                    mutate(pw_dist = map(pdb, bio3d::dm.pdb)) %>%
                    mutate(contacts = pmap(list(pw_dist = pw_dist,
                                  bw = `bw: full_table`,
                                  pdb.xyz = pdb.xyz,
                                  pae = pae), get_contacts)) %>%
                    unnest(contacts),
                    metrics %>%
                           filter(seq_match == "different" | lig1_location == "I")
                  )

    message("number rows before renaming ", nrow(metrics))

    message("directory: ", directory)

    modify_file_names(input_path_models = input_path_models,
                      dir_name = directory,
                      run_name = run_id,
                      algorithm = alg,
                      metrics = metrics)

    message("renameing done")

    }

    metrics <- metrics %>% select(!where(is.list))

    message("number rows of metrics before saving ", nrow(metrics))


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








jobs <- jobs %>%
  mutate(metrics = file.exists(paste0(input_path_models, "/", dir_name, "/metrics_v1.csv")))

res <- bind_rows(
  map(comp_jobs_sub[["og_file_name"]][comp_jobs_sub[["metrics"]]],
      ~readRDS(paste0(input_path_models, "/", ., "/metrics_v1.rds")))
)

file_name <- paste0(run_analysis_dir, "/", run_id, ".rds")

if(!file.exists(file_name)) saveRDS(res, file_name)



