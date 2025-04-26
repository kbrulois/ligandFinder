



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

set_db_path("/home/groups/ebutcher/kevin/ligandFinder")
#demo("afpd_rename_files", package = "ligandFinder")


run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"

run_id <- "test2"
alg <- "AF2v3"

input_path_models <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models", "/", run_id)

setwd("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts")

analysis_name <- run_id
out_file_name <- paste0(analysis_name, "_contacts")

dir.create(analysis_name)
setwd(analysis_name)
dir.create("input")


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

num_of_grps <- 16
comp_jobs_sub <- comp_jobs_sub %>%
  mutate(group = paste0("job", ntile(n = num_of_grps), ".rds"))



comp_jobs_sub %>%
  group_by(group) %>%
  group_walk(~ saveRDS(.x, file = paste0("input/", .y$group)))

jobs <- unique(comp_jobs_sub$group)

future::plan(strategy = future::multicore(workers = 16))

start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(jobs, \(job) {

  to_do <- readRDS(file = paste0("input/", job))

  lapply(1:nrow(to_do), \(sub_job) {

    dir_name <- to_do %>% slice(sub_job) %>% pull(og_file_name)

    tryCatch({

    metrics <- import_raw_metrics(dir_name = dir_name)

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
                         filter(seq_match == "different"))

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

    modify_file_names(input_path_models = input_path_models,
                      dir_name = dir_name,
                      run_name = run_id,
                      algorithm = "AF2v3")

    }

    metrics <- metrics %>% select(!where(is.list))

    data.table::fwrite(metrics,
                       file = paste0(input_path_models, "/", dir_name, "/", "metrics_v1.csv"))


  }, error = function(e) message("problem with ", job, " ", file_name))

})

  to_do <- to_do %>%
    mutate(metrics = file.exists(paste0(input_path_models, "/", og_file_name, "/metrics_v1.rds")))


  res <- bind_rows(
    map(to_do[["og_file_name"]][to_do[["metrics"]]],
           ~readRDS(paste0(input_path_models, "/", ., "/metrics_v1.rds")))
  )

  saveRDS(res, job)

  message("completed ", job)

})

end <- Sys.time()

end - start

yo()

comp_jobs_sub <- comp_jobs_sub %>%
  mutate(metrics = file.exists(paste0(input_path_models, "/", og_file_name, "/metrics_v1.rds")))

res <- bind_rows(
  map(comp_jobs_sub[["og_file_name"]][comp_jobs_sub[["metrics"]]],
      ~readRDS(paste0(input_path_models, "/", ., "/metrics_v1.rds")))
)

file_name <- paste0(run_analysis_dir, "/", run_id, ".rds")

if(!file.exists(file_name)) saveRDS(res, file_name)



