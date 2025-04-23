



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
if (!requireNamespace("arrow", quietly = TRUE)) install.packages("arrow")
install.packages("bio3d", dependencies=TRUE)
install.packages("httr", dependencies=TRUE)

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")

library(ligandFinder)
set_db_path("/home/groups/ebutcher/kevin/ligandFinder")
demo("afpd_rename_files", package = "ligandFinder")

paste0(get_db_path(), "/residue_db")


run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"

run_id <- "MC4R_CART"

input_path_models <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models", "/", run_id)

setwd("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts")

analysis_name <- run_id
out_file_name <- paste0(analysis_name, "_contacts")

dir.create(analysis_name)
setwd(analysis_name)
dir.create("input")

dists_to_comp <- tibble(receptor = c("EC", "IC", "mid", "mid"),
                        ligand = c("mid", "mid", "CT", "NT"))

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))



comp_jobs <- parse_dirname(run_dir = input_path_models) %>%
             make_new_dirname(input = .) %>%
             rename_dir(run_dir = input_path_models,
                        input = .,
                        from = "afpd_dir_name",
                        to = "new_dir_name") %>%
              parse_afpd_files(input = .,
                               dir_name = "new_dir_name",
                               run_dir = input_path_models) %>%
               make_new_file_names(input = .,
                                   dir_name = "new_dir_name",
                                   run_name = run_id)


comp_jobs <- parse_afpd_files(input = tibble(new_dir_name = list.files(input_path_models)),
                 dir_name = "new_dir_name",
                 run_dir = input_path_models) %>%
  make_new_file_names(input = .,
                      dir_name = "new_dir_name",
                      run_name = NA)

comp_jobs %>%
  rename_files(run_dir = input_path_models,
               input = .,
               from = "og_file_name",
               to = "new_file_name")


comp_jobs_sub <- comp_jobs %>%
  mutate(num_files = map_int(og_file_name, ~length(list.files(paste0(input_path_models, "/", .))))) %>%
  filter(num_files == get_mode(num_files)) %>%
  filter(p1_id %in% (gpcr_list %>% filter(bw_avail == "available") %>% pull(uniprot_name)))

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
                               rename(p1_id = uniprot_name) %>%
                               select(p1_id, `bw: full_table`),
                             by = "p1_id")

    metrics <- process_metrics(input_data = metrics)

    metrics <- bind_rows(
                        compute_RLdists(input_data = metrics %>%
                                    filter(seq_match == "match")),
                        metrics %>%
                         filter(seq_match == "different"))

    metrics <- metrics %>%
                  mutate(pw_dist = map(pdb, bio3d::dm.pdb))



    metrics <- bind_rows(
                  metrics %>%
                    filter(seq_match == "match") %>%
                     mutate(contacts = pmap(list(pw_dist = pw_dist,
                                  bw = `bw: full_table`,
                                  pdb.xyz = pdb.xyz,
                                  pae = pae), get_contacts)) %>%
                                  unnest(contacts),
                    metrics %>%
                           filter(seq_match == "different"))


    saveRDS(metrics %>% select(!where(is.list)),
            file = paste0(input_path_models, "/", dir_name, "/", "metrics_v1.rds"))


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


res <- res %>%
  mutate(ligand_location = case_when(EC_lig1_mid < IC_lig1_mid ~ "EC",
                                     EC_lig1_mid > IC_lig1_mid ~ "IC",
                                     TRUE ~ "bw_not_available"), .after = "model") %>%
  mutate(totalCP = replace_na(totalCP, 0)) %>%
  mutate(run_name = run_id, .after = "pdb_files")

file_name <- paste0(run_analysis_dir, "/", run_id, ".rds")

if(!file.exists(file_name)) saveRDS(res, file_name)



