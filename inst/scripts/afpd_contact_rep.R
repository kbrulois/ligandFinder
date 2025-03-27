



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)

devtools::install_github("kbrulois/ligandFinder", auth_token = "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")

library(ligandFinder)

input_path_models <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/mxrun"

getwd()
setwd("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts")

analysis_name <- "mxrun_analysis"
out_file_name <- paste0(analysis_name, "_contacts")


dir.create(analysis_name)
setwd(analysis_name)
dir.create("input")

dists_to_comp <- tibble(receptor = c("EC", "IC", "mid", "mid"),
                        ligand = c("mid", "mid", "CT", "NT"))

#gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
gpcr_list <- readRDS("/oak/stanford/groups/ebutcher/deorphan-AI-ze/uniprot/gpcrdb_receptor_list_KB_250304_final.rds")

bw_align <- summarize_bw(gpcr_list = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/uniprot/gpcrdb_receptor_list_KB_250304_final.rds")

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))


comp_jobs <- tibble(og_file_name = list.files(input_path_models))

comp_jobs <- comp_jobs %>%
mutate(parse_proteins(og_file_name,
                      delim_proteins = "_and_",
                      delim_ranges = "_",
                      delim_start_end = "-")) %>%
  mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
  mutate(new_file_name = paste(paste0("h", p1_id, p1_range_type), paste(p2_id, p2_range, sep = "x"), sep = "_"), .after = "og_file_name")

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

    file_name <- to_do %>% slice(sub_job) %>% pull(og_file_name)

    tryCatch({

    metrics <- import_raw_metrics(input_data = to_do %>% slice(sub_job))

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
            file = paste0(input_path_models, "/", file_name, "/", "metrics_v1.rds"))


  }, error = function(e) message("problem with ", job, file_name))

})

  message("completed ", job)

})

end <- Sys.time()

end - start

yo()

nums <- list(1:4, 5:8, 9:12, 13:16)

for(x in seq_along(nums)) {

test <- lapply(paste0("job", x, ".rds"), readRDS)

test2 <- lapply(test, \(x) {x %>% select(!where(is.list))})

saveRDS(bind_rows(test2), paste0("consolidated", x, ".rds"))

}

res2 <- bind_rows(
  lapply(paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/pocForGrant7/consolidated", 0:3, ".rds"), readRDS)
)


res <- bind_rows(
  lapply(jobs, readRDS)
)

res <- bind_rows(res, res2)



