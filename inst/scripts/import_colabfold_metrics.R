






.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)

set_db_path("/home/groups/ebutcher/kevin/ligandFinder")
#demo("afpd_rename_files", package = "ligandFinder")


run_analysis_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/run_analyses"

run_id <- "X12"
alg <- "AF2v3"

input_path_models <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models", "/", "benchmarking_APACE")
input_path_models <- "~/peptide_alg/testing_set"
input_path_models <- "/scratch/groups/ebutcher/deorphan/models/cxcl12"

pq_path <- "~/ligandFinder_data/residue_db"
pq_path <- "/home/groups/ebutcher/kevin/ligandFinder/residue_db"
voronota_path <- "/usr/local/bin/voronota-contacts"
voronota_path <- "/home/groups/ebutcher/programs/voronota/bin/voronota-contacts"


res_db <- arrow::open_dataset(source = pq_path)


all_files = fs::dir_ls(input_path_models)

pat <- ".done.txt$"

dat <- tibble(complex_name = all_files %>%
                             fs::path_file(.) %>%
                             stringr::str_subset(., pat) %>%
                             stringr::str_remove(., pat))

dat <- dat %>%
        mutate(files = map(complex_name, ~stringr::str_subset(all_files, paste0("/",.))))

dat <- dat %>%
        mutate(pdb_files = map(files, ~grep(".pdb$", ., value = TRUE)),
               pae_files = map(files, ~grep("predicted_aligned_error_v1.json$", ., value = TRUE)),
               score_files = map(files, ~grep("_scores_.*.json$", ., value = TRUE)))


dat <- dat %>%
        mutate(pdb_files = map(pdb_files,
        ~tibble(pdb_files = .,
                rank = stringr::str_extract(., "_rank_\\d+") %>% stringr::str_remove(., "^_rank_"),
                model_num = stringr::str_extract(., "_model_\\d+") %>% stringr::str_remove(., "^_model_"),
                seed = stringr::str_extract(., "_seed_\\d+") %>% stringr::str_remove(., "^_seed_"),
                model_e = paste0("rank_", rank, "_alphafold2_multimer_v3_", "model_", model_num, "_seed_", seed))
        )) %>%
  mutate(score_files = map(score_files,
                         ~tibble(score_files = .,
                                 rank = stringr::str_extract(., "_rank_\\d+") %>% stringr::str_remove(., "^_rank_"),
                                 model_num = stringr::str_extract(., "_model_\\d+") %>% stringr::str_remove(., "^_model_"),
                                 seed = stringr::str_extract(., "_seed_\\d+") %>% stringr::str_remove(., "^_seed_"),
                                 model_e = paste0("rank_", rank, "_alphafold2_multimer_v3_", "model_", model_num, "_seed_", seed))
  ))

dat %>%
  mutate(is_match = map2_lgl(pdb_files, score_files, \(x, y) {
    identical(x[["model_e"]], y[["model_e"]])
  })) %>%
  filter(!is_match)

dat1 <- dat %>%
  unnest(pdb_files)

dat2 <- dat %>%
  unnest(score_files)

dat <- left_join(dat1 %>% select(-score_files), dat2 %>% select(complex_name, score_files, model_e), by = join_by(complex_name, model_e))

dat <- dat %>%
  mutate(scores = map(score_files, ~jsonlite::read_json(.)))


dat <- dat %>%
          mutate(iptm = map_dbl(scores, \(x) x[["iptm"]]))


test <- data.table::fread("~/AF2_analysis/cxcl12.csv") %>% as_tibble

test2 <- tibble(complex_name = clipr::read_clip())

test <- test[match(test2$complex_name, test$complex_name), ]

clipr::write_clip(test$iptm)


gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

dirs <- fs::dir_ls(input_path_models) %>% stringr::str_subset(., ".tar$")

if(length(dirs) > 0) {
  lapply(dirs, \(x) untar(tarfile = x, exdir = path.expand("~/peptide_alg/testing_set/")))
}

num_of_grps <- 16

future::plan(strategy = future::multicore(workers = num_of_grps))


runs <- parse_dirname(run_dir = input_path_models,
                      delim_proteins = "_",
                      delim_ranges = "x",
                      delim_start_end = "x") %>%
  mutate(parsed_pair = furrr::future_map(parsed_pair, ~pivot_wider(., names_from=c("protein", "annotation"), values_from=value))) %>%
  unnest(parsed_pair) %>%
  mutate(num_files = furrr::future_map_int(afpd_dir_name, ~length(list.files(paste0(input_path_models, "/", .))))) %>%
  mutate(complete = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "ranking_debug.json", sep = "/"))
  })) %>%
  mutate(complete2 = furrr::future_map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "metrics_v1.csv", sep = "/"))
  }))

jobs <- runs %>%
  filter(complete)


jobs <- jobs %>%
  mutate(data_files = furrr::future_map(afpd_dir_name, ~fs::dir_ls(paste0(input_path_models, "/", .)) %>% basename()))

jobs <- jobs %>%
  mutate(num_E_models = map_int(data_files, \(x) {
    grep("_ark_.*_[A-Z]{4}[a-z]{7}_s\\d+m\\d+p\\d+_r\\d+", x, value = TRUE) %>%
      stringr::str_detect(., "_E_") %>%
      sum})) %>%
  mutate(v2c_present = map_chr(data_files, \(x) {
    ifelse("metrics_v2c.rds" %in% x, "yes", "no")
  }))

jobs <- jobs %>%
  mutate(contact_raw = furrr::future_map(afpd_dir_name, \(x) {

    file <- paste(input_path_models, x, "metrics_v2c.rds", sep = "/")
    if(file.exists(file)) {
      return(readRDS(file))
    } else {
      return("none")
    }
  }))

jobs <- jobs %>%
  mutate(contact_good = furrr::future_map_int(contact_raw, \(x) {
    sum(!sapply(x, is.null))
  }))





jobs <- jobs %>%
  mutate(group = paste0("job", ntile(n = num_of_grps)))

future::plan(strategy = future::multicore(workers = num_of_grps))


start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(unique(jobs[["group"]]), \(job) {

  to_do <- jobs %>%
    filter(group == job)

  dirs2 <- to_do %>% pull(afpd_dir_name)


  proteins <- c(to_do[["p1_name"]], to_do[["p2_name"]]) %>% unique

  residue_data <- res_db %>%
    filter(uni_gene %in% proteins) %>%
    collect()

  tryCatch({
    lapply(dirs2, \(x) do_metrics(directory = x, job = job, residue_data = residue_data))
  }, error = function(e) e)

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



runs %>%
  mutate(v2_computed = file.exists(paste0(input_path_models, "/", afpd_dir_name, "/metrics_v2.csv"))) %>%
  pull(v2_computed) %>%
  sum(.)


test2 <- test %>%
  group_by(afpd_dir_name) %>%
  summarize(numE = sum(lig1_location == "E" & location != "APPP"))

test3 <- left_join(test, test2, by = "afpd_dir_name")
