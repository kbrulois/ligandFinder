

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

fs::dir_copy("/home/groups/ebutcher/kevin/ligandFinder/residue_db",
             "/scratch/groups/ebutcher/deorphan/ligandFinder",
             overwrite = TRUE)

fs::dir_copy("/home/groups/ebutcher/programs/voronota",
             "/scratch/groups/ebutcher/deorphan/ligandFinder",
             overwrite = TRUE)



set_db_path("/scratch/groups/ebutcher/deorphan/ligandFinder")
pq_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db"
voronota_path <- "/home/groups/ebutcher/programs/voronota/bin/voronota-contacts"


gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
gpcr_sub <- gpcr_list %>%
  filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))

gpcr_cols <- c("p1_id",
               "model")

bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

oak_models <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models"
scratch_models <- "/scratch/groups/ebutcher/deorphan/models"


alg <- "AF2v3"
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

#num_cores <- 60
future::plan(strategy = future::multicore(workers = num_cores  %/% 2))


run_dirs <- list.files(scratch_models)

run_dirs <- "known_pairs"

run_dirs <- c("bm", "add_bm", "bm_more_rec", "top200NC")



###extract from OAK

if(FALSE) {

start <- Sys.time()

fs::dir_create(fs::path(scratch_models, paste0(run_dirs)), mode = "u=rwx,g=rwx")

tar_files <- Map(\(x) fs::dir_ls(fs::path(oak_models, x), glob = "*.tar"), run_dirs)

for(run_dir in run_dirs) {

furrr::future_walk2(tar_files[[run_dir]], run_dir, \(f, rd) {

  outdir <- fs::path(scratch_models, rd, stringr::str_remove(fs::path_file(f), ".tar$"))

  tryCatch({
    utils::untar(f, exdir = outdir, extra = "--strip-components=1")
  }, error = function(e) {
    message("❌ Failed to extract ", f, ": ", conditionMessage(e))
  })
})

}

end <- Sys.time()

end - start

yo()

}

####import runs


tmp <- map(run_dirs, ~fs::dir_ls(fs::path(scratch_models, .))) %>% do.call(c, .)

runs <- tibble(afpd_dir_name = fs::path_file(tmp),
               afpd_dir = tmp,
               run_dir = fs::path_dir(tmp)) %>%
  mutate(run_name = fs::path_file(run_dir))


rm(tmp)

runs <- runs %>%
  mutate(file_name_type = case_when(stringr::str_detect(afpd_dir_name, "\\w+_and_\\w+_\\d+-\\d+") ~ "raw_afpd",
                                    stringr::str_detect(afpd_dir_name, "h\\w+_\\w+x\\d+x\\d+") ~ "renamed_dir",
                                    TRUE ~ "unknown")) %>%
  split(., f = .[["file_name_type"]])

sapply(runs, nrow)



runs <- afpd_check_metrics(runs[["renamed_dir"]])

table(runs[["metrics_good"]])

runs <- afpd_check_contacts(runs)

table(runs[["contacts_good"]])


####rename some....

if(FALSE) {

purrr::walk(unique(runs[["raw_afpd"]][["run_dir"]]),

~do_renaming(run_dir = .,
            run_name = fs::path_file(.),
            pairing_dir = NULL,
            afpd_raw = TRUE,
            delim_proteins = "_",
            delim_ranges = "x",
            delim_start_end = "x",
            p1_prefix = "h",
            p1_suffix = NA,
            p2_prefix = "h",
            exclude_p1_range = TRUE,
            site = "SU",
            submitter = "KB",
            algorithm = "AF2v3",
            random_seed = 42)
)


}

###check file name codes
if(FALSE) {


runs <- runs %>%
          mutate(codes = map(file_names, ~stringr::str_extract(., "_[A-Z]{4}[a-z0-9]{7}_") %>%
                                            stringr::str_remove(., "^_[A-Z]{4}") %>%
                                            stringr::str_remove(., "_$") %>%
                                            .[!is.na(.)] %>%
                                            unique))


runs <- runs %>%
          select(-file_names) %>%
          unnest(codes)

codes <- data.table::fread("/home/groups/ebutcher/kevin/ligandFinder/random_codes.csv")

codes <- codes %>% as_tibble()

used_ids <- codes %>%
  dplyr::filter(usage == "used") %>%
  dplyr::pull(id) %>%
  stringr::str_remove("^[A-Z]{4}")

job_codes <- runs$codes %>%
  stringr::str_remove(".{2}$")

sum(used_ids %in% job_codes)

#fs::file_copy("/scratch/groups/ebutcher/deorphan/ligandFinder/random_codes.csv", "/home/groups/ebutcher/kevin/ligandFinder", overwrite = TRUE)


}

####tar existing

if(FALSE) {

purrr::walk(run_dirs, tar_run_dir)

}

####do metric extraction

runs <- runs %>%
  filter(has_json_debug & !contacts_good) %>%
  filter(num_xtr == 5 | num_ark == 5)


runs <- runs %>%
  mutate(group = paste0("job", ntile(n = num_cores)))

runs <- runs %>%
  mutate(run_id = fs::path_file(run_dir)) %>%
  select(p1_name, p2_name, afpd_dir, group, run_dir, run_id)

gc()

future::plan(strategy = future::multicore(workers = num_cores  %/% 2))

start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(unique(runs[["group"]]), \(job) {

  to_do <- runs %>%
    filter(group == job)

  dirs2 <- to_do %>% pull(afpd_dir)

  proteins <- c(to_do[["p1_name"]], to_do[["p2_name"]]) %>% unique

  res_db <- arrow::open_dataset(source = pq_path)

  residue_data <- res_db %>%
    filter(uni_gene %in% proteins) %>%
    collect()

  tryCatch({
    purrr::pwalk(
      list(directory = to_do[["afpd_dir"]],
           job = to_do[["group"]],
           run_name = to_do[["run_id"]]),
      \(directory, job, run_name) {
        do_metrics(directory = directory,
                   job = job,
                   res_dat = residue_data,
                   run_name = run_name)
      }
    )
  }, error = function(e) conditionMessage(e))

  message("completed ", job)

})

end <- Sys.time()

end - start

yo()
