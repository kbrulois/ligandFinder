




.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

set_db_path("/scratch/groups/ebutcher/deorphan/ligandFinder")
pq_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db"
voronota_path <- "/scratch/groups/ebutcher/kevin/voronota/bin/voronota-contacts"

res_db <- arrow::open_dataset(source = pq_path)
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
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) %/% 1.5

#num_cores <- 60
future::plan(strategy = future::multicore(workers = num_cores))


run_id <- "cxc17_gp15l"

run_dirs <- c("CXCL14_jh_w", "CXCL14_jh_wo", "CXCL14_mm_w", "CXCL14_mm_wo")
run_dirs <- "bm_sep28"
run_dirs <- "GPCRvCXCL14_oct7"
run_dirs <- c("top200NC_Nov12", "top200NC_Oct23_cleanup")
run_dirs <- "brinp"
run_dirs <- "neuro"

###extract from OAK

start <- Sys.time()

time <- format(Sys.time(), "%b%e")

run_dirs <- c("top200NC")

fs::dir_create(fs::path(scratch_models, paste0(run_dirs, "_", time)), mode = "u=rwx,g=rwx")

tar_files <- Map(\(x) fs::dir_ls(fs::path(oak_models, x), glob = "*.tar"), run_dirs)

for(run_dir in run_dirs) {

furrr::future_walk2(tar_files[[run_dir]], run_dir, \(f, rd) {

  outdir <- fs::path(scratch_models, paste0(rd, "_", time), stringr::str_remove(fs::path_file(f), ".tar$"))

  tryCatch({
    utils::untar(f, exdir = outdir)
  }, error = function(e) {
    message("❌ Failed to extract ", f, ": ", conditionMessage(e))
  })
})

}

end <- Sys.time()

start - end

yo()



####Check run dirs

if(exists("time")) {
  run_dirs <- stringr::str_c(run_dirs, "_", time)
}


tmp <- furrr::future_map(run_dirs, ~fs::dir_ls(fs::path(scratch_models, .))) %>% do.call(c, .)

jobs <- tibble(afpd_dir_name = fs::path_file(tmp),
               afpd_dir = tmp,
               run_dir = fs::path_dir(tmp))

rm(tmp)

jobs <- jobs %>%
  mutate(file_name_type = case_when(stringr::str_detect(afpd_dir_name, "\\w+_and_\\w+_\\d+-\\d+") ~ "raw_afpd",
                                    stringr::str_detect(afpd_dir_name, "h\\w+_\\w+x\\d+x\\d+") ~ "renamed_dir",
                                    TRUE ~ "unknown")) %>%
  split(., f = .[["file_name_type"]])

sapply(jobs, nrow)



####rename some....

if(FALSE) {

if("renamed_dir" %in% names(jobs)) {
  jobs[["renamed_dir"]] <- jobs[["renamed_dir"]] %>%
    mutate(parse_proteins(file = afpd_dir_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
    #mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
    {left_join(., gpcr_sub %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
               by = "p1_id")} %>%
    mutate(model = paste0(model,
                          ";",
                          p2_id,
                          ",",
                          stringr::str_replace(p2_range, "x", "-")))
}

if("raw_afpd" %in% names(jobs)) {
  jobs[["raw_afpd"]] <- jobs[["raw_afpd"]] %>%
    mutate(parse_proteins(file = afpd_dir_name)) %>%
    mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
    {left_join(., gpcr_sub %>% rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
               by = "p1_id")} %>%
    mutate(model = paste0(model,
                          ";",
                          p2_id,
                          if_else(p2_range == "", "", paste0(",", p2_range))))


}

to_rename <- jobs[["raw_afpd"]] %>%
              filter(!model %in% !!jobs[["renamed_dir"]][["model"]])

rename_dir <- fs::path(scratch_models, "top200NCnewnew")

fs::dir_create(rename_dir)

to_rename <- to_rename %>%
  mutate(rename_dir = rename_dir)


furrr::future_walk(to_rename[["afpd_dir"]], \(d) {
  fs::dir_copy(d, fs::path(rename_dir, fs::path_file(d)))
})

yo()

test <- fs::dir_ls(rename_dir)

run_dir_rename <- jobs[["raw_afpd"]][["run_dir"]][[1]]

rename_data <- parse_dirname(run_dir = run_dir_rename,
                             afpd_raw = TRUE)

rename_data <- make_new_dirname(input = rename_data,
                                delim_proteins = "_",
                                delim_ranges = "x",
                                delim_start_end = "x",
                                p1_prefix = "h",
                                p1_suffix = NA,
                                p2_prefix = "h",
                                exclude_p1_range = TRUE)

rename_dir(run_dir = run_dir_rename,
           input = rename_data,
           from = "afpd_dir_name",
           to = "new_dir_name")


rename_data <- parse_afpd_files(input = rename_data,
                                dir_name = "new_dir_name",
                                run_dir = run_dir_rename)

rename_data <- make_new_file_names(input = rename_data,
                                   dir_name = "new_dir_name",
                                   run_name = run_id,
                                   site = "SU",
                                   submitter = "KB",
                                   algorithm = "AF2v3",
                                   random_seed = 42)

rename_data <- rename_files(run_dir = run_dir_rename,
                            input = rename_data,
                            from = "og_file_name",
                            to = "new_file_name")







}




jobs <- jobs[["renamed_dir"]]

jobs <- jobs %>%
  mutate(parse_proteins(afpd_dir_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  dplyr::rename(p1_name = p1_id, p2_name = p2_id) %>%
  mutate(file_names = map(afpd_dir, ~fs::dir_ls(.) %>% fs::path_file(.)))




####check curent metrics

jobs <- jobs %>%
  mutate(furrr::future_map_dfr(file_names, \(x) {
        tibble(has_json_debug = "ranking_debug.json" %in% x,
               has_metrics_csv = "metrics_v2.csv" %in% x,
               has_contact_rds = "metrics_v2c.rds" %in% x,
               num_models = sum(grepl("_xtr_.*.pdb$", x)))
    })) %>%
  mutate(metrics = furrr::future_map(afpd_dir, \(x) {
    file <- fs::path(x, "metrics_v2.csv")
    if(file.exists(file)) {
      return(data.table::fread(file) %>% as_tibble)
    } else {
      return(NULL)
    }
  }))

metrics_good <- jobs %>%
  filter(has_metrics_csv) %>%
  filter(!map_lgl(metrics, is.null)) %>%
  filter(map_int(metrics, ~nrow(.)) == num_models) %>%
  filter(map_lgl(metrics, ~"afpd_dir" %in% colnames(.))) %>%
  pull(afpd_dir_name)


jobs <- jobs %>%
  mutate(metrics_good = afpd_dir_name %in% !!metrics_good)

table(jobs[["metrics_good"]], jobs[["run_dir"]])




###check contacts

jobs <- jobs %>%
  mutate(contacts = furrr::future_map(afpd_dir, \(x) {
    file <- fs::path(x, "metrics_v2c.rds")
    if(file.exists(file)) {
      return(tryCatch({readRDS(file)}, error = function(e) NULL))
    } else {
      return(NULL)
    }
  }))



contacts_good <- jobs %>%
  filter(has_contact_rds) %>%
  #filter(!map_lgl(contacts, is.null)) %>%
  filter(map2_lgl(metrics, contacts, \(x, y) {
      if(is.null(y)) {test <- rep("irrelevant", nrow(x))} else {
      test <- map_chr(y, \(z) {
                if(is.null(z)) {return("irrelevant")} else
                if(nrow(z) == 0) {return("irrelevant")} else
                if(!"area" %in% colnames(z)) {return("irrelevant")} else
                {return("relevant")}
      })
      }
      return(identical(x[["location"]], test))
    })) %>%
  pull(afpd_dir_name)

jobs <- jobs %>%
  mutate(contacts_good = afpd_dir_name %in% !!contacts_good)

table(jobs[["contacts_good"]])




jobs <- jobs %>%
  filter(has_json_debug & !contacts_good)


jobs <- jobs %>%
  mutate(group = paste0("job", ntile(n = num_cores)))

jobs <- jobs %>% select(p1_name, p2_name, afpd_dir, group)

gc()

future::plan(strategy = future::multicore(workers = num_cores))


start <- Sys.time()

options(future.globals.maxSize = 10e9)

furrr::future_map(unique(jobs[["group"]]), \(job) {

  to_do <- jobs %>%
    filter(group == job)

  dirs2 <- to_do %>% pull(afpd_dir)

  proteins <- c(to_do[["p1_name"]], to_do[["p2_name"]]) %>% unique

  res_db <- arrow::open_dataset(source = pq_path)

  residue_data <- res_db %>%
    filter(uni_gene %in% proteins) %>%
    collect()

  tryCatch({
    lapply(dirs2, \(x) do_metrics(directory = x, job = job, res_dat = residue_data))
  }, error = function(e) conditionMessage(e))

  message("completed ", job)

})

end <- Sys.time()

end - start

yo()


