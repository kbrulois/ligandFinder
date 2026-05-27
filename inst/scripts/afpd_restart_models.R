



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
future::plan(strategy = future::multicore(workers = num_cores))

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/top200NC"

out_dir <- c("/scratch/groups/ebutcher/deorphan/models/top200NC_Oct23",
             "/scratch/groups/ebutcher/deorphan/models/top200NC_Oct23_cleanup")


out_dir <- c("/scratch/groups/ebutcher/deorphan/models/top200NCnew",
             "/scratch/groups/ebutcher/deorphan/models/top200NCnewnew")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/brinp"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/brinp_final")


job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/add_bm"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/add_bm_Oct28")


job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/neuro"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/neuro")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/CXCL14peptides"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/CXCL14peptides")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/sven_fig_peps"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/sven_fig_peps")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/new_pep"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/new_peps")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/bm_more_rec"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/bm_more_rec")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/jan30"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/jan30")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/uni_pep"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/uni_pep")

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/top50dbc"
out_dir <- c("/scratch/groups/ebutcher/deorphan/models/top50dbc")

run_name <- "dTbm"

job_dir <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/", run_name)
out_dir <- paste0("/scratch/groups/ebutcher/deorphan/models/", run_name)


gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
gpcr_sub <- gpcr_list %>%
  #filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  #filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  #filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160 | is.na(`bw: length N-term`), model_name_dNT, model_name))

gpcr_cols <- c("p1_id",
               "model")

gpcr_sub


job_dat <- tibble(jobs = list.files(job_dir) %>% stringr::str_subset(., ".txt$"))

job_dat <- job_dat[gtools::mixedorder(job_dat[["jobs"]]), ]

job_dat <- job_dat %>%
  mutate(model = map(jobs, ~data.table::fread(paste0(job_dir, "/", .), sep = NULL, header = FALSE, col.names = "model")))

job_dat <- job_dat %>%
  unnest(model)

test <- Map(list.files, out_dir)

out_dat <- tibble(afpd_dir_name = do.call(c, test),
                  out_dir = do.call(c, lapply(names(test), \(x) rep(x, length(test[[x]])))))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

out_dat <- out_dat %>%
            mutate(file_name_type = case_when(stringr::str_detect(afpd_dir_name, "\\w+_and_\\w+_\\d+-\\d+") ~ "raw_afpd",
                                              stringr::str_detect(afpd_dir_name, "h\\w+_\\w+x\\d+x\\d+") ~ "renamed_dir",
                                              TRUE ~ "unknown")) %>%
            split(., f = .[["file_name_type"]])


sapply(out_dat, nrow)

if("renamed_dir" %in% names(out_dat)) {
out_dat[["renamed_dir"]] <- out_dat[["renamed_dir"]] %>%
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

if("raw_afpd" %in% names(out_dat)) {
  out_dat[["raw_afpd"]] <- out_dat[["raw_afpd"]] %>%
    mutate(parse_proteins(file = afpd_dir_name)) %>%
    mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
    {left_join(., gpcr_sub %>% dplyr::rename(p1_id = uniprot_name) %>% select(all_of(gpcr_cols)),
               by = "p1_id")} %>%
    mutate(model = paste0(model,
                          ";",
                          p2_id,
                          if_else(p2_range == "", "", paste0(",", p2_range))))

}

out_dat <- bind_rows(out_dat)


out_dat <- out_dat %>%
  mutate(has_json_debug = map2_lgl(afpd_dir_name, out_dir, \(x, y) {
    file.exists(paste(y, x, "ranking_debug.json", sep = "/"))
  })) %>%
  {print(table(.[["has_json_debug"]])); .} %>%
  filter(has_json_debug)

sum(job_dat[["model"]] %in% out_dat[["model"]])

job_dat <- job_dat %>%
  mutate(complete = model %in% out_dat[["model"]])


table(job_dat$complete)

out_dat %>%
  filter(!model %in% job_dat$model) %>%
  select(p1_id, p2_id) %>%
  print(n = 1000)

job_dat %>%
  filter(!complete) %>%
  print(n = 1000)


#####rename current models

run_dir <- unique(out_dat$out_dir)

do_renaming(run_dir = run_dir,
            run_name = fs::path_file(run_dir),
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




####tar existing models

run_dir <- unique(out_dat$out_dir)


purrr::walk(run_dirs, tar_run_dir)







####set up clean up run


afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)

job_dat <- job_dat %>%
    mutate(parse_proteins(model, delim_proteins = ";", delim_ranges = ",", delim_start_end = "-")) %>%
    mutate(in_afpd_db = p1_id %in% afpd_db & p2_id %in% afpd_db)

table(job_dat[["in_afpd_db"]])

#job_names <- system("squeue -h -o %j -t RUNNING -u kbrulois", intern = TRUE) %>% stringr::str_subset(., "^job")

to_still_run <- job_dat %>%
  filter(!model %in% out_dat[["model"]]) %>%
  filter(in_afpd_db) %>%
  #filter(!jobs %in% paste0(job_names, ".txt")) %>%
  distinct(model, .keep_all = TRUE)

table(to_still_run$jobs)
nrow(to_still_run)
list.files(job_dir)

job_dir <- paste0(job_dir, "/cleanup")

unlink(job_dir, recursive = TRUE)

dir.create(job_dir)


group_size <- nrow(to_still_run) %/% 16 + 1

group_size <- min(c(group_size, 48))

to_still_run <- to_still_run %>%
  slice_sample(prop = 1) %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))
  #%>%
  #mutate(group = if_else(grepl("^NPY", model) & !group %in% paste0("job", 1:4, ".txt"), paste0("job", sample(1:4, 1), ".txt"), group)) %>%
  #group_by(group) %>%
  #arrange(if_else(grepl("^NPY", model), 0, 1), .by_group = TRUE)

to_still_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]],
                           file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE,
                           col.names = FALSE,
                           quote = FALSE))

message(job_dir)
out_dir
unique(to_still_run[["group"]] %>% stringr::str_remove(., "^job") %>% stringr::str_remove(., ".txt$") %>% as.numeric(.) %>% max(.))












