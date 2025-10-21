






.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

set_db_path("/scratch/groups/ebutcher/deorphan/ligandFinder")
pq_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db"
voronota_path <- "/home/groups/ebutcher/programs/voronota/bin/voronota-contacts"

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
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))



run_dirs <- c("bm_sep28", "GPCRvCXCL14_oct7", "CXCL14_jh_w", "CXCL14_jh_wo", "CXCL14_mm_w", "CXCL14_mm_wo", "brinp_final", "top200NCnew")

run_dirs <- c("bm_sep28", "GPCRvCXCL14_oct7", "brinp_final", "cxc17_gp15l")


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



runs <- runs[["renamed_dir"]]

runs <- runs %>%
  mutate(parse_proteins(afpd_dir_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  rename(p1_name = p1_id, p2_name = p2_id) %>%
  mutate(file_names = map(afpd_dir, ~fs::dir_ls(.) %>% fs::path_file(.)))


runs <- runs %>%
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

metrics_good <- runs %>%
  filter(has_metrics_csv) %>%
  filter(!map_lgl(metrics, is.null)) %>%
  filter(map_int(metrics, ~nrow(.)) == num_models) %>%
  filter(map_lgl(metrics, ~"depth" %in% colnames(.))) %>%
  pull(afpd_dir_name)


runs <- runs %>%
  mutate(metrics_good = afpd_dir_name %in% !!metrics_good)

table(runs[["metrics_good"]])









runs <- runs %>%
  mutate(contacts = furrr::future_map(afpd_dir, \(x) {
    file <- fs::path(x, "metrics_v2c.rds")
    if(file.exists(file)) {
      return(tryCatch({readRDS(file)}, error = function(e) NULL))
    } else {
      return(NULL)
    }
  }))



contacts_good <- runs %>%
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

runs <- runs %>%
  mutate(contacts_good = afpd_dir_name %in% !!contacts_good)

table(runs[["contacts_good"]])


####add critical columns

ligand_list <- readRDS("/oak/stanford/groups/ebutcher/kevin/ligand_list.rds")

runs <- runs %>%
  {left_join(., ligand_list %>% select(afpd_dir_name, known_pair), by = "afpd_dir_name")} %>%
  mutate(known_pair = if_else(is.na(known_pair), "unknown", known_pair), .after = "afpd_dir")



gpcr_cols <- c("p1_name",
               "ecb: Class or type",
               "ecb: Prioritization Notes",
               "gtp: Family name",
               "gpcrdb: receptor_class",
               "gpcrdb: receptor_family",
               "gpcrdb: subfamily")

runs <- runs %>%
  {left_join(., gpcr_list %>% rename(p1_name = uniprot_name) %>% select(all_of(gpcr_cols)),
             by = "p1_name")} %>%
  relocate(all_of(gpcr_cols[-1]), .before = "p1_name")


runs <- runs %>%
  mutate(gpcr_family = if_else(`gpcrdb: receptor_family` == "", `gtp: Family name`, `gpcrdb: receptor_family`), .after = "gpcrdb: receptor_family") %>%
  mutate(complex_type = paste0(p1_name, "_", p2_name), .after = "afpd_dir_name")



####expand out contacts

runs <- runs %>%
  filter(metrics_good) %>%
  mutate(metrics = map2(metrics, contacts, \(x, y) {
    if(nrow(x) != length(y)) {y <- lapply(1:nrow(x), \(x) {return(NULL)})}
    x[["contacts"]] <- y
    x
  }))

runs <- runs %>%
        select(-contacts)

run_cols <- colnames(runs)

runs_m <- runs %>%
  filter(metrics_good) %>%
  filter(!map_lgl(metrics, is.null)) %>%
  mutate(metrics = map(metrics, ~ select(.x, !any_of(run_cols)))) %>%
  unnest(metrics)

runs_m <- runs_m %>%
  mutate(lig1_end_clean = if_else(lig1_end %in% c("1C", "1N"), lig1_end, "loop"))



runs_c <- runs_m %>%
    mutate(run_name = fs::path_file(run_dir)) %>%
    filter(location == "relevant") %>%
    rename(pLDDT_lig1_og = pLDDT_lig1,
           pLDDT_rec_og = pLDDT_rec) %>%
    unnest(contacts)




####save analysis


out_dir <- "/oak/stanford/groups/ebutcher/kevin"
local_dir <- "~/AF2_analysis"
file_name <- "all_metrics_oct17.csv"

data.table::fwrite(runs_m %>%
                     select(!where(is.list)),
                   paste0(out_dir, "/", file_name))

message("scp kbrulois@dtn.sherlock.stanford.edu:", out_dir, "/", file_name, " ", local_dir, "/", file_name)









