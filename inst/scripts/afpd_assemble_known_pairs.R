


.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)

set_db_path("/scratch/groups/ebutcher/deorphan/ligandFinder")
pq_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db"
voronota_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/voronota/bin/voronota-contacts"

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
run_dirs <- c("bm_sep28", "top200NC_ffinal")
run_dirs <- c("bm_sep28", "brinp_final", "top200NC_Nov12", "top200NC_Oct23_cleanup")

run_dirs <- c("bm_sep28", "brinp_Oct28", "CXCL14vGPCRs")

run_dirs <- c("CXCL14peptides")

run_dirs <- "new_peps"

run_dirs <- "bm_sep28"

run_dirs <- list.files(scratch_models)


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






####add critical columns


ligand_list <- data.table::fread(system.file("extdata/GPCRdb_known_pairings_human_plus2more_unique.csv",
                                             package = "ligandFinder")) %>% as_tibble

ligand_list <- ligand_list %>%
  mutate(afpd_dir_name = paste0(rec, "_", lig),
         known_pair = "known")

runs <- runs %>%
  {left_join(., ligand_list %>% select(afpd_dir_name, known_pair), by = "afpd_dir_name")} %>%
  mutate(known_pair = if_else(is.na(known_pair), "unknown", known_pair), .after = "afpd_dir")

kp <- runs %>%
        filter(known_pair == "known") %>%
        distinct(fs::path_file(afpd_dir), .keep_all = TRUE)

kp_dir <- fs::path(scratch_models, "known_pairs")

fs::dir_create(kp_dir, mode = "ug=rwx")

purrr::walk(kp[["afpd_dir"]], ~fs::dir_copy(., kp_dir))






