






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


list.files(scratch_models)

run_dirs <- c("CXCL14peptides", "cxcl14spep", "jan30", "uni_pep", "gdf5", "new_peps")

run_dirs <- c("top50dbc")



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


####add critical columns


ligand_list <- data.table::fread(system.file("extdata/GPCRdb_known_pairings_human_plus2more_unique.csv",
                                             package = "ligandFinder")) %>% as_tibble

ligand_list <- ligand_list %>%
  mutate(afpd_dir_name = paste0(rec, "_", lig),
         known_pair = "known")

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



all_ligands <- runs[["p2_name"]] %>% unique(.)

res_db <- arrow::open_dataset(source = pq_path)

residue_data <- res_db %>%
  filter(uni_gene %in% all_ligands) %>%
  select(uni_gene, sequence_uni) %>%
  collect()



ligs <- runs %>%
  distinct(p2_name, p2_range) %>%
  separate_wider_delim(p2_range, delim = "x", names = c("start", "end")) %>%
  mutate(across(all_of(c("start", "end")), as.numeric)) %>%
  {left_join(., residue_data, by = join_by(p2_name == uni_gene))}


find_adjacent_dibasic <- function(start,
                                  end,
                                  sequence_uni) {

  start_seq <- stringr::str_sub(sequence_uni, start = max(c(1, start - 5)), end = max(c(1, start - 1)))

  start_gap <- stringr::str_locate(stringi::stri_reverse(start_seq), "KK|KR|RK|RR")[1] - 1

  term <- nchar(sequence_uni)

  end_seq <- stringr::str_sub(sequence_uni, start = min(c(term, end + 1)), end = min(c(term, end + 5)))

  end_gap <- stringr::str_locate(end_seq, "KK|KR|RK|RR")[1] - 1

  tibble(N_term_db = start_gap,
         C_term_db = end_gap)
}

ligs <- ligs %>%
  mutate(dibasic = pmap(.l = list(start, end, sequence_uni), .f = find_adjacent_dibasic)) %>%
  mutate(bind_rows(dibasic)) %>%
  mutate(p2_range = stringr::str_c(start, end, sep = "x")) %>%
  select(-dibasic, -sequence_uni, -start, -end)

runs <- runs %>%
  {left_join(., ligs, by = join_by(p2_name, p2_range))}






####expand out models

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
  mutate(lig1_end_clean = if_else(lig1_end %in% c("1C", "1N"), stringr::str_remove(lig1_end, "^1"), "L"))

runs_m <- runs_m %>%
  mutate(dibasic = case_when(!is.na(N_term_db) & !is.na(C_term_db) ~ "Bphs",
                             is.na(N_term_db) & !is.na(C_term_db) ~ "Cphs",
                             !is.na(N_term_db) & !is.na(C_term_db) ~ "B",
                             is.na(N_term_db) & !is.na(C_term_db) ~ "C",
                             !is.na(N_term_db) & is.na(C_term_db) ~ "N",
                             TRUE ~ "X")) %>%
  mutate(lig1_insert = paste0(lig1_end_clean, dibasic))










####save analysis


out_dir <- "/oak/stanford/groups/ebutcher/kevin"
local_dir <- "~/AF2_analysis"
file_name <- "top50dbc"

data.table::fwrite(runs_m %>%
                     select(!where(is.list)),
                   paste0(out_dir, "/", file_name, ".csv"))

message("scp kbrulois@dtn.sherlock.stanford.edu:", out_dir, "/", file_name, ".csv", " ", local_dir, "/", file_name, ".csv")






######archive metrics

archive_metrics(runs)

yo("archive complete")





####expand out contacts

runs_c <- runs_m %>%
  mutate(run_name = fs::path_file(run_dir)) %>%
  filter(location == "relevant") %>%
  rename(pLDDT_lig1_og = pLDDT_lig1,
         pLDDT_rec_og = pLDDT_rec) %>%
  unnest(contacts)



saveRDS(runs_m, paste0(out_dir, "/", file_name, ".rds"))

message("scp kbrulois@dtn.sherlock.stanford.edu:", out_dir, "/", file_name, " ", local_dir, "/", file_name, ".rds")






