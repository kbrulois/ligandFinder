

dat <- data.table::fread("~/AF2_analysis/most_recent/bm_update_3_subset_lig_features_coexpression_lig_type_clusters_receptor_features_depcod_new_metrics.csv") %>% as_tibble


dat %>%
  filter(`p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new` == "GP152" & lig1_location == "E") %>%
  arrange(desc(iptm)) %>%
  select(`p2_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new`, p2_range, iptm) %>%
  print(n = 200)


receptors <- c("GP152", "GP132", "GPR31", "MAS", "DRD2", "DRD3", "5HT1E", "HRH3")




pdb_files <- dat %>%
  filter(`p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new` %in% receptors & lig1_location == "E") %>%
  group_by(`p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new`) %>%
  mutate(rank = rank(-iptm, ties.method = "first"),
         rank_rev = rank(iptm, ties.method = "first")) %>%
  filter(rank <= 10 | rank_rev <= 10) %>%
  arrange(rank) %>%
  mutate(iptm_range = c(rep("top10", 10), rep("bottom10", 10)))

pdb_files2 <- dat %>%
  filter(`p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new` == "GP152" & `p2_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new` == "ANGT") %>%
  mutate(iptm_range = "GP152_ANGT")

pdb_files <- bind_rows(pdb_files, pdb_files2)

pdb_files %>% select(iptm, pdb_files, afpd_dir_name, rank, iptm_range) %>% arrange(`p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new`) %>% print(n = 200)


to_download <- pdb_files %>% distinct(iptm_range, `p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new`)

sunet_id <- "kbrulois"

remote_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking"

local_dir_og <- "~/AF2_analysis/top200NC"

dir.create(local_dir_og)


for(i in 1:nrow(to_download)) {

local_dir <- paste0(local_dir_og, "/", paste0(to_download[i, ], collapse = "_"))

dir.create(local_dir)

to_fetch <- tibble(pdb_file = pdb_files %>% filter(iptm_range == to_download[i, "iptm_range"][[1]] &
                                                     `p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new` == to_download[i, "p1_name of bm_update_3_subset_lig_features_coexpression_lig_type_clusters_new"][[1]]) %>% pull(pdb_files),
                   code = stringr::str_extract(pdb_file, "_[A-Z]{4}[a-z0-9]{7}_") %>% str_remove_all(., "_"),
                   directory = stringr::str_extract(pdb_file, "^[^_]*_[^_]*")) %>%
            nest(.by = directory)


scp_commands <- paste("scp",
                      paste0(sunet_id,
                             "@dtn.oak.stanford.edu:",
                             remote_dir,
                             "/",
                             unique(to_fetch[["directory"]]),
                             ".tar"),
                      local_dir)

lapply(scp_commands, system)


res <- to_fetch %>%
  mutate(
    tarfile = file.path(local_dir, paste0(directory, ".tar")),
    matches = map2(tarfile, data, ~{
      tarfile <- .x
      codes   <- .y[["code"]]
      if(length(codes) == 0) return(character(0))

      files_in_tar <- system2("tar", c("-tf", tarfile), stdout = TRUE)
      if(length(files_in_tar) == 0) return(character(0))

      pat <- paste0("ark_.*(", paste0(codes, collapse = "|"), ")")
      files_in_tar[grepl(pat, files_in_tar, perl = TRUE)]
    })
  )

res %>% transmute(directory, n_matches = map_int(matches, length))

res %>% pwalk(function(directory, tarfile, matches, ...) {
  if(length(matches) == 0) return(invisible(NULL))
  outdir <- file.path(local_dir)
  system2("tar", args = c("-xf", tarfile, "--strip-components=1", "-C", outdir, matches))
})

all_items <- list.files(path = local_dir, full.names = TRUE, recursive = TRUE)
pdb_files_to_keep <- all_items[grepl("\\.pdb$", all_items)]
items_to_remove <- setdiff(all_items, pdb_files_to_keep)
unlink(items_to_remove, recursive = TRUE)

}







dat <- data.table::fread("~/AF2_analysis/all_but_top200NC.csv") %>% as_tibble()

out_file <- "/scratch/groups/ebutcher/deorphan/analysis/brinp_cons_cons.csv"

runs_c %>%
  mutate(pep = paste0(p2_name, "_", p2_range)) %>%
  filter((p1_name %in% c("NPY2R", "NPY5R") & pep %in% c("BRNP1_354x368", "BRNP2_386x397", "BRNP3_369x380", "PYY_29x64", "NPY_29x64")) |
         (p1_name %in% c("AGTR1", "AGTR2") & pep %in% c("ANGT_25x32", "ANGT_25x31")) |
         (p1_name %in% c("AGTR1", "AGTR2", "BKRB1", "BKRB2") & pep %in% c("KNG1_380x388")) |
         (p1_name %in% c("BKRB1") & pep %in% c("CXL14_35x102"))) %>%
  filter(rank == 0) %>%
  select(!where(is.list)) %>%
  data.table::fwrite(. , out_file)

message("scp kbrulois@dtn.sherlock.stanford.edu:", out_file, " ", "~/Desktop")




to_download <- clipr::read_clip()

to_download <-stringr::str_extract(to_download, "h\\w+x\\d+x\\d+")

to_download <- unique(to_download[!is.na(to_download)])


#run_dir <- "/scratch/groups/ebutcher/deorphan/models/bm_sep28"
run_dir <- "/scratch/groups/ebutcher/deorphan/models/cxc17_gp15l"

local_dir <- " ~/Desktop"

clipr::write_clip(paste0("scp -r kbrulois@dtn.sherlock.stanford.edu:", run_dir,
"/",
to_download,
#".tar",
local_dir))



test <- clipr::read_clip_tbl()

to_download <- c("hGPR25_hCXL11x22x94")

brnp123 npy2R with or with g_alpha

c("hNPY2R_hBRNP1x356x368", "hNPY2R_hBRNP2x386x397", "hNPY2R_hBRNP3x369x380",
  "hNPY5R_hBRNP1x356x368")

test %>%
  filter(afpd_dir_name...1 == "hNPY5R_hBRNP1x356x368")





