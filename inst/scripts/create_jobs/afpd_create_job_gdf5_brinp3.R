




.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)



gpcr_list_new <- data.table::fread("~/R_projects/ligandFinder/inst/extdata/gpcr_selectiong_Sun_Feb_8_for_dendrograms_ECB_hGPCR_GN_perc_ortho.csv") %>% as_tibble

gpcr_list_new <- gpcr_list_new %>%
  filter(decision_Feb8_ecb %in% c("include", "add", "add(check_size)", "include temp", "include_temp", "?? Add for Irina> OR TOO LARGE? Run separatelY"))


gpcr_sub <- gpcr_list %>%
  filter(gene_name_primary %in% gpcr_list_new$Gene) %>%
  #filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  #filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))



afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)

ligand_list <- tibble(uniprot_name = c("GDF5", "GDF5", "GDF5", "BRNP3", "BRNP3", "BRNP3"),
                      start = c(351, 351, 351, 160, 160, 160),
                      end = c(376, 375, 374, 184, 185, 186),
                      run = "gdf5")


ligand_list <- ligand_list %>%
  mutate(model = paste0(uniprot_name, ",", start, "-", end))


to_run <- expand.grid(ligand = ligand_list[["model"]],
                      receptor = gpcr_sub[["model"]],
                      stringsAsFactors = FALSE) %>%
  as_tibble %>%
  drop_na() %>%
  mutate(model = paste(receptor, ligand, sep = ";"))


to_run <- to_run %>%
  mutate(parse_proteins(model, delim_proteins = ";", delim_ranges = ",", delim_start_end = "-")) %>%
  mutate(in_afpd_db = p1_id %in% afpd_db & p2_id %in% afpd_db)

table(to_run[["in_afpd_db"]])

to_run %>%
  filter(!in_afpd_db) %>% print(n =100)

group_size <- 48

#receptors_first <- c("AGTR1", "AGTR2", "BKRB1", "BKRB2", "APJ", "GPR25", "GPR15", "RXFP1", "RXFP2", "RL3R1", "RL3R2")

to_run <- to_run %>%
  filter(in_afpd_db) %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))



job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/gdf5"

dir.create(job_dir)


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/gdf5"

dir.create(out_dir)






