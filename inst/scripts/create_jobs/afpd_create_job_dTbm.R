






.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)


test <- openxlsx::read.xlsx(xlsxFile = "~/Desktop/ligand_receptor_table.xlsx", sheet = 2) %>% as_tibble()




gpcr_sub <- gpcr_list %>%
  filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))

gpcr_sub <- left_join(test, gpcr_sub %>% select(uniprot_name, dNT_index, first_AA, signal_peptide, last_AA), by = "uniprot_name")

gpcr_sub <- gpcr_sub %>%
  mutate(start = ifelse(start == "NT", 1, start)) %>%
  mutate(end = ifelse(end == "CT", last_AA,end))


afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)

ligand_list <- openxlsx::read.xlsx(xlsxFile = "~/Desktop/ligand_receptor_table.xlsx", sheet = 1) %>% as_tibble()


ligand_list <- ligand_list %>%
  mutate(model = paste0(uniprot_name, ",", start, "-", end))

gpcr_sub <- gpcr_sub %>%
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



job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/dTbm"

dir.create(job_dir)


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


job_dir

to_run$group %>% unique

out_dir <- "/scratch/groups/ebutcher/deorphan/models/dTbm"

dir.create(out_dir)






