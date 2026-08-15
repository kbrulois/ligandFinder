



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)

genes <- c("BRS3", "FSHR", "G37L1", "GPR37", "GRPR", "LGR4",
           "LGR5", "LGR6", "LSHR", "NMBR", "NPY42", "OPSG2", "OPSG3", "OPSX",
           "RGR", "RXFP1", "RXFP2", "TSHR")

gpcr_sub <- gpcr_list %>%
  #filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  #filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  #filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  #mutate(model = ifelse(`bw: length N-term` > 160 & !is.na(`bw: length N-term`), model_name_dNT, model_name)) %>%
  filter(uniprot_name %in% genes) %>%
  #mutate(has_bw = map_lgl(`bw: full_table`, \(x) nrow(.) > 0)) %>%
  select(model_name, signal_peptide, model_name_dCT)



afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)


ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))


ligand_list <- ligand_list %>%
  distinct(uniprot_name, start, end, .keep_all = TRUE)

ligand_list <- ligand_list %>%
  filter(ecb_cull == "y" & database %in% c("both", "gpcrdb") & !is.na(start))

ligand_list <- ligand_list %>%
                  mutate(model = paste0(uniprot_name, ",", start, "-", end))

#bm_dat <- data.table::fread("~/AF2_analysis/bm_sep28.csv")
#bm_ligs <- unique(paste0(bm_dat$p2_id, "_", bm_dat$p2_range))
#ll_ligs <- unique(paste0(ligand_list$uniprot_name, "_", ligand_list$start, "x", ligand_list$end))
#setdiff(bm_ligs, ll_ligs)
#intersect(bm_ligs, ll_ligs)

to_run <- expand.grid(ligand = ligand_list[["model"]],
                      receptor = gpcr_sub[["model_name"]],
                      stringsAsFactors = FALSE) %>%
  as_tibble %>%
  drop_na() %>%
  mutate(model = paste(receptor, ligand, sep = ";"))


to_run <- to_run %>%
  mutate(parse_proteins(model, delim_proteins = ";", delim_ranges = ",", delim_start_end = "-")) %>%
  mutate(in_afpd_db = p1_id %in% afpd_db & p2_id %in% afpd_db)

to_run <- to_run %>%
          filter(!p1_id %in% c("RXFP1", "RXFP2"))

table(to_run[["in_afpd_db"]])

group_size <- 48

#receptors_first <- c("AGTR1", "AGTR2", "BKRB1", "BKRB2", "APJ", "GPR25", "GPR15", "RXFP1", "RXFP2", "RL3R1", "RL3R2")

to_run <- to_run %>%
  #slice_sample(prop = 1) %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))
  #group_by(group) %>%
  #arrange(if_else(p1_id %in% receptors_first, 0, 1))


job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/bm_more_rec"

dir.create(job_dir)


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/bm_more_rec"

dir.create(out_dir)






