





run_name <- "igtest"

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

gpcr_list_new <- data.table::fread("~/R_projects/ligandFinder/inst/extdata/gpcr_selectiong_Sun_Feb_8_for_dendrograms_ECB_hGPCR_GN_perc_ortho.csv") %>% as_tibble

gpcr_list_new <- gpcr_list_new %>%
  filter(decision_Feb8_ecb %in% c("include", "add", "add(check_size)", "include temp", "include_temp", "?? Add for Irina> OR TOO LARGE? Run separatelY"))

gpcr_list_new <- list()

gpcr_list_new$Gene <- c("ADCYAP1R1", "ADGRA1", "ADORA2A", "AGTR1", "AGTR2", "APLNR",
                        "AVPR1A", "AVPR1B", "AVPR2", "BDKRB1", "BDKRB2", "BRS3", "C3AR1",
                        "C5AR1", "C5AR2", "CALCR", "CALCRL", "CCKAR", "CCKBR", "CCR1",
                        "CCR10", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8",
                        "CCR9", "CCRL2", "CHRM1", "CHRM2", "CHRM4", "CHRM5", "CMKLR1",
                        "CMKLR2", "CRHR1", "CRHR2", "CX3CR1", "CXCR1", "CXCR2", "CXCR3",
                        "CXCR4", "CXCR5", "CXCR6", "CYSLTR1", "CYSLTR2", "DRD2", "DRD3",
                        "DRD4", "EDNRA", "EDNRB", "F2R", "F2RL1", "F2RL2", "F2RL3", "FPR1",
                        "FPR2", "FPR3", "GALR1", "GALR2", "GALR3", "GCGR", "GHRHR", "GHSR",
                        "GIPR", "GLP1R", "GLP2R", "GNRHR", "GPBAR1", "GPER1", "GPR132",
                        "GPR139", "GPR141", "GPR142", "GPR146", "GPR148", "GPR149", "GPR15",
                        "GPR150", "GPR151", "GPR152", "GPR160", "GPR17", "GPR171", "GPR173",
                        "GPR176", "GPR18", "GPR183", "GPR20", "GPR22", "GPR25", "GPR27",
                        "GPR31", "GPR32", "GPR33", "GPR34", "GPR35", "GPR37", "GPR37L1",
                        "GPR39", "GPR4", "GPR42", "GPR50", "GPR55", "GPR65", "GPR68",
                        "GPR75", "GPR82", "GPR83", "GPR84", "GPR87", "GPR88", "GPRC5A",
                        "GPRC5B", "GPRC5C", "GPRC5D", "GRPR", "HCAR1", "HCAR2", "HCAR3",
                        "HCRTR1", "HCRTR2", "HRH3", "HRH4", "HTR1B", "HTR1D", "HTR2B",
                        "KISS1R", "LPAR4", "LPAR5", "LPAR6", "LTB4R", "LTB4R2", "MAS1",
                        "MAS1L", "MC1R", "MC2R", "MC3R", "MC4R", "MC5R", "MCHR1", "MCHR2",
                        "MLNR", "MRGPRD", "MRGPRE", "MRGPRF", "MRGPRX1", "MRGPRX2", "MRGPRX3",
                        "MRGPRX4", "MTNR1A", "MTNR1B", "NMBR", "NMUR1", "NMUR2", "NPBWR1",
                        "NPBWR2", "NPFFR1", "NPFFR2", "NPSR1", "NPY1R", "NPY2R", "NPY4R",
                        "NPY5R", "NTSR1", "NTSR2", "OPRD1", "OPRK1", "OPRL1", "OPRM1",
                        "OXER1", "OXGR1", "OXTR", "P2RY1", "P2RY10", "P2RY11", "P2RY12",
                        "P2RY13", "P2RY14", "P2RY2", "P2RY4", "P2RY6", "P2RY8", "PRLHR",
                        "PROKR1", "PROKR2", "PTAFR", "PTGDR2", "PTGER4", "PTH1R", "PTH2R",
                        "QRFPR", "RXFP3", "RXFP4", "SCTR", "SSTR1", "SSTR2", "SSTR3",
                        "SSTR4", "SSTR5", "SUCNR1", "TAAR2", "TAAR5", "TAAR9", "TACR1",
                        "TACR2", "TACR3", "TRHR", "TSHR", "UTS2R", "VIPR1", "VIPR2",
                        "XCR1", "OR51E2", "OR51E1", "OR14C36")

gpcr_sub <- gpcr_list %>%
  filter(gene_name_primary %in% gpcr_list_new$Gene) %>%
  #filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  #filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))

manifest <- data.table::fread("/oak/stanford/groups/ebutcher/deorphan-AI-ze/open_models/open_receptors_manifest.csv") %>% as_tibble

gpcr_sub <- gpcr_sub %>% filter(uniprot_name %in% manifest$receptor)

afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)


ligand_list <- list()


ligand_list[["model"]] <- c("CXL14,88-102", "BRNP2,386-397", "PYY,31-64")

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
  arrange(ligand) %>%
  filter(in_afpd_db) %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))



job_dir <- paste0("/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/", run_name)

dir.create(job_dir)


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


job_dir

out_dir <- paste0("/scratch/groups/ebutcher/deorphan/models/", run_name)

dir.create(out_dir)

total_jobs <- unique(to_run$group) %>% stringr::str_extract(., "\\d+") %>% as.numeric %>% max

total_jobs


