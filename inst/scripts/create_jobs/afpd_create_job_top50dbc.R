



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)




gpcr_sub <- gpcr_list %>%
  filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
  filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
  filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
  mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))


afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)

ligand_list <- structure(list(uniprot_name = c(TECTA = "TECTA", BRINP1 = "BRNP1",
                                               GRIA2 = "GRIA2", IFNW1 = "IFNW1", HMGB1 = "HMGB1", `HLA-DRB5` = "DRB5",
                                               HSPA1A = "HS71A", KRT10 = "K1C10", SFN = "1433S", GPC6 = "GPC6",
                                               GRIA1 = "GRIA1", SPINK5 = "ISK5", HMGB1 = "HMGB1", INSL6 = "INSL6",
                                               HSP90AB1 = "HS90B", VGF = "VGF", C4B = "A0A140TA29", LAD1 = "LAD1",
                                               VASH1 = "VASH1", ENOX1 = "ENOX1", GRID2 = "GRID2", GRID1 = "GRID1",
                                               HMGB2 = "HMGB2", GLG1 = "GSLG1", NLRP3 = "A0A7I2R3P8", CCDC126 = "CC126",
                                               GRM8 = "GRM8", SAMD1 = "SAMD1", SFTPA1 = "SFTA1", LAD1 = "LAD1",
                                               YBX1 = "YBOX1", LTBP1 = "LTBP1", TMTC4 = "TMTC4", TMTC3 = "TMTC3",
                                               CHST9 = "CHST9", FGFBP3 = "FGFP3", PIBF1 = "PIBF1", NUCB1 = "NUCB1",
                                               ENOX1 = "ENOX1", ADAMTS18 = "ATS18", HMGB1 = "HMGB1", ITPRIPL1 = "IPIL1",
                                               PI16 = "PI16", ARSI = "ARSI", SIDT2 = "SIDT2", IGFBP5 = "IBP5",
                                               ANXA2 = "ANXA2", PHEX = "PHEX", TMTC2 = "TMTC2", PKDCC = "PKDCC"
), start = c(1709L, 291L, 294L, 128L, 108L, 103L, 242L, 308L,
             115L, 188L, 287L, 188L, 156L, 184L, 431L, 291L, 724L, 15L, 120L,
             221L, 316L, 324L, 14L, 436L, 59L, 61L, 359L, 22L, 157L, 3L, 222L,
             720L, 681L, 663L, 295L, 89L, 548L, 170L, 354L, 436L, 14L, 47L,
             60L, 529L, 577L, 168L, 221L, 572L, 710L, 65L), end = c(1725L,
                                                                    314L, 322L, 156L, 126L, 121L, 270L, 332L, 139L, 204L, 315L, 216L,
                                                                    171L, 198L, 447L, 307L, 742L, 29L, 144L, 241L, 344L, 352L, 42L,
                                                                    456L, 87L, 82L, 375L, 44L, 178L, 31L, 245L, 741L, 705L, 682L,
                                                                    310L, 117L, 570L, 198L, 369L, 453L, 42L, 61L, 83L, 557L, 596L,
                                                                    192L, 244L, 590L, 734L, 83L)), row.names = c(NA, -50L), class = c("tbl_df",
                                                                                                                                      "tbl", "data.frame")) %>%
  mutate(run = "top50dbc")


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



job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/top50dbc"

dir.create(job_dir)


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/top50dbc"

dir.create(out_dir)






