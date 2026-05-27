




.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)



gpcr_list_new <- data.table::fread("~/R_projects/ligandFinder/inst/extdata/gpcr_selectiong_Sun_Feb_8_for_dendrograms_ECB_hGPCR_GN_perc_ortho.csv") %>% as_tibble

gpcr_list_new <- data.table::fread(system.file("/inst/extdata/gpcr_selectiong_Sun_Feb_8_for_dendrograms_ECB_hGPCR_GN_perc_ortho.csv", package = "ligandFinder")) %>% as_tibble

gpcr_list_new <- gpcr_list_new %>%
                  filter(decision_Feb8_ecb %in% c("include", "add", "add(check_size)", "include temp", "include_temp", "?? Add for Irina> OR TOO LARGE? Run separatelY"))


gpcr_sub <- gpcr_list %>%
              filter(gene_name_primary %in% gpcr_list_new$Gene) %>%
              #filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
              #filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
              filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
              mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))


ligand_list <- data.table::fread("~/R_projects/ligandFinder/inst/extdata/Candidate_peps_KB_from_known_pep_unknownGPCR_run#1.txt") %>% as_tibble

ligand_list <- ligand_list %>%
                  filter(!ALso_run %in% c("no", "not_now", "??", "did we run?"))

ligand_list <- ligand_list %>%
                  mutate(uniprot_name = setNames(id_map$`Entry Name`, id_map$Entry)[Accession]) %>%
                  mutate(model = paste0(uniprot_name, ",", Start, "-", End), .before = everything())

ligand_list %>%
  select(model) %>%
  dput()

gpcr_sub %>%
  select(model) %>%
  dput()

afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)



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



job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/uni_pep"

dir.create(job_dir)


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))


job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/uni_pep"

dir.create(out_dir)














secretome <- readRDS("~/AF2_analysis/secretome_latest.rds")

extract_roi <- \(x, y) {
  rois <- x %>%
    filter(source %in% c("gpcrdb_gtp", "sven", "phs", "top200NC", "disha", "uni_pep")) %>%
    dplyr::rename("roi_name" = "type", "roi_type" = "source")

  dummy_table <- tibble(accession = y, roi_name = NA, roi_type = NA)
  if(nrow(rois) > 0) {
    return(bind_cols(accession = rep(y, nrow(rois)), rois))
  } else {
    return(dummy_table)
  }
}

secretome <- secretome %>%
                mutate(features = map(features, .f = ~define_ol(x = ., a = "uni_pep", b = "gpcrdb_gtp")))


left_side <- secretome %>%
  rowwise() %>%
  reframe(extract_roi(features, accession))

left_side <- left_side %>%
          mutate(pep_id = paste0(accession, "_", start, "-", end)) %>%
          mutate(pep_id_N = paste0(accession, "_", start)) %>%
          mutate(pep_id_C = paste0(accession, "_", end))

kps <- left_side %>%
        filter(roi_type %in% c("top200NC", "gpcrdb_gtp"))

test <- left_side %>%
          filter(roi_type  == "uni_pep") %>%
          mutate(known_similarity = case_when(pep_id %in% kps[["pep_id"]] ~ "identical",
                                              pep_id_N %in% kps[["pep_id_N"]] ~ "known_start",
                                              pep_id_C %in% kps[["pep_id_C"]] ~ "known_end",
                                              TRUE ~ ""))

test <- test %>%
  filter(known_similarity == "")


test <- test %>%
          distinct(accession, start, end, .keep_all = TRUE)



test <- left_join(test, secretome %>% select(accession, gene), by = "accession")


data.table::fwrite(test %>% select(285, 249, 250, 243, 211, 1:7), file = "~/Desktop/uniprot_peptides.csv")
























