

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)



gpcr_sub <- gpcr_list %>%
              filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
              filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
              filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
              mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/new_pep"

dir.create(job_dir)


id_map <- system.file("data/id_mapping.rds", package = "ligandFinder")
id_map <- readRDS(id_map)
sym2name <- setNames(id_map$`Entry Name`, id_map$`Gene Names (primary)`)
sym2id <- setNames(id_map$`Entry`, id_map$`Gene Names (primary)`)


afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)


ligand_list <- peps %>%
              filter(selection %in% c("secretoglobin", "tm_focused")) %>%
              mutate(uniprot_name = sym2name[gene]) %>%
              filter(uniprot_name %in% afpd_db) %>%
              mutate(model_name = paste0(uniprot_name, ",", start, "-", end))

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


to_run <-  to_run %>%
  filter(in_afpd_db)

group_size <- 48

to_run <- to_run %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))

to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE, col.names = FALSE, quote = FALSE))

job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/new_pep"

dir.create(out_dir)












####Katrin Svensson peptides

id_map <- system.file("data/id_mapping.rds", package = "ligandFinder")
id_map <- readRDS(id_map)
name2id <- setNames(id_map$Entry, id_map$`Entry Name`)
id2name <- setNames(id_map$`Entry Name`, id_map$Entry)

ks_peps <- system.file("extdata/Katrin Svensson THRILL BRP Nature peptide list 2023-08-14576C-Supplementary Table 2.xlsx", package = "ligandFinder")

ks_peps <- openxlsx::read.xlsx(ks_peps, startRow = 3) %>% as_tibble


ks_peps <- ks_peps %>%
  mutate(tibble_split(Protein, "\\|", names = c("sp", "accession", "gene"))) %>%
  mutate(uniprot_name = id2name[accession])

uniprot_goi <- c("FSTL4", "EDIL3", "SCG1", "FGF5", "FGF3", "CHGB")
peptide_num <- c(3, 5, 9, 5, 4, 9)
pep_names <- paste0(uniprot_goi, "_pep", peptide_num)

ks_peps %>%
  filter(uniprot_name %in% uniprot_goi) %>%
  arrange(desc(PeptideRanker_score)) %>%
  group_by(uniprot_name) %>%
  mutate(peptide_num = row_number()) %>%
  mutate(peptide_name = paste0(uniprot_name, "_pep", peptide_num)) %>%
  filter(peptide_name %in% pep_names)

ks_peps2 <- uniprot_t %>%
            filter(gene %in% c(uniprot_goi, "CHGB")) %>%
            select(features, gene) %>%
            mutate(bind_rows(map2(features, gene, \(x, y) {
              x %>%
                filter(source == "sven") %>%
                arrange(desc(PeptideRanker_score)) %>%
                mutate(rank = row_number()) %>%
                mutate(peptide_name = paste0(y, "_pep", rank)) %>%
                filter(peptide_name %in% pep_names) %>%
                select(peptide_name, start, end)
            })))

ks_peps2[4,"gene"][[1]] <- "SCG1"

ligand_list <- ks_peps2 %>%
                dplyr::rename(uniprot_name = gene) %>%
                select(uniprot_name, start, end)










ligand_list <- tibble(uniprot_name = c("GP15L", "GP15L", "CXL17", "CCL25", "RARR2", "RARR2", "RARR2", "RARR2"),
                      start = c(25, 71, 64, 24, 21, 21, 21,21),
                      end = c(81, 81, 119, 150, 163, 157, 158, 156))

ligand_list <- tibble(end = c(95:111),
                      uniprot_name = "CXL14",
                      start = 88)

ligand_list <- tibble(uniprot_name = c("CART", "CART", "CART", "CART", "CQ067", "CQ067", "SDF1", "SDF1", "CO061", "CO061", "CO061", "RARR2", "CCL27", "DMKN", "DMKN", "DMKN"),
                      start = c(28, 64, 76, 93, 19, 38, 71, 71, 29, 31, 45, 29, 77, 460, 355, 395),
                      end = c(61, 73, 90, 116, 35, 90, 82, 88, 42, 42, 157, 45, 88, 476, 376, 408),
                      run = "mxrun")

ligand_list <- tibble(uniprot_name = c("LUZP2", "LUZP2", "LUZP2", "LUZP2"),
                      start = c(157, 123, 22, 327),
                      end = c(167, 154, 62, 346),
                      run = "luzp2")

ligand_list <- tibble(uniprot_name = c("BRNP2", "BRNP2", "BRNP1", "BRNP3", "BRNP2", "BRNP1", "BRNP3"),
                      start = c(372, 347, 356, 369, 386, 356, 369),
                      end = c(397, 397, 368, 380, 397, 368, 380),
                      run = "brinp")

ligand_list <- tibble(uniprot_name = c("BRNP2", "BRNP2", "BRNP1", "BRNP3", "BRNP2", "BRNP1", "BRNP3"),
                      start = c(372, 347, 356, 369, 386, 356, 369),
                      end = c(397, 397, 368, 380, 397, 368, 380),
                      run = "brinp")

ligand_list <- tibble(model = c("BRNP1,317-368", "BRNP1,318-368", "BRNP1,342-368", "BRNP1,354-368",
                                "BRNP3,330-380", "BRNP3,331-380", "BRNP3,355-380", "BRNP3,364-380",
                                "BRNP2,348-397"))

ligand_list <- test %>%
  mutate(V1 = stringr::str_extract(V1, "h\\w+x\\d+x\\d+")) %>%
  mutate(parse_proteins(file_name = V1, delim_ranges = "x", delim_start_end = "x")) %>%
  mutate(p1_range = stringr::str_replace(p1_range, "x", "-")) %>%
  mutate(model = paste0(p1_id, ",", p1_range))

add_bm <- readRDS(system.file("extdata/aug4_ligands.rds", package = "ligandFinder"))
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")

add_bm <- tibble(ligs = add_bm) %>%
          mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
          mutate(uniprot_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
          mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-")) %>%
          mutate(start = stringr::str_extract(p2_range, "^[^-]+") %>% as.numeric) %>%
          mutate(end = stringr::str_remove(p2_range, paste0(start, "-")) %>% as.numeric) %>%
          mutate(aug4_list = "yes") %>%
          mutate(model_name = paste0(uniprot_name, ",", start, "-", end), .before = everything())


ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))


ligand_list <- full_join(ligand_list, add_bm, by = join_by(uniprot_name, start, end))

ligand_list %>%
  filter(ecb_cull == "y" | is.na(ecb_cull)) %>%
  select(uniprot_name, start, end) -> ligand_list_CZ


ligand_list %>%
  filter()

ligand_list <- ligand_list %>%
  mutate(included_in_bm = if_else(ecb_cull == "y" & database %in% c("both", "gpcrdb") & !is.na(start),
                                  "yes", "no"), .after = "final_name")

ligand_list <- ligand_list %>%
                  distinct(uniprot_name, start, end, .keep_all = TRUE)

ligand_list <- ligand_list %>%
                  filter(ecb_cull == "y" & database %in% c("both", "gpcrdb") & !is.na(start))


ligand_list <- left_join(ligand_list, gpcr_list %>% dplyr::rename(receptor = uniprot_name), by = "receptor") %>%
                  mutate(known_model = paste0(model_name, ";", model))


ligand_truncs <- expand.grid(terminus = c("C", "N"),
                             direction = c("minus", "plus"),
                             size = c("1", "2"), stringsAsFactors = FALSE) %>%
                    rowwise %>%
                    mutate(trunc = stringr::str_flatten(c_across(everything()), collapse = "_"))

trunc_funcs <- c(`-`, `+`)
names(trunc_funcs) <- c("minus", "plus")

trunc_type <- c(setNames("start", "N"),
                setNames("end", "C"))

ligand_truncs <- bind_rows(
  map(1:nrow(ligand_truncs), \(x) {
  trunc <- ligand_truncs[x,]
  trunc_func <- function(z, col_name) {trunc_funcs[[trunc[["direction"]]]](z, if(col_name == trunc_type[trunc[["terminus"]]]) {as.numeric(trunc[["size"]])} else {0})}
    ligand_list %>%
      mutate(model_trunc = paste0(uniprot_name, ",", trunc_func(start, "start"), "-", trunc_func(end, "end")), .after = "end") %>%
      mutate(trunc_term = trunc[["terminus"]],
             trunc_dir = trunc[["direction"]],
             trunc_size = trunc[["size"]])
  })
)

to_run <- ligand_truncs %>%
            mutate(model_trunc = paste0(model_name, ";", model_trunc), .after = "model_trunc")

to_run <- to_run %>%
            filter(!grepl(",-1-", model_trunc)) %>%
            filter(!grepl(",0-", model_trunc))

to_run <- left_join(to_run, id_map %>% rename(uniprot_name = `Entry Name`), by = "uniprot_name")

to_run <- to_run %>%
  filter(Length.y >= end) %>%
  filter(start > 0)





neuro <- openxlsx::read.xlsx(system.file("extdata/MultiPep_Scores_Top1068_EXTRA_BMP8B_to_COL16A1_updated_peptides_2025-10-20T20-06-57-002Z_html orig.xlsx",
                                         package = "ligandFinder"))

colnames(neuro)[c(9, 32)] <- c("Peptide2", "GeneName2")

ligand_list <- neuro %>%
  as_tibble %>%
  mutate(p1_range = stringr::str_extract(AF2.Code, "x\\d+x\\d+") %>%
                    stringr::str_remove(., "^x") %>%
                    stringr::str_replace(., "x", "-"),
         .before = everything()) %>%
  mutate(p1_name = stringr::str_remove(AF2.Code, "x\\d+x\\d+") %>%
                   stringr::str_remove(., "^h"), .before = everything()) %>%
  drop_na(p1_name) %>%
  mutate(model = paste0(p1_name, ",", p1_range)) %>%
  mutate(order = 1:nrow(.), .before = everything())


afpd_db <- tibble(files = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens"))

afpd_db <- afpd_db %>% filter(grepl(".pkl.xz$", files)) %>% mutate(uniprot_name = stringr::str_remove(files, ".pkl.xz$")) %>% pull(uniprot_name)




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

to_run <-  to_run %>%
  filter(in_afpd_db)

to_run <- left_join(to_run, ligand_list %>% select(model, order) %>% rename(ligand = model), by = "ligand") %>%
  arrange(order)


group_size <- 48

receptors_first <- c("AGTR1", "AGTR2", "BKRB1", "BKRB2", "APJ", "GPR25", "GPR15", "RXFP1", "RXFP2", "RL3R1", "RL3R2")

to_run <- to_run %>%
  #slice_sample(prop = 1) %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n())) %>%
  #group_by(group) %>%
  #arrange(if_else(p1_id %in% receptors_first, 0, 1))

to_run %>%
      group_by(group) %>%
      group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                               row.names = FALSE, col.names = FALSE, quote = FALSE))

saveRDS(to_run, "to_run.rds")

job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/new_pep"

dir.create(out_dir)






input_path_models <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/benchmarking"

comp_jobs <- parse_dirname(run_dir = input_path_models,
                           delim_proteins = "_",
                           delim_ranges = "x",
                           delim_start_end = "x") %>%
  mutate(parsed_pair = map(parsed_pair, ~pivot_wider(., names_from=c("protein", "annotation"), values_from=value))) %>%
  unnest(parsed_pair)

comp_jobs <- comp_jobs %>%
  mutate(complete = map_lgl(afpd_dir_name, \(x) {
    file.exists(paste(input_path_models, x, "ranking_debug.json", sep = "/"))
  }))

yo()


comp_jobs %>%
  filter(complete)

to_run <- to_run %>%
  filter(!model_trunc %in% comp_jobs$model_name)


comp_jobs <- tibble(afpd_file_name = list.files("/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/CXCL14vGPCRs"))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

comp_jobs <- comp_jobs %>%
  mutate(parse_proteins(afpd_file_name)) %>%
  mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
  mutate(model = paste0(p1_id,
                        if_else(p1_range == "", "", paste0(",", p1_range)),
                        ";",
                        p2_id,
                        if_else(p2_range == "", "", paste0(",", p2_range))))

sum(to_run[["model"]] %in% comp_jobs[["model"]])


to_still_run <- to_run %>%
                  filter(!model %in% (comp_jobs %>% filter(complete) %>% pull(model_name)))


left_join(to_still_run, id_map %>% rename(receptor_id = `Entry Name`), by = "receptor_id") %>%
  select(receptor_id, Length) %>%
  print(n = 1000)



group_size <- 48
to_still_run <- to_still_run %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))

to_still_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model"]],
                           file = .y[["group"]],
                           row.names = FALSE,
                           col.names = FALSE,
                           quote = FALSE))



job_script <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh"

submit_model_jobs(script = job_script,
                  input_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/CXCL14vGPCRs/cleanup",
                  output_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/CXCL14vGPCRs",
                  jobs = unique(to_still_run$group),
                  start_from = 1,
                  max_jobs = 16)












