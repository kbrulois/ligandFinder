

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
library(ligandFinder)




gpcr_sub <- gpcr_list %>%
              filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
              filter(`ecb: Prioritization Notes` != "Small organic molecule") %>% #########caution
              filter(map_lgl(`bw: full_table`, ~nrow(.) > 0)) %>%
              mutate(model = ifelse(`bw: length N-term` > 160, model_name_dNT, model_name))

job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/brinp_w_gai"

dir.create(job_dir)

urls <- c("https://stacks.stanford.edu/file/druid:sc075gg6264/c_terminal.csv",
          "https://stacks.stanford.edu/file/druid:sc075gg6264/n_terminal.csv")

map(urls, ~download_roi_data(url = .))

termini <- c("c", "n")

dat <- bind_rows(map(termini,
                     \(x) data.table::fread(paste0(get_db_path(), "/", x, "_terminal.csv")) %>%
                       as_tibble %>%
                       mutate(terminus = x)))

ligand_list <- dat %>%
  group_by(terminus) %>%
  filter(is.na(!!rlang::sym("percent_ol_phs_hsr:gtp"))) %>%
  mutate(stringr::str_remove(!!rlang::sym("overlap_region_phs:phs_hsr"), "^phs_") %>%
           stringr::str_split_fixed(., "-", n = 2) %>%
           as.data.frame %>%
           rename(start = V1, end = V2) %>%
           mutate(across(everything(), as.integer)), .before = everything()) %>%
  mutate(length2 = end - start + 1, .after = "end") %>%
  filter(length2 > 5 & length2 < 50) %>%
  mutate(uniprot_name = setNames(id_map[["Entry Name"]], id_map[["Entry"]])[accession], .before = everything()) %>%
  mutate(model_id = paste0(accession, ",", start, "-", end), .before = everything()) %>%
  mutate(model_name = paste0(uniprot_name, ",", start, "-", end), .before = everything()) %>%
  slice_max(n = 200, order_by = score_nn10c__entire) %>%
  relocate(terminus, score_nn10c__entire, .after = "length2")


to_run <- expand.grid(ligand = ligand_list[["model_name"]],
                      receptor = gpcr_sub[["model_name"]],
                      stringsAsFactors = FALSE) %>%
          as_tibble %>%
          mutate(model_name = paste(receptor, ligand, sep = ";"))

to_run2 <- expand.grid(ligand = ligand_list[["model_id"]],
                      receptor = gpcr_sub[["model_id"]],
                      stringsAsFactors = FALSE) %>%
  as_tibble %>%
  mutate(model_id = paste(receptor, ligand, sep = ";"))

to_run <- bind_cols(to_run, to_run2 %>% select(-ligand, -receptor))


group_size <- 48

to_run <- to_run %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))


to_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["model_name"]], file = .y[["group"]],
                           row.names = FALSE, col.names = FALSE, quote = FALSE))

data.table::fwrite(to_run, "~/Desktop/n_c_term_candidates_May_2025.csv")



ligand_list <- tibble(uniprot_name = c("GP15L", "GP15L", "CXL17", "CCL25", "RARR2", "RARR2", "RARR2", "RARR2"),
                      start = c(25, 71, 64, 24, 21, 21, 21,21),
                      end = c(81, 81, 119, 150, 163, 157, 158, 156))

ligand_list <- tibble(end = c(87,95:104),
                      uniprot_name = "CXL14",
                      start = 35)

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


ligand_list <- ligand_list %>%
  mutate(model = paste0(uniprot_name, ",", start, "-", end))


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




to_run <- expand.grid(ligand = ligand_list[["model"]],
                      receptor = gpcr_sub[["model"]],
                      stringsAsFactors = FALSE) %>%
          as_tibble %>%
          drop_na() %>%
          mutate(model = paste(receptor, ligand, sep = ";"))


to_run <- to_run %>%
              mutate(known = if_else(model %in% ligand_list$known_model, "known", "unknown")) %>%
              arrange(known)

group_size <- 48

to_run <- to_run %>%
  slice_sample(prop = 1) %>%
    mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n())) %>%
    group_by(group) %>%
  arrange(if_else(grepl("^NPY", model), 0, 1), .by_group = TRUE)

to_run %>%
      group_by(group) %>%
      group_walk(~ write.table(.x[["model"]], file = paste0(job_dir, "/", .y[["group"]]),
                               row.names = FALSE, col.names = FALSE, quote = FALSE))

saveRDS(to_run, "to_run.rds")

job_dir

out_dir <- "/scratch/groups/ebutcher/deorphan/models/brinp_more"

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












