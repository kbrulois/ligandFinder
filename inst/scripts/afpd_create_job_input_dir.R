

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
library(ligandFinder)


gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

gpcr_sub <- gpcr_list %>%
              filter(grepl("^#", `ecb: Order of runs (priority)`)) %>%
              filter(`ecb: Exclude due to N term >160AA` != "Exclude due to N term >160AA") %>%
              filter(map_lgl(`bw: full_table`, ~nrow(.) > 0))






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



ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))

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
                      receptor = gpcr_sub[["model_name"]],
                      stringsAsFactors = FALSE) %>%
          as_tibble %>%
          mutate(model = paste(receptor, ligand, sep = ";"))


to_run <- to_run %>%
              mutate(known = if_else(model %in% ligand_list$known_model, "known", "unknown")) %>%
              arrange(known)

group_size <- 48

to_run <- to_run %>%
    mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))


to_run %>%
      group_by(group) %>%
      group_walk(~ write.table(.x[["model_trunc"]], file = .y[["group"]],
                               row.names = FALSE, col.names = FALSE, quote = FALSE))

saveRDS(to_run, "to_run.rds")









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

gpcr_list %>% mutate(lut = setNames(model_name, uniprot_name)) %>% pull(lut) -> gpcr_lut

comp_jobs <- comp_jobs %>%
                mutate(model)

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
                  filter(!model %in% comp_jobs[["model"]]) %>%
                  mutate(parse_proteins(model, delim_proteins = ";", delim_ranges = ","))

left_join(to_still_run, id_map %>% rename(receptor_id = `Entry Name`), by = "receptor_id") %>%
  select(receptor_id, Length) %>%
  print(n = 1000)



group_size <- 48
to_still_run <- to_still_run %>%
  mutate(group = rep(paste0("job", 1:ceiling(n() / group_size), ".txt"), each = group_size, length.out = n()))

to_still_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x$model, file = .y$group,
                           row.names = FALSE, col.names = FALSE, quote = FALSE))



job_script <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh"

submit_model_jobs(script = job_script,
                  input_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/CXCL14vGPCRs/cleanup",
                  output_path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/CXCL14vGPCRs",
                  jobs = unique(to_still_run$group),
                  start_from = 1,
                  max_jobs = 16)












