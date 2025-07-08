


run_dir <- "/scratch/groups/ebutcher/deorphan/models/benchmarking_AF3"

tmp <- tibble(files = fs::dir_ls(run_dir) %>% basename(),
              file_parts = map(files, ~stringr::str_split(., "_", simplify = TRUE))) %>%
  mutate(file_part_len = map_int(file_parts, length)) %>%
  mutate(parse_proteins(files, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  mutate(data_files = furrr::future_map(files, ~fs::dir_ls(paste0(run_dir, "/", .)) %>% basename()))


tmp <- tmp %>%
  mutate(num_files = map_int(data_files, length))

tmp

tmp <- tmp %>%
  mutate(num_models = map_int(data_files, \(x) {stringr::str_detect(x, "^model_seed-\\d+-sample-\\d+.cif$") %>% sum}))

tmp2 <- tmp %>%
          filter(num_models == 10)

tmp <- tmp %>%
        mutate(ligand_id = paste0(p2_id, "x", p2_range))

pairing_lib <- expand.grid(receptors = unique(tmp$p1_id),
                           ligands = unique(tmp$ligand_id)) %>%
                mutate(pairing_id = paste0(receptors, "_", ligands))

tmp <- tmp %>%
        mutate(pairing_id = paste0(p1_id, "_", ligand_id))

sum(!pairing_lib$pairing_id %in% tmp$pairing_id)
