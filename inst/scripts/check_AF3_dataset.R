


run_dir

tmp <- tibble(files = fs::dir_ls(run_dir) %>% basename(),
              file_parts = map(files, ~stringr::str_split(., "_", simplify = TRUE))) %>%
  mutate(file_part_len = map_int(file_parts, length)) %>%
  mutate(parse_proteins(files, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
  mutate(data_files = furrr::future_map(files, ~fs::dir_ls(.) %>% basename()))


tmp <- tmp %>%
  mutate(num_files = map_int(data_files, length))

tmp

tmp <- tmp %>%
  mutate(num_models = map_int(data_files, \(x) {stringr::str_detect(x, "^model_seed-\\d+-sample-\\d+.cif$") %>% sum}))
