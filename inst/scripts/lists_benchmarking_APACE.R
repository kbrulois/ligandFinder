
id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

test <- readRDS("~/oak/deorphan-AI-ze/scripts/benchmarking/to_run.rds")

test2 <- test %>%
  mutate(parse_proteins(model,
                        delim_proteins = ";",
                        delim_ranges = ",",
                        delim_start_end = "-")) %>%
  mutate(across(ends_with("_id"), ~setNames(id_map[["Entry"]], id_map[["Entry Name"]])[.])) %>%
  mutate(model_id = case_when(p1_range == "" ~ paste0(p1_id, ";", p2_id, ",", p2_range),
                              p1_range != "" ~ paste0(p1_id, ",", p1_range, ";", p2_id, ",", p2_range))) %>%
  dplyr::rename(model_name = model)


data.table::fwrite(test2, "~/Desktop/benchmarking_models.csv")
