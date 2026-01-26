
id_map <- readRDS(system.file("/data/id_mapping.rds", package = "ligandFinder"))

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
