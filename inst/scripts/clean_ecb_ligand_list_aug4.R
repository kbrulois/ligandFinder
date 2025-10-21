


id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

known_pairs <- get_known_pairs()

known_pairs_tbl <- do.call(rbind, known_pairs) %>% as_tibble



dat <- tibble(og_ligand = dat)


dat <- dat %>%
  mutate(p2_name = stringr::str_remove(og_ligand, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
  mutate(p2_range = stringr::str_extract(og_ligand, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-"))



dat <- dat %>%
  mutate(p2_id = setNames(id_map$Entry, id_map$`Entry Name`)[p2_name]) %>%
  mutate(p2_gene = setNames(id_map$`Gene Names (primary)`, id_map$`Entry Name`)[p2_name]) %>%
  rowwise %>%
  mutate(known_receptors = paste0(known_pairs_tbl[known_pairs_tbl$V1 == p2_name, "V2"], collapse = "; "))


full_tbl <- clipr::read_clip_tbl()

full_tbl <- as_tibble(full_tbl)

full_tbl <- bind_cols(full_tbl, dat %>% select(-og_ligand))
data.table::fwrite(full_tbl, "~/Desktop/ligand_tbl.csv")



add_bm
ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))
ligand_list
add_bm
saveRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
saveRDS(add_bm, "~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- readRDS(system.file("extdata/aug4_ligands.rds", package = "ligandFinder"))
devtools::load_all(".")
add_bm <- readRDS(system.file("extdata/aug4_ligands.rds", package = "ligandFinder"))
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm
add_bm <- tibble(ligs = add_bm) %>%
mutate(p2_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-"))
add_bm
View(add_bm)
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == "hANF104x151", "hANFx104x151", ligs)) %>%
mutate(p2_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-"))
add_bm
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
mutate(p2_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-"))
add_bm
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
mutate(p2_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-")) %>%
mutate(aug4_list = "yes")
se_known2
se_known <- data.table::fread("~/AF2_analysis/most_recent/2025_07_23_UCSD508_pairings.csv") %>% as_tibble()
se_receptors <- se_known$rec %>% stringr::str_remove(., "^h")
se_known <- se_known %>%
mutate(parse_proteins(paste0(rec, "_", lig), delim_proteins = "_", delim_ranges = "x", delim_start_end = "x", num_proteins = 2))
library(ligandFinder)
se_known <- se_known %>%
mutate(parse_proteins(paste0(rec, "_", lig), delim_proteins = "_", delim_ranges = "x", delim_start_end = "x", num_proteins = 2))
se_known2 <- lapply(1:nrow(se_known), \(x) c(se_known[x, "p1_id"][[1]], se_known[x, "p2_id"][[1]]))
se_known2
se_known
ligand_list <- full_join(ligand_list, add_bm, by = join_by(uniprot_name, start, end))
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
mutate(uniprot_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-")) %>%
mutate(start = stringr::str_extract(p2_range, "^[^-]+")) %>%
mutate(end = stringr::str_extract(p2_range, "^[^-]+-(.*)$")) %>%
mutate(aug4_list = "yes")
ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))
ligand_list <- full_join(ligand_list, add_bm, by = join_by(uniprot_name, start, end))
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
mutate(uniprot_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-")) %>%
mutate(start = stringr::str_extract(p2_range, "^[^-]+") %>% as.numeric) %>%
mutate(end = stringr::str_extract(p2_range, "^[^-]+-(.*)$") %>% as.numeric) %>%
mutate(aug4_list = "yes")
ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))
ligand_list <- full_join(ligand_list, add_bm, by = join_by(uniprot_name, start, end))
View(ligand_list)
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
mutate(uniprot_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-")) %>%
mutate(start = stringr::str_extract(p2_range, "^[^-]+") %>% as.numeric) %>%
mutate(end = stringr::str_extract(p2_range, "^[^-]+-(.*)$") %>% as.numeric) %>%
mutate(aug4_list = "yes")
add_bm <- readRDS("~/R_projects/ligandFinder/inst/extdata/aug4_ligands.rds")
add_bm <- tibble(ligs = add_bm) %>%
mutate(ligs = if_else(ligs == ">hANF104x151", ">hANFx104x151", ligs)) %>%
mutate(uniprot_name = stringr::str_remove(ligs, "^>h") %>% stringr::str_remove(., "x\\d+x\\d+")) %>%
mutate(p2_range = stringr::str_extract(ligs, "x\\d+x\\d+") %>% stringr::str_remove(., "^x") %>% stringr::str_replace(., "x", "-")) %>%
mutate(start = stringr::str_extract(p2_range, "^[^-]+") %>% as.numeric) %>%
mutate(end = stringr::str_remove(p2_range, paste0(start, "-")) %>% as.numeric) %>%
mutate(aug4_list = "yes")
ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))
ligand_list <- full_join(ligand_list, add_bm, by = join_by(uniprot_name, start, end))
View(ligand_list)
ligand_list %>%
filter(ecb_cull == "yes" | is.na(ecb_cull))
ligand_list %>%
filter(ecb_cull == "yes")
table(ligand_list$ecb_cull)
ligand_list %>%
filter(ecb_cull == "y" | is.na(ecb_cull))
ligand_list %>%
filter(ecb_cull == "y" | is.na(ecb_cull)) %>%
select(uniprot_name, start, end) -> ligand_list_CZ
write.csv(ligand_list_CZ, "~/Desktop/ligand_list_cz.csv")




