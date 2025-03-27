

parse_ranges <- function(y, 
                         delim_ranges = "_",
                         delim_start_end = "-") {
  
  if(delim_ranges == delim_start_end) {
    delim_ranges_old <- delim_ranges
    delim_ranges <- "ZSE76FJ"
    y <- sub(delim_ranges_old, delim_ranges, y)
  }
  
  y2 <- stringr::str_split(y, delim_ranges, simplify = TRUE) %>% as_tibble
  logi <- apply(y2, 2, \(z) any(grepl(paste0("\\d+", delim_start_end, "\\d+"), z)))
  to_collapse <- names(logi)[!logi]
  y2 <- y2 %>% 
    rowwise %>%
    mutate(id = paste(!!!rlang::syms(to_collapse), sep = "_")) %>%
    select(-all_of(to_collapse)) %>%
    select(id, everything()) %>%
    ungroup
  colnames(y2) <- c("id", "range")[1:ncol(y2)]
  y2
  
}

parse_proteins <- function(file_name, 
                           delim_proteins = "_and_", 
                           delim_ranges = "_",
                           delim_start_end = "-",
                           protein_names = paste0("p", 1:(1 + max(stringr::str_count(file_name, delim_proteins))))) {
  
  #if(delim_proteins == delim_ranges) {
  delim_proteins_old <- delim_proteins
  delim_proteins <- "CBFHDHE"
  file_name <- sub(delim_proteins_old, delim_proteins, file_name) # only works if wi
}

tmp <- stringr::str_split(file_name, delim_proteins, simplify = TRUE)
colnames(tmp) <- protein_names
tmp %>%
  as_tibble %>%
  mutate(across(everything(), ~parse_ranges(., delim_ranges, delim_start_end), .unpack = TRUE)) %>%
  select(-all_of(protein_names))

}



split_name <- function(x) {
  tmp <- stringr::str_split(x, "_", simplify = TRUE)
  colnames(tmp) <- c("protein", "annotation")
  tmp %>% as_tibble
}

tmp <- tibble(file_name = list.files("~/Desktop/250226_APACE_hGPCRs_hPep_fasta"))


tmp2 <- tmp %>%
  mutate(parse_proteins(file_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>% 
  mutate(across(matches("p\\d+_id"), ~stringr::str_replace(., "^h", ""))) %>%
  mutate(across(matches("p\\d+_range"), \(x) {stringr::str_replace(x, ".fasta$", "") %>% 
                                         stringr::str_replace_all(., "x", "-") })) %>%
  mutate(p1_range = if_else(stringr::str_detect(p1_id, "dSP$"), "dSP", ""), .after = "p1_id") %>%
  mutate(p1_range = if_else(stringr::str_detect(p1_id, "dT$"), "dT", p1_range)) %>% 
  mutate(p1_id =  stringr::str_replace(p1_id, "dT$", "") %>% stringr::str_replace(., "dSP$", "")) %>% 
  {left_join(., gpcr_list %>% rename(p1_id = uniprot_name), by = "p1_id")} %>%
  mutate(p1_range = if_else(p1_range == "dT", paste0((`bw: indices 1.50` - 32), "-", last_AA),
                            if_else(p1_range == "dSP", paste0(signal_peptide, "-", last_AA), ""))) %>% 
  mutate(across(matches("p\\d+_id"), ~setNames(id_map[["Entry"]], id_map[["Entry Name"]])[.])) %>%
  mutate(afpd_models = paste0(p1_id, ",", p1_range, ";", p2_id, ",", p2_range, ";", p3_id, ",", p3_range)) %>%
  mutate(afpd_models = stringr::str_remove(afpd_models, ";$")) %>%
  mutate(afpd_models = stringr::str_remove(afpd_models, ",$")) %>%
  mutate(afpd_models = stringr::str_remove(afpd_models, ";NA$")) %>%
  mutate(afpd_models = stringr::str_remove(afpd_models, ",(?=;)"))

tmp2$afpd_models %>%
  writeLines(., con = "~/Desktop/250226_APACE_hGPCRs_hPep_AFPD.txt")

