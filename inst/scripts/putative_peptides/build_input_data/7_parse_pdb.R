###parse pdb files ~ 1.5 h

pdb_dir <- paste0(s_localDir, "/raw/UP000005640_9606_HUMAN_v4")

parse_pdb <- function(pdb_path) {
  bio3d::read.pdb(pdb_path)[["atom"]] %>%
    mutate(AA = bio3d::aa321(resid)) %>%
    select(all_of(c("AA", "x", "y", "z", "resno", "b", "chain", "elety", "elesy"))) %>%
    nest_by(chain, resno, .key = "atom_level") %>%
    ungroup %>%
    mutate(AA = map_chr(atom_level, ~unique(.[["AA"]])), .before = everything()) %>%
    mutate(pLDDT = map_dbl(atom_level, ~mean(.[["b"]])), .after = "AA") %>%
    mutate(map_df(atom_level, \(x) {bind_cols(x %>%
                                                filter(elety == "CA") %>%
                                                select(all_of(c("x", "y", "z"))) %>%
                                                rename_with(.fn = ~paste0("CA_", .x)),
                                              x %>%
                                                filter(!elety %in% c("CA", "N", "O", "C")) %>%
                                                summarise(across(all_of(c("x", "y", "z")), mean))
    )})) %>%
    mutate(
        x = coalesce(na_if(x, NaN), CA_x),
        y = coalesce(na_if(y, NaN), CA_y),
        z = coalesce(na_if(z, NaN), CA_z)
      )
  
}


start <- Sys.time()

future::plan(strategy = future::multisession(workers = 10))

pdb_dat <- tibble(files = list.files(pdb_dir),
                 accession = str_extract(files, "(?<=-)[^-]+(?=-)")) %>%
          filter(str_detect(files, ".pdb$")) %>%
  mutate(af_dat = furrr::future_map(.x = files, parse_pdb(pdb_path = paste0(pdb_dir, "/", .x))))

end <- Sys.time()

end - start

pdb_dat <- pdb_dat %>%
  mutate(af_dat = map(af_dat, \(x) {
    x %>%
      mutate(
        x = coalesce(na_if(x, NaN), CA_x),
        y = coalesce(na_if(y, NaN), CA_y),
        z = coalesce(na_if(z, NaN), CA_z)
      )
  }))

pdb_dat <- pdb_dat %>%
  rename(af_xyz = af_dat)

saveRDS(pdb_dat, paste0(s_localDir, "/processed/pdb_dat.rds"))










