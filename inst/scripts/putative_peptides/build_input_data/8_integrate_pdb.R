
#~10min

secretome <- readRDS(paste0(s_localDir, "/processed/secretome_4.rds"))

af_xyz <- readRDS(paste0(s_localDir, "/processed/pdb_dat.rds"))


af_xyz <- af_xyz %>%
  mutate(sequence_afxyz = map_chr(af_xyz, ~paste(.[["AA"]], collapse = "")))

af_xyz <- af_xyz %>%
  mutate(af_xyz = map(af_xyz, \(x) {x[["atom_level"]] <- as.list(x[["atom_level"]]); return(x)}))

secretome <- left_join(secretome, af_xyz, by = "accession")

secretome %>%
  mutate(same_seq = sequence_uni == sequence_afxyz) %>% pull(same_seq) %>% sum(., na.rm = TRUE)

 

start <- Sys.time()

secretome <- secretome %>%
  mutate(af_xyz = pmap(.l = list(seq1 = sequence_afxyz, 
                                 seq2 = sequence_uni, 
                                 to_map = af_xyz),
                       .f = map_table))

end <- Sys.time()
end - start

saveRDS(secretome, paste0(s_localDir, "/processed/secretome_5.rds"))



