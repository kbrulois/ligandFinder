
library(tidyverse)

secretome <- readRDS(paste0(s_localDir, "/processed/secretome_2.rds"))

known_peps <- as_tibble(data.table::fread(paste0(s_localDir, "/raw/guidetopharmacology_peptides.csv"))

get_accession <- function(accession, Name, id) {
  
  x <- strsplit(accession, split = "\\|")[[1]]
  x <- x[x %in% secretome[["accession"]]]
  if(length(x) == 0) {
    return(tibble(accession_new = "nothing", id = id))
  } else {
    return(tibble(accession_new = x,
                  id = rep(id, length(x))))
  }
}

left_side <- known_peps %>%
  dplyr::rename(accession = `UniProt id`, sequence = `Single letter amino acid sequence`) %>%
  rowwise() %>%
  reframe(get_accession(accession, Name, id = `PubChem SID`))

known_peps <- left_join(left_side, known_peps %>% dplyr::rename(id = `PubChem SID`), by = 'id')

known_peps <- known_peps %>%
  dplyr::rename(sequence = `Single letter amino acid sequence`) %>%
  dplyr::filter(grepl("Human", Species) & accession_new %in% secretome[["accession"]] & sequence != "") %>%
  mutate(ligand_length = nchar(sequence)) %>%
  filter(ligand_length < 60) %>%
  dplyr::rename(accession = accession_new) %>%
  dplyr::select(accession, Name, sequence, ligand_length)


secretome <- secretome %>%
  mutate(gene_containing_known_ligand = FALSE)


for(x in unique(known_peps[["accession"]])) {
  
  feats <- secretome[secretome[["accession"]] == x, "features"][[1]][[1]]
  sequence_ref <- secretome[secretome[["accession"]] == x, "sequence_uni"][[1]]
  sequence_peps <- known_peps[known_peps[["accession"]] == x, "sequence"][[1]]
  
  new_feats <- bind_rows(feats,
                         lapply(seq_along(sequence_peps), \(y) {
                           alignment <- Biostrings::pairwiseAlignment(sequence_peps[y], sequence_ref, type = "local")
                           
                           
                           tibble(type = known_peps[known_peps[["accession"]] == x, "Name"][[1]][y],
                                  evidence = "gtp",
                                  start = Biostrings::start(Biostrings::subject(alignment)),
                                  end = Biostrings::end(Biostrings::subject(alignment)),
                                  source = "gtp")
                         })
  )
  
  secretome[secretome[["accession"]] == x, "features"][[1]][[1]] <- new_feats
  secretome[secretome[["accession"]] == x, "gene_containing_known_ligand"][[1]] <- TRUE
  
}


saveRDS(secretome, paste0(s_localDir, "/processed/secretome_3.rds"))



