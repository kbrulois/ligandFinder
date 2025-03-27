
library(tidyverse)

secretome <- readRDS("~/peptide_alg/secretome_4.rds")

known_peps <- as_tibble(data.table::fread("~/peptide_alg/guidetopharmacology_peptides.csv"))

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

secretome <- secretome %>%
                distinct(accession, .keep_all = TRUE)

problem_peptides <- list()

for(x in unique(final_ligand[["accession"]])) {
  
  inds <- !is.na(final_ligand[["accession"]]) & final_ligand[["accession"]] == x
  
  message("integrating ", x, " ", final_ligand[inds, "final_name"][[1]])
  
  if(!x %in% secretome[["accession"]] | is.na(x)) {
    
    message(x,  " ", final_ligand[inds, "final_name"][[1]], " not found")
    problem_peptides[[x]] <- final_ligand[inds, "final_name"][[1]]
    next
  }
  
  feats <- secretome[secretome[["accession"]] == x, "features"][[1]][[1]] %>% filter(!source %in% c("gtp", "both", "gpcrdb"))
  sequence_ref <- secretome[secretome[["accession"]] == x, "sequence_uni"][[1]]
  sequence_peps <- final_ligand[inds, "sequence"][[1]]
  
  new_feats <- bind_rows(feats,
  lapply(seq_along(sequence_peps), \(y) {
    alignment <- Biostrings::pairwiseAlignment(sequence_peps[y], sequence_ref, type = "local")
    
    
    tibble(type = final_ligand[inds, "final_name"][[1]][y],
           evidence = final_ligand[inds, "database"][[1]][y],
           start = Biostrings::start(Biostrings::subject(alignment)),
           end = Biostrings::end(Biostrings::subject(alignment)),
           source = final_ligand[inds, "database"][[1]][y],
           potency = final_ligand[inds, "potency_ranking"][[1]][y])
  })
  )
  
  secretome[secretome[["accession"]] == x, "features"][[1]][[1]] <- new_feats
  secretome[secretome[["accession"]] == x, "gene_containing_known_ligand"][[1]] <- TRUE
}


saveRDS(secretome, "~/peptide_alg/secretome_4.5.rds")



