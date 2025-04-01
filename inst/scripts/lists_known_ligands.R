



library(dplyr)
library(tidyr)
library(purrr)

pq_path <- paste0(ligandFinder:::get_db_path(), "/residue_db")




gpcrdb_ligands <- jsonlite::read_json("https://gpcrdb.org/services/ligands/endogenousligands")

gpcrdb_ligands <- bind_rows(gpcrdb_ligands) %>%
  filter(ligand_type != "small-molecule") %>%
  distinct(sequence, .keep_all = TRUE) %>%
  mutate(species = sub("^_", "", stringr::str_extract(receptor, "_[^_]*$"))) %>%
  filter(species == "human")




ligand_ids <- map_AA_sequence_to_uniprot(gpcrdb_ligands$sequence, species = "human")

ligand_ids <- ligand_ids %>%
  select(sequence, Subject) %>%
  nest(.by = sequence)

ligand_ids <- ligand_ids %>%
  mutate(num_matches = map_int(data, nrow))

gpcrdb_ligands <- left_join(gpcrdb_ligands, ligand_ids, by = "sequence")

gpcrdb_ligands %>%
  arrange(desc(num_matches)) %>% print(n = 11)

gpcrdb_ligands <- gpcrdb_ligands %>%
  drop_na(sequence)

gpcrdb_ligands[gpcrdb_ligands$sequence == "VYIHPF", "data"][[1]][[1]] <- gpcrdb_ligands[gpcrdb_ligands$sequence == "VYIHPF", "data"][[1]][[1]][1,]
gpcrdb_ligands[gpcrdb_ligands$sequence == "DRVYIHP", "data"][[1]][[1]] <- gpcrdb_ligands[gpcrdb_ligands$sequence == "DRVYIHP", "data"][[1]][[1]][2,]
gpcrdb_ligands[gpcrdb_ligands$sequence == "IPYIL", "data"][[1]][[1]] <- gpcrdb_ligands[gpcrdb_ligands$sequence == "IPYIL", "data"][[1]][[1]][1,]
gpcrdb_ligands[gpcrdb_ligands$sequence == "WMDF", "data"][[1]][[1]] <- gpcrdb_ligands[gpcrdb_ligands$sequence == "WMDF", "data"][[1]][[1]][1,]
gpcrdb_ligands[gpcrdb_ligands$sequence == "LVVYPWT", "data"][[1]][[1]] <- gpcrdb_ligands[gpcrdb_ligands$sequence == "LVVYPWT", "data"][[1]][[1]][1,]
gpcrdb_ligands[gpcrdb_ligands$sequence == "YGGFL", "data"][[1]][[1]] <- gpcrdb_ligands[gpcrdb_ligands$sequence == "YGGFL", "data"][[1]][[1]][4,]

gpcrdb_ligands <- gpcrdb_ligands %>%
  mutate(num_matches = map_int(data, \(x) if(is.null(x)) {return(NA)} else {nrow(x)})) %>%
  mutate(data = map_chr(data, \(x) if(is.null(x)) {return(NA)} else if(nrow(x) > 1) {return(NA)} else {return(x[[1]])}))

split_id <- function(x) {
  setNames(stringr::str_split(x, "\\|")[[1]][2:3],
           c("accession", "uni_name"))
}

gpcrdb_ligands <- gpcrdb_ligands %>%
  mutate(map_df(data, split_id)) %>%
  select(-data)





known_peps <- as_tibble(data.table::fread("~/peptide_alg/guidetopharmacology_peptides.csv")) %>%
  filter(grepl("human", Species, ignore.case = TRUE)) %>%
  dplyr::rename(sequence = `Single letter amino acid sequence`,
                accession = `UniProt id`) %>%
  mutate(ligand_length = nchar(sequence)) %>%
  select(accession, Name, sequence, ligand_length) %>%
  filter(sequence != "" & accession != "") %>%
  mutate(accession = sub("\\|.*", "", accession)) %>%
  mutate(database = if_else(sequence %in% gpcrdb_ligands$sequence, "both", "gtp"))

gpcrdb_ligands <- gpcrdb_ligands %>%
  mutate(database = if_else(sequence %in% known_peps$sequence, "both", "gpcrdb"))

final_ligand <- full_join(known_peps, gpcrdb_ligands, by = "sequence", suffix = c("_gtp", "_gpcrdb"))

final_ligand <- final_ligand %>%
  mutate(database = coalesce(database_gtp, database_gpcrdb), .before = everything()) %>%
  mutate(accession = if_else(!is.na(accession_gpcrdb), accession_gpcrdb, accession_gtp), .before = everything())

final_ligand <- final_ligand %>%
  mutate(final_name = if_else(!is.na(Name), Name, ligand_name))

data.table::fwrite(final_ligand,
                   "~/peptide_alg/gpcrdb_ligand_list.csv")


saveRDS(final_ligand, "~/peptide_alg/final_ligand.rds")














ligands <- readRDS('~/peptide_alg/final_ligand.rds')

ligands <- ligands %>%
              distinct(sequence, .keep_all = TRUE)


ecb_ligands <- data.table::fread("~/peptide_alg/gpcrdb_ligands_gtp and both Kevin KB ECB selected prelim 4 APACE.txt", header = FALSE) %>%
                    as_tibble

ecb_ligands[1,"V4"][[1]] <- "y"

ecb_ligands <- ecb_ligands %>%
  rename(ecb_cull = V4, sequence = V6)


ligands <- left_join(ligands, ecb_ligands %>% select(sequence, ecb_cull), by = "sequence")

replacements <- c("&alpha;" = "alpha",
                  "&beta;" = "beta",
                  "&gamma;" = "gamma",
                  "&delta;" = "delta",
                  "&epsilon;" = "epsilon",
                  "&phi;" = "phi",
                  "&kappa;" = "kappa",
                  "&lambda;" = "lambda",
                  "&omega;" = "omega",
                  "<sub>" = "",
                  "</sub>" = "",
                  "<sup>" = "",
                  "</sup>" = "",
                  " " = "-")

ligands <- ligands %>%
  mutate(in_uniprot = if_else(accession == "", "n", "y"), .after = "database") %>%
  relocate(accession, .before = everything()) %>%
  relocate(final_name, .before = everything()) %>%
  relocate(ecb_cull, .after = "final_name") %>%
  relocate(sequence, .after = "in_uniprot") %>%
  relocate(receptor, .after = "sequence") %>%
  mutate(ecb_cull = if_else(ecb_cull == "gpcrdb" & in_uniprot == "y", "y", ecb_cull)) %>%
  mutate(ecb_cull = factor(ecb_cull,
                           levels = c("y", "n", "gpcrdb"))) %>%
  mutate(final_name = stringr::str_replace_all(final_name, replacements)) %>%
  mutate(database = factor(database,
                           levels = c("both", "gpcrdb", "gtp"))) %>%
  mutate(receptor = toupper(stringr::str_remove(receptor, "_human$"))) %>%
  arrange(ecb_cull, database, desc(in_uniprot), final_name)


id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder", mustWork = TRUE))

ligands <- ligands %>%
  mutate(uniprot_name = setNames(id_map[["Entry Name"]], id_map[["Entry"]])[accession], .before = everything())

res_db <- arrow::open_dataset(source = pq_path)

residue_anno <- res_db %>%
    filter(uni_gene %in% ligands[["uniprot_name"]]) %>%
    group_by(gene_grp) %>%
    collect() %>%
    ungroup()

residue_anno

ligands <- left_join(ligands, residue_anno %>% rename(uniprot_name = uni_gene) %>% select(uniprot_name,
                                                                                          sequence_uni,
                                                                                          secreted_final,
                                                                                          gene_containing_known_ligand), by = "uniprot_name")

ligands <- ligands %>%
  mutate(sequence = toupper(sequence)) %>%
  mutate(map_input = map(sequence, ~tibble::tibble(AA = c(stringr::str_split(., "", simplify = TRUE))))) %>%
  filter(!is.na(sequence_uni)) %>%
  mutate(mapped_peptide_sequence = map2(.x = sequence_uni, .y = map_input, .f = ~map_table(seq2 = .x, to_map = .y)))



ligands <- ligands %>%
  mutate(map_df(mapped_peptide_sequence, ~get_nonNA_ranges(.[["ms"]])), .after = "uniprot_name")


saveRDS(ligands, "~/R_projects/ligandFinder/inst/extdata/ligand_list.rds")




















