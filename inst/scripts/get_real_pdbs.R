
library(tidyverse)
library(ligandFinder)

real_pdbs <- jsonlite::fromJSON("https://gpcrdb.org/services/structure/", )


real_pdbs <- real_pdbs %>%
  as_tibble %>%
  dplyr::rename(method = type) %>%
  unnest_wider(ligands, names_sep = "_") %>%
  unnest_wider(signalling_protein) %>%
  unnest_wider(data) %>%
  unnest_wider(c(entity1, entity2, entity3), names_sep = "_")


table(real_pdbs$species)

real_pdbs <- real_pdbs %>%
  filter(species == "Homo sapiens")

table(real_pdbs$class)

real_pdbs <- real_pdbs %>%
  mutate(peptide_ligand = map_lgl(ligands_type, \(x) {
    any(x %in% c("protein", "peptide"))
  }))

real_pdb2 <- tibble(files = fs::dir_ls("~/AF2_analysis/PDB_wPeptides") %>% basename)

real_pdb2 <- real_pdb2 %>%
  mutate(map_df(files, \(x) {
    x <- stringr::str_remove(x, "_pdb$")
    y <- stringr::str_split(x, "_", simplify = TRUE)
    colnames(y) <- c("receptor", "pdb_code", "chain_SE")
    as.data.frame(y)
  })) %>%
  mutate(pdb_code = toupper(pdb_code))

real_pdbs <- left_join(real_pdbs, real_pdb2, by = "pdb_code")

real_pdbs %>%
  filter(peptide_ligand) %>%
  View()

real_pdb2$pdb_code[!real_pdb2$pdb_code %in% real_pdbs$pdb_code]



pdb_dir <- "~/AF2_analysis/real_pdbs"

dir.create(pdb_dir)

real_pdbs %>%
  filter(peptide_ligand) %>%
  mutate(pdb = map(pdb_code, ~bio3d::get.pdb(ids = ., path = pdb_dir)))



pep_pdbs <- real_pdbs %>%
  mutate(pdb_file = stringr::str_c(pdb_dir, "/", toupper(pdb_code), ".pdb")) %>%
  filter(file.exists(pdb_file)) %>%
  filter(peptide_ligand) %>%
  mutate(pdb = map(pdb_file, ~bio3d::read.pdb(file = .)))


pep_pdbs <- pep_pdbs %>%
  filter(pdb_code != "7XOX")

num_of_grps <- 10

future::plan(strategy = future::multisession(workers = num_of_grps))
options(future.globals.maxSize = +Inf)

pep_pdbs <- pep_pdbs %>%
  mutate(pdb_parsed = furrr::future_map(pdb, parse_pdb))

pep_pdbs <- pep_pdbs %>%
  mutate(chain_names = map(pdb_parsed, \(x) {

    x %>%
      group_by(chain) %>%
      summarise(sequence = stringr::str_c(AA, collapse = "")) %>%
      mutate(uniprot = map(sequence, ~map_AA_sequence_to_uniprot(x = .))) %>%
      unnest(uniprot, names_sep = "_")

  }))


pep_pdbs %>%
  unnest(chain_names) %>%
  select(pdb_code, protein, preferred_chain, chain, chain_SE, uniprot_Subject) %>%
  View()



