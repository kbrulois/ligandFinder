

library(tidyverse)

s_localDir <- "~/peptide_alg/build_residue_db"

uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_3.rds"))

id_map <- system.file("data/id_mapping.rds", package = "ligandFinder")
id_map <- readRDS(id_map)
name2id <- setNames(id_map$Entry, id_map$`Entry Name`)
id2name <- setNames(id_map$`Entry Name`, id_map$Entry)



####integrate Katrin Svensson peptides

ks_peps <- system.file("extdata/Katrin Svensson THRILL BRP Nature peptide list 2023-08-14576C-Supplementary Table 2.xlsx", package = "ligandFinder")

ks_peps <- openxlsx::read.xlsx(ks_peps, startRow = 3) %>% as_tibble


ks_peps <- ks_peps %>%
              mutate(tibble_split(Protein, "\\|", names = c("sp", "accession", "gene"))) %>%
              mutate(uniprot_name = id2name[accession])

ks_accessions <- unique(ks_peps[["accession"]])

ks_accessions <- ks_accessions[ks_accessions %in% uniprot_t[["accession"]]]

for(x in ks_accessions) {

  feats <- uniprot_t %>% filter(accession == x) %>% pull(features) %>% `[[`(1)
  sequence_ref <- uniprot_t %>% filter(accession == x) %>% pull(sequence_uni)
  sequence_peps <- ks_peps %>% filter(accession == x) %>% pull(Peptide)
  uni_name_pep <- ks_peps %>% filter(accession == x) %>% pull(uniprot_name) %>% unique
  pep_score <- ks_peps %>% filter(accession == x) %>% pull(PeptideRanker_score)

  new_feats <- bind_rows(feats,
                         lapply(seq_along(sequence_peps), \(y) {

                           alignment <- Biostrings::pairwiseAlignment(sequence_peps[y], sequence_ref, type = "local")

                           tibble(evidence = "sven",
                                  start = Biostrings::start(Biostrings::subject(alignment)),
                                  end = Biostrings::end(Biostrings::subject(alignment)),
                                  description = paste0(uni_name_pep, "x", start, "x", end),
                                  source = "sven",
                                  PeptideRanker_score = pep_score[y]) %>%
                             mutate(type = paste0("sven_pep", y), .before = everything())


                         })
  )

  uniprot_t[uniprot_t[["accession"]] == x, "features"][[1]][[1]] <- new_feats

}




known_peps <- system.file("extdata/GPCRdb_known_pairings_human_plus2more_unique.csv", package = "ligandFinder")

known_peps <- data.table::fread(known_peps) %>% as_tibble

known_peps <- known_peps %>%
  distinct(lig) %>%
  mutate(stringr::str_remove(lig, "^h") %>% tibble_split(., "x", names = c("uniprot_name", "start", "end"))) %>%
  mutate(accession = name2id[uniprot_name]) %>%
  filter(!str_detect(uniprot_name, "^m")) %>%
  mutate(across(all_of(c("start", "end")), as.numeric)) %>%
  mutate(type = paste0("gpcr_pep", row_number()), .by = accession)

known_peps[known_peps[["uniprot_name"]] == "Q9Y494", "accession"] <- "Q9Y494"

kp_accessions <- unique(known_peps[["accession"]])

kp_accessions <- kp_accessions[kp_accessions %in% uniprot_t[["accession"]]]

for(x in kp_accessions) {

    feats <- uniprot_t %>% filter(accession == x) %>% pull(features) %>% `[[`(1)
    pep <- known_peps %>% filter(accession == x)

    new_feats <- bind_rows(feats,
                           tibble(type = pep[["type"]],
                                  evidence = "gpcrdb_gtp",
                                  start = pep[["start"]],
                                  end = pep[["end"]],
                                  description = pep[["lig"]] %>% stringr::str_remove(., "^h"),
                                  source = "gpcrdb_gtp")
    )

    uniprot_t[uniprot_t[["accession"]] == x, "features"][[1]][[1]] <- new_feats

}



####integrate top200NC peptides


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

ligand_list <- ligand_list %>%
  ungroup() %>%
  mutate(type = paste0("top200NC_pep", row_number()), .by = accession) %>%
  mutate(lig = paste0(uniprot_name, "x", start, "x", end))


top200_accessions <- unique(ligand_list[["accession"]])

top200_accessions <- top200_accessions[top200_accessions %in% uniprot_t[["accession"]]]

for(x in top200_accessions) {

  feats <- uniprot_t %>% filter(accession == x) %>% pull(features) %>% `[[`(1)
  pep <- ligand_list %>% filter(accession == x)

  new_feats <- bind_rows(feats,
                         tibble(type = pep[["type"]],
                                evidence = "top200NC",
                                start = pep[["start"]],
                                end = pep[["end"]],
                                description = pep[["lig"]],
                                source = "top200NC")
  )

  uniprot_t[uniprot_t[["accession"]] == x & !is.na(uniprot_t[["accession"]]), "features"][[1]][[1]] <- new_feats

}

uniprot_t %>%
  mutate(test = map_lgl(features, \(x) "top200NC" %in% x[["source"]])) %>%
  pull(test) %>%
  sum


saveRDS(uniprot_t, paste0(s_localDir, "/processed/uniprot_4.rds"))














