

uniprotid_gpcrdb <- data.table::fread("https://files.gpcrdb.org/uniprot_mapping.txt", header = FALSE) %>%
                    as_tibble

gpcrdb_receptors <- jsonlite::read_json("https://gpcrdb.org/services/receptorlist/")

gpcrdb_receptors <- bind_rows(lapply(gpcrdb_receptors, \(x) {
    x %>%
    as_tibble %>%
    mutate(endogenous_ligands = map_chr(endogenous_ligands, \(y) {paste(unlist(y), collapse = "; ")}))
}))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

gpcrdb_receptors <- left_join(x = gpcrdb_receptors, y = id_mapping %>% dplyr::rename(accession = Entry), by = "accession")

gpcrdb_receptors <- gpcrdb_receptors %>%
  group_by(accession) %>%
  summarise(across(everything(), get_unique_or_collapse)) %>%
  filter(species == "Homo sapiens")

gpcrdb_receptors <- gpcrdb_receptors %>%
  dplyr::rename(uniprot_id = accession,
         uniprot_name = `Entry Name`,
         gene_name_primary = `Gene Names (primary)`,
         gene_name_all = `Gene Names`) %>%
  select(-species, -Reviewed, -Organism) %>%
  rename_with(~paste0("gpcrdb: ", .)) %>%
  mutate(source = "gpcrdb")


gtp <- data.table::fread("https://www.guidetopharmacology.org/DATA/GPCRTargets.csv") %>%
  as_tibble

gtp <- gtp %>%
  dplyr::rename(gene_name = `HGNC symbol`)

gtp <- gtp %>%
  select(-c(18:33)) %>%
  select(-`Target systematic name`, -`Target abbreviated name`, -`Human genetic localisation`) %>%
  rename_with(~paste0("gtp: ", .)) %>%
  dplyr::rename(uniprot_id = `gtp: Human SwissProt`) %>%
  select(uniprot_id, `gtp: gene_name`, `gtp: Family name`, everything()) %>%
  mutate(source = "gtp", .after = `gtp: gene_name`)

gtp <- gtp %>%
        filter(uniprot_id != "")

gtp2 <- left_join(gtp %>% distinct(uniprot_id, .keep_all = TRUE),
                  id_map %>% dplyr::rename(uniprot_id = Entry,
                                           uniprot_name = `Entry Name`,
                                           gene_name_primary = `Gene Names (primary)`,
                                           gene_names_all = `Gene Names`),
                  by = "uniprot_id")













czList <- data.table::fread(system.file("extdata/Catherines_gpcrs_List_ballesteros_numbering.txt", package = "ligandFinder")) %>%
  as_tibble

ecbList <- data.table::fread(system.file("extdata/GPCR list draft March 2025 for APACE ecb orphan.txt", package = "ligandFinder")) %>%
  as_tibble

ecbList <- left_join(ecbList %>% dplyr::rename(gene = Gene_human), czList, by = "gene")

ecb_gene_sym_to_uni <- id_map %>%
  mutate(`Gene Names (primary)` = stringr::str_extract(pattern = "^([^;]+)", string = `Gene Names (primary)`)) %>%
  group_by(`Gene Names (primary)`) %>%
  filter(n() == 1 | (n() > 1 & Reviewed == "reviewed")) %>%
  filter(`Gene Names (primary)` %in% ecbList$gene)

ecb_gene_sym_to_uni %>% filter(duplicated(`Gene Names (primary)`)) %>% pull(`Gene Names (primary)`) -> dup_genes

ecb_gene_sym_to_uni %>%
  filter(`Gene Names (primary)` %in% dup_genes)

ecb_gene_sym_to_uni2 <- ecb_gene_sym_to_uni %>%
  distinct(`Gene Names (primary)`, .keep_all = TRUE)



problem_genes <- ecbList$gene[!ecbList$gene %in% ecb_gene_sym_to_uni2$`Gene Names (primary)`]
problem_genes <- tibble(problem_genes = problem_genes)
problem_genes$primary_gene <- NA

problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[1]] <- "OR2A1"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[2]] <- "OR2I1"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[3]] <- "OR3A4P"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[4]] <- "OR4F3"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[5]] <- "OR5BS1"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[6]] <- "OR51C1"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[7]] <- "OR56B2"
problem_genes$primary_gene[problem_genes$problem_genes == problem_genes$problem_genes[8]] <- "OR8G3"

ecbList$gene_change_tracker <- NA
for(problem_gene in problem_genes$problem_genes) {
  old_name <- problem_genes %>% filter(problem_genes == problem_gene) %>% pull(problem_genes)
  new_name <- problem_genes %>% filter(problem_genes == problem_gene) %>% pull(primary_gene)
  ecbList[ecbList$gene == problem_gene, "gene_change_tracker"] <- paste(old_name, "->", new_name)
  ecbList[ecbList$gene == problem_gene, "gene"] <- new_name
}

ecbList <- ecbList %>%
  group_by(gene) %>%
  filter(n() == 1 | (n() > 1 & !is.na(gene_change_tracker))) %>%
  dplyr::rename(`gene: old_name -> new_name` = gene_change_tracker) %>%
  distinct(gene, .keep_all = TRUE)

ecbList <- ecbList %>%
  dplyr::rename(`uniprot_id_old` = `uniprot_id`,
                `ballesteros_numbering_old` = ballesteros_numbering) %>%
  mutate(uniprot_id_old = toupper(uniprot_id_old))

ecbList_c <- left_join(ecbList %>% ungroup, ecb_gene_sym_to_uni %>% dplyr::rename(gene = `Gene Names (primary)`), by = "gene")

ecbList_c <- ecbList_c %>%
  dplyr::rename(uniprot_name = `Entry Name`,
                uniprot_id = Entry,
                gene_name_primary = gene,
                gene_names_all = `Gene Names`) %>%
  select(-Organism) %>%
  select(uniprot_name,
         `Order of runs (priority)`,
         `Exclude due to N term >160AA`,
         `Prioritization Notes`,
         `Class or type`,
         gene_name_primary,
         `Protein names`,
         gene_names_all,
         `Length`,
         uniprot_id,
         Reviewed,
         uniprot_id_old,
         `gene: old_name -> new_name`,
         ballesteros_numbering_old
  )

ecbList_c <- ecbList_c %>%
  dplyr::rename(uniprot_name_old = uniprot_id_old) %>%
  mutate(`uniprot_name: old_name -> new_name` = if_else(uniprot_name != uniprot_name_old, paste(uniprot_name_old, "->", uniprot_name), NA)) %>%
  mutate(source = "ecb")











comb_gpcr <- full_join(ecbList_c %>% rename_with(~paste0("ecb: ", .), .cols = -c(uniprot_id,
                                                                            uniprot_name,
                                                                            Reviewed,
                                                                            gene_names_all,
                                                                            gene_name_primary,
                                                                            `Protein names`,
                                                                            source)),
                  gtp2, by = "uniprot_id")

cols_to_combine <- sub(".x$", "", grep(".x$", colnames(comb_gpcr), value = TRUE))


for(x in cols_to_combine) {
  comb_gpcr <- comb_gpcr %>%
              mutate(!!sym(x) := get_unique_or_collapse_v(!!sym(paste0(x, ".x")), !!sym(paste0(x, ".y")))) %>%
              select(-all_of(c(paste0(x, ".x"), paste0(x, ".y"))))
}




comb_gpcr <- full_join(comb_gpcr, gpcrdb_receptors %>%
                              dplyr::rename(`gpcrdb: gene_names_all` = `gpcrdb: gene_name_all`) %>%
                              rename_with(~sub("gpcrdb: ", "", .),
                                          .cols = all_of(c(paste("gpcrdb:", c("uniprot_id",
                                                                            "uniprot_name",
                                                                            "gene_names_all",
                                                                            "gene_name_primary",
                                                                            "Protein names")),
                                                                            "source"))),
                   by = "uniprot_id")

cols_to_combine <- sub(".x$", "", grep(".x$", colnames(comb_gpcr), value = TRUE))

for(x in cols_to_combine) {
  comb_gpcr <- comb_gpcr %>%
    mutate(!!sym(x) := get_unique_or_collapse_v(!!sym(paste0(x, ".x")), !!sym(paste0(x, ".y")))) %>%
    select(-all_of(c(paste0(x, ".x"), paste0(x, ".y"))))
}


comb_gpcr <- comb_gpcr %>%
  select(all_of(c("uniprot_name",
                  "gene_name_primary",
                  "Protein names",
                  "source")), everything())

comb_gpcr <- left_join(comb_gpcr, uniprotid_gpcrdb %>% dplyr::rename(uniprot_id = V1,
                                             gpcrdb_id = V2), by = "uniprot_id") %>%
  relocate(gpcrdb_id, .after = uniprot_name) %>%
  mutate(uniprot_gpdcrdb_id_match = if_else(uniprot_name == toupper(gpcrdb_id), "same", NA), .after = "gpcrdb_id")


bw_res_un <- map(comb_gpcr %>%
                   pull(uniprot_name),
                 ~getBallesterosFromGPCRDB(gene_name = tolower(.)), .progress = TRUE)

bw_res_gpcrdb <- map(comb_gpcr %>% pull(gpcrdb_id),
                     ~getBallesterosFromGPCRDB(gene_name = tolower(.)), .progress = TRUE)

comb_gpcr <- comb_gpcr %>%
  mutate(`bw: full_table` = if_else(map_lgl(bw_res_un, ~nrow(.) > 0), bw_res_un, bw_res_gpcrdb), .after = source)


comb_gpcr <- comb_gpcr %>%
  mutate(map_df(`bw: full_table`, .f = get_bw_inds), .after = `bw: full_table`) %>%
  mutate(map_df(`bw: full_table`, .f = get_lengths), .after = `bw: full_table`) %>%
  mutate(`bw: Long_nterm?` = if_else(`bw: length N-term` > 160, "long", NA), .after = `bw: full_table`)

comb_gpcr %>%
  filter(is.na(`ecb: Order of runs (priority)`)) %>%
  select(`ecb: Order of runs (priority)`, uniprot_name, `gtp: Family name`, `gpcrdb: receptor_class`, `gpcrdb: receptor_family`)

comb_gpcr$`ecb: Order of runs (priority)`[833:841] <- "#3"

comb_gpcr <- comb_gpcr %>%
                dplyr::slice(-842)


res_db <- arrow::open_dataset(source = pq_path)

residue_anno <- res_db %>%
  filter(uni_gene %in% comb_gpcr[["uniprot_name"]]) %>%
  group_by(gene_grp) %>%
  collect() %>%
  ungroup()

residue_anno

comb_gpcr <- left_join(comb_gpcr, residue_anno %>%
                        dplyr::rename(uniprot_name = uni_gene) %>%
                        select(uniprot_name, features),
                      "uniprot_name")

thing <- c("OR2I1", "OR5BS1", "OR51C1", "OR56B2", "OR8G3")

comb_gpcr[comb_gpcr$gene_name_primary %in% thing, "uniprot_name"][[1]] <- c("OR2I1", "O5BS1", "O51C1", "O56B2", "OR8G3")

fasta_files <- c(system.file("extdata/hGPCRs_classA.fasta", package = "ligandFinder"),
                 system.file("extdata/2024_07_UniProt_GtoP_peptide_and_orphan_GPCRs.fasta", package = "ligandFinder"))

fasta_files <- c('~/R_projects/ligandFinder/inst/extdata/hGPCRs_classA.fasta',
                 '~/R_projects/ligandFinder/inst/extdata/2024_07_UniProt_GtoP_peptide_and_orphan_GPCRs.fasta')

a3m <- readLines(fasta_files[1], skipNul = TRUE)

a3m_len <- length(a3m)

hGPCRs_classA <- tibble(id = a3m[seq(1, a3m_len, by = 2)] %>% sub("^>", "", .),
                        sequence = a3m[seq(2, a3m_len, by = 2)] %>% gsub("[a-z]", "", .)) %>%
  mutate(stringr::str_split(id, " ", simplify = TRUE) %>%
           as_tibble %>%
           rename_with(.fn = function(x) {c("uniprot_name", "uniprot_id", "class", "TM")})) %>%
  mutate(uniprot_name = stringr::str_remove(uniprot_name, "_HUMAN$")) %>%
  mutate(uniprot_id = stringr::str_remove(uniprot_id, "^AC=")) %>%
  mutate(class = stringr::str_remove(class, "class=")) %>%
  mutate(TM = stringr::str_remove(TM, "TM="))

a3m <- readLines(fasta_files[2], skipNul = TRUE)

a3m_len <- length(a3m)

GtoP_peptide_and_orphan_GPCRs <- tibble(uniprot_name = a3m[seq(1, a3m_len, by = 2)] %>% sub("^>", "", .),
                                        sequence = a3m[seq(2, a3m_len, by = 2)] %>% gsub("[a-z]", "", .))



hGPCRs_classA <- hGPCRs_classA %>%
  mutate(GtoP_peptide_and_orphan_GPCR = if_else(uniprot_name %in% GtoP_peptide_and_orphan_GPCRs$uniprot_name, "GtoP_sub_list", ""))

hGPCRs_classA <- hGPCRs_classA %>%
  rename_with(~paste0("UCSD: ", .)) %>%
  dplyr::rename(uniprot_name = `UCSD: uniprot_name`)

hGPCRs_classA$uniprot_name[!hGPCRs_classA$uniprot_name %in% gpcr_list$uniprot_name]

unname(comb_gpcr$uniprot_name[!gpcr_list$uniprot_name %in% hGPCRs_classA$uniprot_name])

hGPCRs_classA <- hGPCRs_classA %>%
  mutate(source = "UCSD")

gpcr_list <- full_join(comb_gpcr, hGPCRs_classA, by = "uniprot_name")

gpcr_list <- gpcr_list %>%
  mutate(across(starts_with("source."), ~replace_na(., ""))) %>%
  mutate(source = paste(source.x, source.y)) %>%
  mutate(source2 = if_else(!is.na(source.x) & ! source.x == "", "SU", "")) %>%
  mutate(source2 = paste(source2, source.y))










gpcr_list <- gpcr_list %>%
  mutate(features = map(features, clean_feats)) %>%
  mutate(signal_peptide = map_int(features, get_sp)) %>%
  mutate(last_AA = map_int(features, get_lastAA)) %>%
  mutate(model_id = if_else(is.na(signal_peptide), uniprot_id, paste0(uniprot_id, ",", signal_peptide, "-", last_AA)), .after = "uniprot_name") %>%
  mutate(model_name = if_else(is.na(signal_peptide), uniprot_name, paste0(uniprot_name, ",", signal_peptide, "-", last_AA)), .after = "uniprot_name") %>%
  mutate(dNT_index = `bw: indices 1.50` - 32) %>%
  mutate(dCT_index = `bw: indices 7.50` + 24) %>%
  mutate(first_AA = if_else(is.na(signal_peptide), 1, signal_peptide))

to_iterate <- expand.grid(c("uniprot_name", "uniprot_id"), c("dNT", "dCT"), stringsAsFactors = FALSE)

for(i in 1:nrow(to_iterate)) {
  uniprot_identifier <- to_iterate[[1]][i]
  trunc <- to_iterate[[2]][i]
  if(trunc == "dNT") {
gpcr_list <- gpcr_list %>%
         mutate(!!paste0("model_",
                stringr::str_remove(uniprot_identifier, "uniprot_"),
                "_",
                trunc) := case_when(!!sym(paste0(trunc, "_index")) > first_AA ~ paste0(!!sym(uniprot_identifier), ",", !!sym(paste0(trunc, "_index")), "-", last_AA),
                             signal_peptide > 1 ~ paste0(!!sym(uniprot_identifier), ",", signal_peptide, "-", last_AA),
                             TRUE ~ !!sym(uniprot_identifier)))
  } else {
gpcr_list <- gpcr_list %>%
      mutate(!!paste0("model_",
                      stringr::str_remove(uniprot_identifier, "uniprot_"),
                      "_",
                      trunc) := case_when(!!sym(paste0(trunc, "_index")) < last_AA ~ paste0(!!sym(uniprot_identifier), ",", first_AA, "-", !!sym(paste0(trunc, "_index"))),
                                          TRUE ~ !!sym(uniprot_identifier)))


  }

}


gpcr_list <- gpcr_list %>%
  mutate(`bw: full_table` = map(`bw: full_table`, extend_bw_notation))

saveRDS(gpcr_list, "~/R_projects/ligandFinder/inst/extdata/gpcr_list.rds")







