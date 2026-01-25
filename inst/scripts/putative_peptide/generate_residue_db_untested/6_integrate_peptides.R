

library(tidyverse)

s_localDir <- "~/peptide_alg/build_residue_db"

uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_4.rds"))

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
  ) %>% distinct(.keep_all = TRUE)

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
                                  source = "known")
    )

    uniprot_t[uniprot_t[["accession"]] == x, "features"][[1]][[1]] <- new_feats

}




saveRDS(uniprot_t, paste0(s_localDir, "/processed/uniprot_4.rds"))














uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_4.rds"))







to_expand <- tibble(type = c("gpcr_pep", "sven_pep", "peptide", "sequence variant", "signal peptide"),
                    source = c("known", "sven", "uniprot", "uniprot", "uniprot"),
                    filter_by = c("source", "source", "type_source", "type_source", "type_source"))

features <- uniprot_t %>% filter(accession == "Q8N5G0") %>% pull(features) %>% `[[`(1)

sequence_uni <- uniprot_t %>% filter(accession == "Q8N5G0") %>% pull(sequence_uni)

fill_type = c("type", "desc")



expand_features <- function(features, sequence_uni) {

  feat_clean <- features %>%
                  filter(!is.na(start)) %>%
                  filter(start <= dat_size) %>%
                  filter(is.na(end) | end <= dat_size)

  dat_size <- nchar(sequence_uni)

  dat_ts <- dat_s <- tibble(.rows = dat_size)

  to_expand_s <- to_expand %>%
                  filter(filter_by == "source")

  if(nrow(to_expand_s) > 0) {

    cols <- to_expand_s %>% pull(source) %>% unique()

    col_names <- expand.grid(cols = cols,
                             fill_type = fill_type) %>%
                      mutate(cols = paste0(cols, "_", fill_type)) %>%
                      pull(cols)

    dat_s[col_names] <- NA

    for(x in cols) {

     feat_sub <- feat_clean %>%
                  filter(source == x)

      for(i in seq_len(nrow(feat_sub))) {

        start_val <- feat_sub[["start"]][i]
        end_val <- feat_sub[["end"]][i]

        if(is.na(end_val)) {
          dat_s[[paste0(x, "_type")]][start_val] <- feat_sub[["type"]][i]
          dat_s[[paste0(x, "_desc")]][start_val] <- feat_sub[["description"]][i]
        } else {
          dat_s[[paste0(x, "_type")]][start_val:end_val] <- feat_sub[["type"]][i]
          dat_s[[paste0(x, "_desc")]][start_val:end_val] <- feat_sub[["description"]][i]
        }
      }

    }

    dat_s <- dat_s %>%
      select(where(~ !all(is.na(.x))))
  }

  #####different type of expansion
  to_expand_ts <- to_expand %>%
    filter(filter_by == "type_source")

  if(nrow(to_expand_ts) > 0) {

    to_expand_ts <- to_expand_ts %>%
                      distinct(source, type)

    cols <- to_expand_ts %>%
              mutate(type2 = paste0(source, "_", type)) %>%
              pull(type2)

    col_names <- expand.grid(cols = cols,
                             fill_type = fill_type) %>%
      mutate(cols = paste0(cols, "_", fill_type)) %>%
      pull(cols)

    dat_ts[col_names] <- NA

    for(x in seq_along(cols)) {

      type_val <- to_expand_ts %>% slice(x) %>% pull(type)
      source_val <- to_expand_ts %>% slice(x) %>% pull(source)
      x <- cols[x]

      feat_sub <- feat_clean %>%
                    filter(type == type_val & source == source_val)

      for(i in seq_len(nrow(feat_sub))) {

        start_val <- feat_sub[["start"]][i]
        end_val <- feat_sub[["end"]][i]

        if(is.na(end_val)) {
          dat_ts[[paste0(x, "_type")]][start_val] <- feat_sub[["type"]][i]
          dat_ts[[paste0(x, "_desc")]][start_val] <- feat_sub[["description"]][i]
        } else {
          dat_ts[[paste0(x, "_type")]][start_val:end_val] <- feat_sub[["type"]][i]
          dat_ts[[paste0(x, "_desc")]][start_val:end_val] <- feat_sub[["description"]][i]
        }
      }

    }

    dat_ts <- dat_ts %>%
      select(where(~ !all(is.na(.x))))
  }

  return(bind_cols(dat_s, dat_ts))

}



uniprot_t <- uniprot_t %>%
  mutate(features_expanded = map2(features, sequence_uni, expand_features))


bind_cols_unique_names <- function(...) {
  dfs <- list(...)

  # Convert each to tibble; replace NULLs with 0-row tibble
  dfs <- lapply(dfs, function(x) {
    if (is.null(x)) return(tibble())
    if (is.vector(x) && !is.list(x)) return(tibble(value = x))
    tibble::as_tibble(x)
  })

  out <- dfs[[1]]

  for (i in 2:length(dfs)) {
    df <- dfs[[i]]

    # Identify cols to keep
    keep <- !names(df) %in% names(out)

    # If df has 0 columns after filtering, skip
    if (!any(keep)) next

    df2 <- df[, keep, drop = FALSE]

    out <- dplyr::bind_cols(out, df2)
  }

  out
}

uniprot_t <- uniprot_t %>%
  mutate(to_expand = pmap(list(features_expanded, cons, alignment_AA, dssp, af_missense), bind_cols_unique_names))


uniprot_t_e <- uniprot_t %>%
              select(-where(is.list), to_expand) %>%
              unnest(to_expand)

uniprot_t_e <- uniprot_t_e %>%
  filter(is.na(`uniprot_signal peptide_type`)) %>%
  filter(!is.na(cons_og)) %>%
  mutate(known = if_else(is.na(known_type), "other", "ligand"))

cols_to_sum <- c("cons_rs", "cons_lrs", "blos_wt_all",
                 "blos_uw_all", "blos_wt_mam", "blos_uw_mam", "gran_wt_all", "gran_uw_all",
                 "gran_wt_mam", "gran_uw_mam", "blos_nr_all", "blos_nr_mam", "gran_nr_all",
                 "gran_nr_mam")

cols_to_sum2 <- c(cols_to_sum, paste0(cols_to_sum, "_n"))

mean_dif <- uniprot_t_e %>%
  filter(relASA > 0.5 & SS != "E") %>%
  mutate(across(all_of(cols_to_sum), ~ .x * species_limit, .names = "{.col}_n")) %>%
  group_by(known) %>%
  summarise(across(all_of(cols_to_sum2), ~mean(., na.rm = TRUE))) %>%
  summarise(across(all_of(cols_to_sum2), ~ (.x[1] - .x[2])))

sd_res <- uniprot_t_e %>%
  filter(relASA > 0.5 & SS != "E") %>%
  mutate(across(all_of(cols_to_sum), ~ .x * species_limit, .names = "{.col}_n")) %>%
  summarise(across(all_of(cols_to_sum2), ~sd(., na.rm = TRUE)))


svglite::svglite(filename = "~/AF2_analysis/conservation_comparison.svg", width = 6, height = 6)
ggplot2::ggplot((mean_dif/sd_res) %>%
                  pivot_longer(everything()) %>%
                  mutate(col = if_else(stringr::str_detect(name, "_n$"), "myo_normalized", "not_normalized"))) +
  geom_col(aes(x = value, y = name, fill = col), color = "black", width = 0.8) +
  xlab("Cohen's D") +
  ylab("") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  ggtitle("residues of known ligands versus other residues\n(gated on relASA > 0.5 and non-strand non-sp)") +
  theme_bw()
dev.off()



ggplot2::ggplot(uniprot_t_e %>% filter(relASA > 0.5)) +
  ggplot2::geom_violin(mapping = ggplot2::aes(x = known, y = blos_nr_all))



uniprot_t_e %>%
  nest(data = -known) %>%
  summarise(
    tidy(t.test(value ~ group, data = df))
  )





