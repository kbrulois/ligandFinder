


uniprot_t <- dplyr::bind_rows(
  map(list.files(uni_dir), ~readRDS(paste0(uni_dir, "/", .)))
)

to_save <- uniprot_t %>%
  select(-annotations, -full_name, -sequence) %>%
  unnest(features)

data.table::fwrite(to_save, "~/Desktop/uniprot_features.csv")

to_save <- uniprot_t %>%
  select(-features, -full_name, -sequence) %>%
  unnest(annotations)

data.table::fwrite(to_save, "~/Desktop/uniprot_annotations.csv")


uniprot_t %>%
  filter(gene %in% c("GPR15", "MARCHF4")) %>%
  select(-annotations, -full_name, -sequence) %>%
  unnest(features) %>%
  View

seqs <- uniprot_t %>%
  filter(gene %in% c("MARCHF4")) %>%
  pull(sequence)

uniprot_t <- uniprot_t %>%
  mutate(has_TM = map_lgl(features, \(x) {"transmembrane region" %in% x[["type"]]})) %>%
  mutate(has_topo = map_lgl(features, \(x) {"topological domain" %in% x[["type"]]})) %>%
  mutate(has_complete_top = map_lgl(features, \(x) {sum(x[["type"]] == "topological domain") == sum(x[["type"]] == "transmembrane region") + 1}))


uniprot_t <- uniprot_t %>%
          mutate(uniprot_location = map_chr(annotations, \(x) { x %>%
                                                        filter(annotation_type == "subcellular location") %>%
                                                        pull(annotation) %>%
                                                        paste0(., collapse = ";")}))

hpa_location <- data.table::fread("~/Downloads/proteinatlas.tsv") %>%
  as_tibble

uniprot_t <- left_join(uniprot_t, hpa_location %>%
            select(c(1, 61, 62,63, 74,75)), by = join_by(gene == Gene))

uniprot_t <- uniprot_t %>%
  mutate(in_pep_csv = if_else(gene %in% unique(pep_csv$gene), "yes", "no"))




data.table::fwrite(uniprot_t %>% select(!where(is.list)), "~/Desktop/uniprot_stuff.csv")

table(uniprot_t$has_TM, uniprot_t$has_topo)
table(uniprot_t$has_topo)
table(uniprot_t$has_TM)
table(uniprot_t$has_TM, uniprot_t$has_complete_top)

c(s = "signal peptide",
  i = "intracellular",
  t = "TM",
  e = "extracellular",
  l = "luminal",
  `-` = "no data available")

topolize <- function(seq,
                     features) {
  if("type" %in% colnames(features)) {
  feats_sub <- features %>%
      filter(type %in% c("transmembrane region", "topological domain", "signal peptide")) %>%
      mutate(final = case_when(type == "transmembrane region" ~ "t",
                               type == "topological domain" & description %in% c("Extracellular") ~ "e",
                               type == "topological domain" & description %in% c("Luminal") ~ "L",
                               type == "topological domain" & description == "Cytoplasmic" ~ "i",
                               type == "signal peptide" ~ "s",
                               TRUE ~ "-"))

    if(nrow(feats_sub) > 0) {
    topo <- rep("-", nchar(seq))

    for(i in 1:nrow(feats_sub)) {
      strt <- feats_sub %>% slice(i) %>% pull(start)
      end <- feats_sub %>% slice(i) %>% pull(end)
      if(is.na(strt) | is.na(end)) {next}
      topo[strt:end] <- feats_sub %>% slice(i) %>% pull(final)
    }

    topo <- stringr::str_c(topo, collapse = "")
    } else {
      topo <- NA
    }
  } else {
    topo <- NA
  }

  return(topo)

}

uniprot_t <- uniprot_t %>%
            mutate(topo = map2_chr(.x = sequence, .y = features,
                                           .f = topolize))



pep_csv <- data.table::fread("~/Desktop/peptides_combined_overlap6_AF2_annotated_UNCLEAR 2.txt") %>% as_tibble

pep_csv %>%
  nest_by(gene)

pep_csv <- left_join(pep_csv, uniprot_t %>% distinct(gene, .keep_all = TRUE), by = "gene")

pep_csv <- pep_csv %>%
  rowwise %>%
  mutate(topo_location = case_when(overlap_descriptions %in% c("Extracellular") ~ 'e_topo',
                                   overlap_descriptions %in% c("Luminal") ~ 'L_topo',
                                    grepl("Extracellular", overlap_descriptions) ~ 'e_topo_partial',
                                    grepl("Cytoplasmic", overlap_descriptions) ~ 'i_topo')) %>%
  mutate(prot_location = case_when(grepl("Secreted", uniprot_location) ~ 'e_location',
                                   grepl("Secreted", `Secretome location`) ~ 'e_location',
                                   TRUE ~ "i_location")) %>%
  mutate(location_final = case_when(topo_location == "e_topo" ~ "e",
                                    prot_location == "e_location" & topo_location %in% c("e_topo_partial", "i_topo") ~ "unclear",
                                    prot_location == "e_location" ~ "e",
                                    prot_location == "i_location" | topo_location == "i_topo" ~ "i",
                                    TRUE ~ "unclear"))

pep_csv <- pep_csv %>%
  mutate(peptide_topo = stringr::str_sub(topo, start = start_pos, end = end_pos))

data.table::fwrite(pep_csv %>% select(!where(is.list)), "~/Desktop/peptides_combined_overlap6_AF2_annotated_kb.csv")




res <- httr::POST("https://biolib.com/api/v1/deeptmhmm/run",
            body = list(sequence = seqs), encode = "json")

httr::content(res)


# install once:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("phobius")
# install.packages(c("httr","jsonlite","stringr","dplyr"))

library(httr); library(jsonlite); library(stringr); library(dplyr); library(phobius)

uniprot_id <- "P07550"  # replace with your protein (e.g., MARCH4 ortholog)

fa <- GET(paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta"))
stop_for_status(fa)
seq <- sub("^>.*\\n", "", content(fa, "text"))
seq <- gsub("\\s", "", seq)

# 2) Run Phobius (predicts TM + inside/outside + signal peptide)
pb <- phobius::phobius(seqs)
pb

# 3) Parse topology string (e.g., "iiiiMMMMMooooMMMMiiiii...")
topo <- pb$topology[[1]]            # per-residue states: i=cytosol, o=non-cytosol, M=TM
runs <- rle(strsplit(topo, "")[[1]])
ends <- cumsum(runs$lengths)
starts <- c(1, head(ends, -1) + 1)

topo_df <- tibble(
  start = starts,
  end   = ends,
  state = runs$values
) %>%
  mutate(region = case_when(
    state == "M" ~ "Transmembrane helix",
    state == "i" ~ "Cytosolic",
    state == "o" ~ "Non-cytosolic (luminal/extracellular)",
    TRUE ~ "Other"
  ))

# 4) Optional sanity checks
# Positive-inside rule: count K/R in loops
kr_count <- function(s) str_count(s, "[KR]")
loops <- topo_df %>% filter(state != "M") %>%
  rowwise() %>%
  mutate(KR = kr_count(substr(seq, start, end))) %>%
  ungroup()

list(
  topology = topo_df,
  loops_with_KR = loops
)


