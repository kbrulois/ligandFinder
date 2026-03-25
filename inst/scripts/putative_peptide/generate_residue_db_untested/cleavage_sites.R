






library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"
uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_4.rds"))



library(reticulate)

reticulate::py_install("iPython")

source_python(file = "~/Python_scripts/cap_peptides.py")





#testing examples
#features <- uniprot_t %>% filter(gene == "GCG") %>% pull(features) %>% `[[`(1)
#sequence_uni <- uniprot_t %>% filter(gene == "GCG") %>% pull(sequence_uni)
#features <- uniprot_t$features[[442]]
#sequence_uni <- uniprot_t$sequence_uni[[442]]



get_known_db_sites <- function(features, sequence_uni, gene, wN = 10, wC = 10) {


shift <- max(c(wN, wC))

n_trunc <- features %>%
                filter(type == "signal peptide") %>%
                pull(end)

if(length(n_trunc) == 0) {
  n_trunc <- 1
}
if(is.na(n_trunc)) {
  n_trunc <- 1
}
if(length(n_trunc) > 1) {
  n_trunc <- max(n_trunc, na.rm = TRUE)
}

c_trunc <- nchar(sequence_uni)

AA_seq <- stringr::str_sub(string = sequence_uni, start = n_trunc + 1, end = c_trunc) %>%
  stringr::str_c(stringr::str_c(rep("-", n_trunc + shift), collapse = ""), ., collapse = "") %>%
  stringr::str_c(., stringr::str_c(rep("-", shift), collapse = ""))

tmp <- features %>%
  filter(source %in% c("gpcrdb_gtp", "uni_pep"))

if(nrow(tmp) > 0) {

tmp <- tmp %>%
  mutate(pep_id = paste0(gene, "_", start, "x", end)) %>%
  pivot_longer(cols = c("start", "end"), values_to = "position") %>%
  group_by(position) %>%
  summarise(terminus = first(name),
            source = paste0(unique(source), collapse = "."),
            pep_id = first(pep_id),
            .groups = "drop")

tmp <- tmp %>%
  mutate(case_when(terminus == "start" ~ tibble(start = position - (wN - shift - 1),
                                                end = position + (shift + wC)),
                   terminus == "end" ~ tibble(start = position - (wN - shift - 1),
                                              end = position + (shift + wC)))) %>%
  mutate(seq_cleavage = map2_chr(start, end, \(x, y) {stringr::str_sub(AA_seq, start = x, end = y)})) %>%
  filter(seq_cleavage != "-----") %>%
  mutate(has_db = map_lgl(seq_cleavage, ~stringr::str_detect(., "KR|RK|RR|KK"))) %>%
  select(-start, -end) %>%
  mutate(gene = gene) %>%
  mutate(db_loc = stringr::str_locate_all(seq_cleavage, "KR|RK|RR|KK"))


center_db <- function(db_loc) {

  if(nrow(db_loc) == 0) {
    return(NA)
  } else {

  return(db_loc %>%
          as_tibble %>%
          mutate(offset = min(abs(((shift / 2) + 0.5) - start), abs(((shift / 2) + 0.5) - end))) %>%
          filter(offset == min(offset)) %>%
          mutate(nudge = start - shift) %>%
          slice(1) %>%
          pull(nudge)
  )

  }
}

tmp <- tmp %>%
  mutate(nudge = map_dbl(db_loc, center_db)) %>%
  mutate(case_when(terminus == "start" ~ tibble(start = position - (wN - shift - 1) + nudge,
                                                end = position + (shift + wC) + nudge),
                   terminus == "end" ~ tibble(start = position - (wN - shift - 1) + nudge,
                                              end = position + (shift + wC) + nudge))) %>%
  mutate(seq_dbcenter = map2_chr(start, end, \(x, y) {stringr::str_sub(AA_seq, start = x, end = y)}))

} else {
  tmp <- tibble(position = numeric(),
                terminus = character(),
                source = character(),
                seq_cleavage = character(),
                has_db = logical(),
                gene = character())
}

return(tmp)

}

known_peps <- bind_rows(pmap(list(features = secretome$features,
                            sequence_uni = secretome$sequence_uni,
                            gene = secretome$gene),
                       get_known_db_sites))




bm_to_join <- dat %>%
  mutate(gene = setNames(id_map$`Gene Names (primary)`, id_map$Entry)[p2_id]) %>%
  mutate(pep_id = paste0(gene, "_", p2_range)) %>%
  filter(known_pair == "known") %>%
  filter(location == "relevant") %>%
  select(pep_id, lig1_end) %>%
  group_by(pep_id) %>%
  summarise(across(everything(), first))

test <- left_join(known_peps,
                  bm_to_join, by = "pep_id")


test <- test %>%
          mutate(lig1_end_clean = case_when(lig1_end %in% c("1C", "2C") ~ "C",
                                            lig1_end %in% c("1N", "2N") ~ "N",
                                            TRUE ~ "loop"))


strReverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

test <- test %>%
  filter(!map_lgl(seq_cleavage, ~stringr::str_detect(., "^----------"))) %>%
  filter(!map_lgl(seq_cleavage, ~stringr::str_detect(., "----------$"))) %>%
  mutate(cap_seq = if_else(terminus == "start",
                           stringr::str_sub(seq_cleavage, 11, 20),
                           stringr::str_sub(seq_cleavage, 1, 10))) %>%
  rowwise() %>%
  mutate(cap_seq = if(terminus == "end") {strReverse(cap_seq)} else {cap_seq}) %>%
  mutate(bind_rows(map(cap_seq, \(x) {
    compute_cap_features(x)})))




ggplot2::ggplot(data = test, aes(x = terminus, y = cap_score, color = terminus)) +
  ggplot2::geom_violin(
  scale = "width",
  #position = ggplot2::position_dodge(0.9),
  trim = T,
  linewidth = 0.1,
  width = 0.4,
  alpha = 0.5,
  show.legend = F) +

  ggbeeswarm::geom_quasirandom(
                               size = 0.5,
                               width = 0.2,
                               stroke = 0.2,
                               inherit.aes = TRUE,
                               show.legend = TRUE) +

  ggplot2::geom_boxplot(width=0.1,
                        fill = "white",
                        #position = ggplot2::position_dodge(0.9),
                        outlier.shape = NA,
                        coef = 0,inherit.aes = TRUE,
                        show.legend = FALSE) +
  ggplot2::facet_grid(cols = vars(lig1_end_clean)) +

  theme_bw()










data.table::fwrite(test, "~/Desktop/cap_scores.csv")










