






library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"
uniprot_t <- readRDS(paste0(s_localDir, "/processed/uniprot_4.rds"))





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
  filter(source %in% c("gpcrdb_gtp") | (source == "uniprot" & type == "peptide"))

if(nrow(tmp) > 0) {

tmp <- tmp %>%
  pivot_longer(cols = c("start", "end"), values_to = "position") %>%
  group_by(position) %>%
  summarise(terminus = first(name),
            source = paste0(unique(source), collapse = "."),
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

test <- bind_rows(pmap(list(features = secretome$features,
                            sequence_uni = secretome$sequence_uni,
                            gene = secretome$gene),
                       get_known_db_sites))






