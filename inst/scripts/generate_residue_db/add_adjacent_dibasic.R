




library(tidyverse)

files <- c("~/AF2_analysis/all_metrics_oct17.csv")

metrics <- c("iptm", "paeL_mean_in")

dat <- map(files, ~as_tibble(data.table::fread(.))) %>%
  bind_rows

table(dat[["run_dir"]])

dat <- dat %>%
        mutate(run_name = basename(run_dir))

dat <- dat %>%
  filter(run_name == "bm_sep28") %>%
  filter(location == "relevant")


all_ligands <- dat[["p2_name"]] %>% unique(.)

pq_path <- "~/ligandFinder_data/residue_db"

res_db <- arrow::open_dataset(source = pq_path)

residue_data <- res_db %>%
  filter(uni_gene %in% all_ligands) %>%
  select(uni_gene, sequence_uni) %>%
  collect()



ligs <- dat %>%
        distinct(p2_name, p2_range) %>%
        separate_wider_delim(p2_range, delim = "x", names = c("start", "end")) %>%
        mutate(across(all_of(c("start", "end")), as.numeric)) %>%
  {left_join(., residue_data, by = join_by(p2_name == uni_gene))}


find_adjacent_dibasic <- function(start,
                                  end,
                                  sequence_uni) {

  start_seq <- str_sub(sequence_uni, start = max(c(1, start - 5)), end = max(c(1, start - 1)))

  start_gap <- stringr::str_locate(stringi::stri_reverse(start_seq), "KK|KR|RK|RR")[1] - 1

  term <- nchar(sequence_uni)

  end_seq <- str_sub(sequence_uni, start = min(c(term, end + 1)), end = min(c(term, end + 5)))

  end_gap <- stringr::str_locate(end_seq, "KK|KR|RK|RR")[1] - 1

  tibble(N_term_db = start_gap,
         C_term_db = end_gap)
}

ligs <- ligs %>%
        mutate(dibasic = pmap(.l = list(start, end, sequence_uni), .f = find_adjacent_dibasic)) %>%
        mutate(bind_rows(dibasic)) %>%
        mutate(p2_range = stringr::str_c(start, end, sep = "x")) %>%
        select(-dibasic, -sequence_uni, -start, -end)


dat <- dat %>%
  {left_join(., ligs, by = join_by(p2_name, p2_range))}


known_dat <- dat %>%
            filter(known_pair == "known" & iptm > 0.5 & !is.na(N_term_db) & !is.na(C_term_db))


















