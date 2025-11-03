




library(tidyverse)

data.table::fread("~/AF2_analysis/brinp_v5.csv") %>%
  as_tibble %>%
  #mutate(run_name = "bm_sep28") %>%
  mutate(depth = if_else(is.na(depth), depth_percent, depth)) %>%
  mutate(radius = if_else(is.na(radius), radius_percent, radius)) %>%
  {data.table::fwrite(., "~/AF2_analysis/brinp_v5.csv")}

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

pq_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db"

res_db <- arrow::open_dataset(source = pq_path)

residue_data <- res_db %>%
  filter(uni_gene %in% all_ligands) %>%
  select(uni_gene, uni_sequence) %>%
  collect()





