



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)



job_dir <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/top200NC"

out_dir <- c("/scratch/groups/ebutcher/deorphan/models/top200NC",
             "/scratch/groups/ebutcher/deorphan/models/top200NCnew",
             "/scratch/groups/ebutcher/deorphan/models/top200NC_tar")

all_comps <- map(out_dir, list.files)

all_comps[[3]] <- stringr::str_remove(all_comps[[3]], ".tar$")

identical(all_comps[[2]], all_comps[[3]])

all_comps[[2]][!all_comps[[2]] %in% all_comps[[1]]]



dat <- tibble(jobs = list.files(job_dir) %>% stringr::str_subset(., ".txt$"))

dat <- dat[gtools::mixedorder(dat$jobs), ]

dat <- dat %>%
          mutate(model = map(jobs, ~data.table::fread(paste0(job_dir, "/", .), sep = NULL, header = FALSE)))

dat <- dat %>%
  unnest(model)

dat2 <- tibble(complex = unique(all_comps[[1]]))

id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

dat2 <- dat2 %>%
  mutate(parse_proteins(complex)) %>%
  mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Entry"]])[.])) %>%
  mutate(model = paste0(p1_id,
                        if_else(p1_range == "", "", paste0(",", p1_range)),
                        ";",
                        p2_id,
                        if_else(p2_range == "", "", paste0(",", p2_range))))

sum(dat[["V1"]] %in% dat2[["model"]])

dat <- dat %>%
  mutate(complete = if_else(V1 %in% dat2[["model"]], "yes", "no"))

dat %>%
  group_by(jobs) %>%
  filter(all(complete == "no")) %>%
  pull(jobs) %>%
  unique






to_still_run <- dat %>%
  filter(!V1 %in% dat2$model) %>%
  distinct(V1, .keep_all = TRUE)


left_join(to_still_run, id_map %>% rename(receptor_id = `Entry Name`), by = "receptor_id") %>%
  select(receptor_id, Length) %>%
  print(n = 1000)



group_size <- 48
to_still_run <- to_still_run %>%
  mutate(group = rep(paste0("job", 3000:(3000 + ceiling(n() / group_size)), ".txt"), each = group_size, length.out = n()))

to_still_run %>%
  group_by(group) %>%
  group_walk(~ write.table(.x[["V1"]],
                           file = paste0(job_dir, "/", .y[["group"]]),
                           row.names = FALSE,
                           col.names = FALSE,
                           quote = FALSE))








