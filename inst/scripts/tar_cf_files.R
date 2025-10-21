



.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder", auth_token = "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")
library(ligandFinder)


num_of_grps <- 16

future::plan(strategy = future::multicore(workers = num_of_grps))


run_dir <- "/scratch/groups/ebutcher/deorphan/models/top200NC"


tmp <- tibble(fs::dir_info(run_dir))

tmp <- tmp %>%
mutate(ranking_debug_exists = furrr::future_map_lgl(path, \(x) {
  file.exists(paste(x, "ranking_debug.json", sep = "/"))
})) %>%
  mutate(data_files = furrr::future_map(path, ~fs::dir_ls(.) %>% basename)) %>%
  mutate(num_files = map_int(data_files, length))


tmp <- tmp %>%
  mutate(complete = num_files == 32 & ranking_debug_exists)


cutoff <- Sys.time() - 60 * 60


tmp <- tmp %>% filter(modification_time < cutoff & complete)

tar_dir <- paste0(run_dir, "_tar")

tmp_tar <- tibble(fs::dir_info(tar_dir))

already_tarred <- tmp_tar[["path"]] %>% basename %>% stringr::str_remove(., ".tar$")

tmp <- tmp %>%
  mutate(comp_name = path %>% basename) %>%
  filter(!comp_name %in% already_tarred)

#dir.create(tar_dir)

tar_one <- function(dir) {
  tarfile <- fs::path(tar_dir, paste0(fs::path_file(dir), ".tar"))
  cmd <- sprintf("tar -cf %s -C %s %s", tarfile, fs::path_dir(dir), fs::path_file(dir))
  system(cmd)
  tarfile
}

result_files <- furrr::future_map(tmp[["path"]], tar_one, .progress = TRUE)

out_dir <- "top200NCnew"

dir.create("top200NCnew")

untar_one <- function(tarfile) {
  cmd <- sprintf("tar -xf %s -C %s", tarfile, out_dir)
  system(cmd)
}
result_files <- furrr::future_map(tmp_tar[["path"]], untar_one, .progress = TRUE)


