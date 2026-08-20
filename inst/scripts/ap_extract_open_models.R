## =============================================================================
## Extract the selected open-pocket model PDBs from the per-pair tar archives.
##
## Layout on Oak is ONE archive per receptor-peptide pair:
##   <models_root>/<run_name>/<afpd_dir_name>.tar
## e.g. /oak/.../models/ends/h5HT1B_hAEBP1x26x66.tar
##
## The manifest already carries run_name and afpd_dir_name, so every archive
## path is computable. That matters: there are ~300k archives, so anything that
## walks or lists them all would dominate the runtime. We touch only the ~300
## archives that actually hold a selected model.
##
## Run AFTER ap_initial_guess_active_state_receptors.R.
## =============================================================================

## tidyverse the meta-package is awkward to install on Sherlock, and nothing
## here needs it: the verbs come from dplyr, and readr/tools/utils calls are
## namespaced. Added to the front of the search path rather than replacing it,
## and only when present, so this still runs off a laptop checkout.
.lf_lib <- "/home/groups/ebutcher/programs/pipeline/R_libs4.1"
if (dir.exists(.lf_lib)) .libPaths(c(.lf_lib, .libPaths()))

library(dplyr)
library(tibble)

## ---- options ---------------------------------------------------------------
models_root  <- "/oak/stanford/groups/ebutcher/deorphan-AI-ze/models"
manifest_csv <- path.expand("/oak/stanford/groups/ebutcher/deorphan-AI-ze/open_models/open_receptors_manifest.csv")
out_dir      <- path.expand("/oak/stanford/groups/ebutcher/deorphan-AI-ze/open_models")
also_pae     <- FALSE   # TRUE = also pull the matching _pae_*.json alongside each pdb
search_lost  <- TRUE    # if an archive isn't at the expected path, do ONE directory
                        # walk to find it by name (run_name can differ between the
                        # scratch paths in the metrics and the layout on Oak)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
man <- readr::read_csv(manifest_csv, show_col_types = FALSE)

## ---- 1. one row per wanted file, with the archive it should live in --------
want <- man %>% transmute(receptor, conformer, run_name, afpd_dir_name, file = pdb_file)
if (also_pae) {
  want <- bind_rows(want,
    man %>% transmute(receptor, conformer, run_name, afpd_dir_name, file = pae_file))
}
want <- want %>%
  filter(!is.na(file)) %>%
  distinct() %>%
  mutate(tar_path = file.path(models_root, run_name, paste0(afpd_dir_name, ".tar")))

want$tar_exists <- file.exists(want$tar_path)
message(sprintf("%d file(s) wanted from %d archive(s); %d archive(s) at the expected path",
                nrow(want), n_distinct(want$tar_path),
                n_distinct(want$tar_path[want$tar_exists])))

## ---- 2. rescue archives that aren't where we expected them -----------------
## One walk for ALL of them, not one search per miss.
if (search_lost && any(!want$tar_exists)) {
  message("searching ", models_root, " for the missing archives (one walk) ...")
  all_tars <- list.files(models_root, pattern = "\\.tar$", recursive = TRUE, full.names = TRUE)
  by_name  <- setNames(all_tars, tools::file_path_sans_ext(basename(all_tars)))
  lost     <- !want$tar_exists
  found    <- unname(by_name[want$afpd_dir_name[lost]])   # NA where still not found
  want$tar_path[lost] <- ifelse(is.na(found), want$tar_path[lost], found)
  want$tar_exists     <- file.exists(want$tar_path)
  message(sprintf("  recovered %d archive(s)", sum(!is.na(found))))
}

## ---- 3. extract, one archive at a time -------------------------------------
## Each archive is listed once to resolve the member's path inside it (we match
## on basename, so the internal nesting doesn't matter), then extracted in a
## single call.
todo <- want %>%
  filter(tar_exists) %>%
  group_by(tar_path) %>%
  summarise(files = list(unique(file)), .groups = "drop")

stage <- file.path(out_dir, ".stage")
got   <- character(0)

for (i in seq_len(nrow(todo))) {
  tp <- todo$tar_path[i]
  fs <- todo$files[[i]]

  members <- tryCatch(utils::untar(tp, list = TRUE),
                      error   = function(e) { warning("cannot list ", tp, call. = FALSE); character(0) },
                      warning = function(w) character(0))
  hit <- members[basename(members) %in% fs]
  if (!length(hit)) next

  unlink(stage, recursive = TRUE)
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  utils::untar(tp, files = hit, exdir = stage)

  src <- list.files(stage, pattern = "\\.(pdb|json)$", recursive = TRUE, full.names = TRUE)
  ok  <- file.copy(src, file.path(out_dir, basename(src)), overwrite = TRUE)
  got <- c(got, basename(src)[ok])

  if (i %% 100 == 0 || i == nrow(todo))
    message(sprintf("  %d / %d archives, %d file(s) extracted", i, nrow(todo), length(got)))
}
unlink(stage, recursive = TRUE)

## ---- 4. log what happened --------------------------------------------------
log_tbl <- want %>%
  mutate(status = case_when(!tar_exists  ~ "archive_not_found",
                            file %in% got ~ "extracted",
                            TRUE          ~ "not_in_archive"))

readr::write_csv(log_tbl, file.path(out_dir, "extraction_log.csv"))

message("\n", paste(capture.output(print(count(log_tbl, status))), collapse = "\n"))
message("\nextracted ", length(unique(got)), " file(s) -> ", out_dir)
message("log: ", file.path(out_dir, "extraction_log.csv"))

failed <- log_tbl %>% filter(status != "extracted")
if (nrow(failed))
  message("!! ", nrow(failed), " file(s) not extracted -- inspect extraction_log.csv; ",
          "first few:\n  ", paste(head(failed$file, 3), collapse = "\n  "))
