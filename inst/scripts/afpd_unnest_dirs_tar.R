


unnest_dirs_strict <- function(root = ".") {
  parents <- fs::dir_ls(root, type = "directory")

  for (p in parents) {
    subs <- fs::dir_ls(p, type = "directory", recurse = FALSE)

    if (length(subs) != 1) next
    if (fs::path_file(subs) != fs::path_file(p)) next

    files <- fs::dir_ls(subs, recurse = FALSE)

    # refuse to overwrite anything
    if (any(fs::file_exists(fs::path(p, fs::path_file(files))))) {
      warning("Name collision in ", p, " — skipping")
      next
    }

    fs::file_move(files, p)
    fs::dir_delete(subs)
  }
}

purrr::walk(run_dirs[-11], unnest_dirs_strict)





  to_rename <- jobs[["unknown"]]

  run_dir_rename <- jobs[["unknown"]][["run_dir"]][[1]]

  rename_data <- parse_dirname(run_dir = run_dir_rename,
                               afpd_raw = FALSE,
                               delim_proteins = "_",
                               delim_ranges = "x",
                               delim_start_end = "x")

  rename_data <- make_new_dirname(input = rename_data,
                                  delim_proteins = "_",
                                  delim_ranges = "x",
                                  delim_start_end = "x",
                                  p1_prefix = "h",
                                  p1_suffix = NA,
                                  p2_prefix = "h",
                                  exclude_p1_range = TRUE)

  rename_dir(run_dir = run_dir_rename,
             input = rename_data,
             from = "afpd_dir_name",
             to = "new_dir_name")


  rename_data <- rename_data %>%
                  mutate(file_names = map(fs::path(run_dir_rename, afpd_dir_name), list.files))

  rename_data <- parse_afpd_files(input = rename_data,
                                  dir_name = "new_dir_name",
                                  run_dir = run_dir_rename)

  rename_data <- make_new_file_names(input = rename_data,
                                     dir_name = "new_dir_name",
                                     run_name = run_id,
                                     site = "SU",
                                     submitter = "KB",
                                     algorithm = "AF2v3",
                                     random_seed = 42)

  rename_data <- rename_files(run_dir = run_dir_rename,
                              input = rename_data,
                              from = "og_file_name",
                              to = "new_file_name")







}












