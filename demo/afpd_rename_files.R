


rename_test_dir <- download_rename_demo()


rename_data <- parse_dirname(run_dir = rename_test_dir,
                             delim_proteins = "_and_",
                             delim_ranges = "_",
                             delim_start_end = "-")

rename_data <- make_new_dirname(input = rename_data,
                                delim_proteins = "_",
                                delim_ranges = "x",
                                delim_start_end = "x",
                                p1_prefix = "h",
                                p1_suffix = NA,
                                p2_prefix = "h",
                                exclude_p1_range = TRUE)

rename_dir(run_dir = rename_test_dir,
           input = rename_data,
           from = "afpd_dir_name",
           to = "new_dir_name")


rename_data <- parse_afpd_files(input = rename_data,
                                dir_name = "new_dir_name",
                                run_dir = rename_test_dir)

rename_data <- make_new_file_names(input = rename_data,
                                   dir_name = "new_dir_name",
                                   run_name = "testRename",
                                   site = "SU",
                                   submitter = "KB",
                                   algorithm = "AF2multV3",
                                   random_seed = 42)

rename_data <- rename_files(run_dir = rename_test_dir,
                            input = rename_data,
                            from = "og_file_name",
                            to = "new_file_name")







