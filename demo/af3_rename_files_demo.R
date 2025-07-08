


#rename_test_dir <- download_rename_demo()

rename_test_dir <- "~/peptide_alg/AF3_test"


rename_data <- parse_dirname(run_dir = rename_test_dir,
                             delim_proteins = "_",
                             delim_ranges = "x",
                             delim_start_end = "x")


rename_data <- parse_af3_files(input = rename_data,
                                dir_name = "afpd_dir_name",
                                run_dir = rename_test_dir)

rename_data <- make_new_af3_file_names(input = rename_data,
                                   dir_name = "afpd_dir_name",
                                   run_name = "testRename",
                                   site = "SD",
                                   submitter = "SE",
                                   algorithm = "AF3")

rename_data <- rename_files(run_dir = rename_test_dir,
                            input = rename_data,
                            dir_name = "afpd_dir_name",
                            from = "og_file_name",
                            to = "new_file_name")







