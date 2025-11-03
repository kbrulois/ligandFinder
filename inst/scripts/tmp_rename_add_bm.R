



rename_data <- tibble(new_dir_name = list.files(run_dir))

run_name = "add_bm"
site = "SU"
submitter = "KB"
algorithm = "AF2v3"
random_seed = 42


rename_data <- parse_afpd_files(input = rename_data,
                                dir_name = "new_dir_name",
                                run_dir = run_dir)

rename_data <- make_new_file_names(input = rename_data,
                                   dir_name = "new_dir_name",
                                   run_name = run_name,
                                   site = site,
                                   submitter = submitter,
                                   algorithm = algorithm,
                                   random_seed = random_seed)

rename_data <- rename_files(run_dir = run_dir,
                            input = rename_data[63:nrow(rename_data), ],
                            from = "og_file_name",
                            to = "new_file_name")
