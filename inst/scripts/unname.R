



rename_data <- data.table::fread(paste0(input_path_models, "/", "hMC4R_hCARTx79x116/file_name_log.csv"))


rename_data <- rename_files(run_dir = input_path_models,
                            input = tibble(new_dir_name = "hMC4R_hCARTx79x116", files = list(as_tibble(rename_data))),
                            from = "new_file_name",
                            to = "og_file_name")
