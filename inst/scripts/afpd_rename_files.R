



tmp <- parse_dirname()


tmp <- make_new_dirname(tmp)

input <- parse_afpd_files(input = tmp,
                         dir_name = "afpd_dir_name",
                         run_dir = "~/peptide_alg/rename_test")

input2 <- make_new_file_names(input = input)

input3 <- rename_files(run_dir = "~/peptide_alg/rename_test",
                       input = input2)





rename_dir(run_dir = "~/peptide_alg/rename_test",
           input = tmp,
           from = "afpd_dir_name",
           to = "new_dir_name")


rename_dir(run_dir = "~/peptide_alg/rename_test",
           input = tmp,
           to = "afpd_dir_name",
           from = "new_dir_name")






tmp2 <- parse_dirname(delim_proteins = "_",
                      delim_ranges = "x",
                      delim_start_end = "x")

tmp3 <- make_new_dirname(input = tmp2)






