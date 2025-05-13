#!/usr/bin/env Rscript

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(ligandFinder)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help="path to input directory", metavar="INPUT"),
  make_option(c("-n", "--name"), type="character", help="name of run", metavar="NAME"),
  make_option(c("-p", "--person"), type="character", help="person submitting job", metavar="PERSON"),
  make_option(c("-l", "--location"), type="character", help="location of run", metavar="LOCATION"),
  make_option(c("-a", "--algorithm"), type="character", help="algorithm", defaul="AF2v3", metavar="ALGORITHM"),
  make_option(c("-s", "--seed"), type = "integer", help="random seed", default=42, metavar="SEED"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print extra output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (opt$verbose) {
  message("Running script with:")
  print(opt)
}

run_directory <- opt[["input"]]

rename_data <- parse_dirname(run_dir = run_directory,
                             afpd_raw = TRUE)

rename_data <- make_new_dirname(input = rename_data,
                                delim_proteins = "_",
                                delim_ranges = "x",
                                delim_start_end = "x",
                                p1_prefix = "h",
                                p1_suffix = NA,
                                p2_prefix = "h",
                                exclude_p1_range = TRUE)

rename_dir(run_dir = run_directory,
           input = rename_data,
           from = "afpd_dir_name",
           to = "new_dir_name")


rename_data <- parse_afpd_files(input = rename_data,
                                dir_name = "new_dir_name",
                                run_dir = run_directory)

rename_data <- make_new_file_names(input = rename_data,
                                   dir_name = "new_dir_name",
                                   run_name = opt[["name"]],
                                   site = opt[["location"]],
                                   submitter = opt[["person"]],
                                   algorithm = opt[["algorithm"]],
                                   random_seed = opt[["seed"]])

rename_data <- rename_files(run_dir = run_directory,
                            input = rename_data,
                            from = "og_file_name",
                            to = "new_file_name")




