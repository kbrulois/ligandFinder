#!/usr/bin/env Rscript

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(ligandFinder)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help="path to input directory", metavar="INPUT"),
  make_option(c("-n", "--name"), type="character", help="name of run", metavar="NAME"),
  make_option(c("-p", "--person"), type="character", help="person submitting job", metavar="PERSON"),
  make_option(c("-l", "--location"), type="character", help="location of run", metavar="LOCATION"),
  make_option(c("-a", "--algorithm"), type="character", help="algorithm", default="AF3", metavar="ALGORITHM"),
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
                             delim_proteins = "_",
                             delim_ranges = "x",
                             delim_start_end = "x",
                             num_proteins = 2)


rename_data <- parse_af3_files(input = rename_data,
                               dir_name = "afpd_dir_name",
                               run_dir = run_directory)

rename_data <- make_new_af3_file_names(input = rename_data,
                                       dir_name = "afpd_dir_name",
                                       run_name = opt[["name"]],
                                       site = opt[["location"]],
                                       submitter = opt[["person"]],
                                       algorithm = opt[["algorithm"]])

rename_data <- rename_files(run_dir = run_directory,
                            input = rename_data,
                            dir_name = "afpd_dir_name",
                            from = "og_file_name",
                            to = "new_file_name")


