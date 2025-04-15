#!/usr/bin/env Rscript

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(ligandFinder)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help="path to input directory (e.g. generated using afpd_create_job_input_dir.R) containing job1.txt, job2.txt, ...", metavar="INPUT"),
  make_option(c("-o", "--output"), type = "character", help="path to output directory to store alpha fold models", metavar="OUTPUT"),
  make_option(c("-s", "--script"), type="character", help="path to alphapulldown script", default = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/scripts/alphapulldown.sh", metavar="SCRIPT"),
  make_option(c("-f", "--first"), type="integer", help="first job number", metavar="FIRST"),
  make_option(c("-l", "--last"), type="integer", help="last job number", metavar="LAST"),
  make_option(c("-m", "--max_jobs"), type = "integer", help="maximum number of jobs running or in queue", default=16, metavar="MAXJOBS"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print extra output")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (opt$verbose) {
  message("Running script with:")
  print(opt)
}

submit_model_jobs(script = opt[["script"]],
                  protein_input_path = opt[["input"]],
                  output_path = opt[["output"]],
                  jobs = paste0("job", c(opt[["first"]]:opt[["last"]]), ".txt"),
                  max_jobs = opt[["max_jobs"]])




