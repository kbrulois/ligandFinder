

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
library(ligandFinder)

args <- commandArgs(trailingOnly = TRUE)
job_dir <- args[1]
task_id <- args[2]

group_file <- file.path(job_dir, paste0("group_", task_id, ".txt"))

if(!file.exists(group_file)) {
  stop("Group file not found: ", group_file)
}

tasks <- data.table::fread(group_file, header = FALSE, sep = "\t",
                           col.names = c("afpd_dir", "p1_name", "p2_name", "run_id"))

set_db_path("/scratch/groups/ebutcher/deorphan/ligandFinder")
pq_path <- "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db"
pq_path_src <- "/home/groups/ebutcher/kevin/ligandFinder/residue_db"
voronota_path <- "/home/groups/ebutcher/programs/voronota/bin/voronota-contacts"

resolve_pq_path <- function(p) {
  if(!dir.exists(p)) return(NA_character_)
  if(dir.exists(file.path(p, "residue_db"))) p <- file.path(p, "residue_db")
  schema_cols <- tryCatch(names(arrow::open_dataset(p)$schema), error = function(e) character(0))
  if("uni_gene" %in% schema_cols) p else NA_character_
}

pq_path_use <- resolve_pq_path(pq_path)
if(is.na(pq_path_use)) {
  message("residue_db at ", pq_path, " missing/empty/wrong-schema; falling back to ", pq_path_src)
  pq_path_use <- resolve_pq_path(pq_path_src)
}
if(is.na(pq_path_use)) {
  stop("No usable residue_db found at ", pq_path, " or ", pq_path_src)
}
message("Using residue_db at ", pq_path_use)

alg <- "AF2v3"

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))
bw_align <- summarize_bw(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

proteins <- unique(c(tasks$p1_name, tasks$p2_name))

res_db <- arrow::open_dataset(source = pq_path_use)
residue_data <- res_db %>%
  filter(uni_gene %in% proteins) %>%
  collect()

log_file <- file.path(job_dir, paste0("log_", task_id, ".tsv"))

message("Processing ", nrow(tasks), " directories (task ", task_id, ")")

results <- character(nrow(tasks))

for(i in seq_len(nrow(tasks))) {
  message("[", i, "/", nrow(tasks), "] ", tasks$afpd_dir[i])
  status <- tryCatch({
    do_metrics(directory = tasks$afpd_dir[i],
               job = paste0("task_", task_id),
               res_dat = residue_data,
               run_name = tasks$run_id[i])
    "success"
  }, error = function(e) paste0("error: ", conditionMessage(e)))
  results[i] <- status
}

writeLines(
  paste(tasks$afpd_dir, results, sep = "\t"),
  log_file
)

n_fail <- sum(results != "success")
message("Task ", task_id, " complete. ", nrow(tasks) - n_fail, "/", nrow(tasks), " succeeded.")
if(n_fail > 0) message("Failures logged to ", log_file)
