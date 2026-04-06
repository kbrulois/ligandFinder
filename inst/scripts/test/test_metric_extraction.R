
# test_metric_extraction.R
# Run metric extraction on a single directory locally (no sbatch / SLURM needed)
# Usage: set `test_dir` below and source this script in R / RStudio

# --- 1. Load package --------------------------------------------------------
devtools::load_all(".")

library(dplyr)
library(purrr)

# --- 2. User settings -------------------------------------------------------

# Path to a single renamed AFPD directory containing:
#   - ranking_debug.json
#   - 5 *_ark_*.pdb (or *_xtr_*.pdb) files
#   - 5 *_pae_*.json files
test_dir <- ""   # <-- SET THIS (e.g. "/path/to/models/bm/hAGTR2_hANGTx25x31")

# Voronota binary (needed for contact extraction)
voronota_path <- "/home/groups/ebutcher/programs/voronota/bin/voronota-contacts"

# Algorithm label
alg <- "AF2v3"

# --- 3. Validate inputs -----------------------------------------------------

stopifnot(
  "test_dir must be set" = nzchar(test_dir),
  "test_dir does not exist" = dir.exists(test_dir),
  "ranking_debug.json not found" = file.exists(file.path(test_dir, "ranking_debug.json")),
  "voronota binary not found" = file.exists(voronota_path)
)

# --- 4. Load reference data --------------------------------------------------

gpcr_list <- readRDS(system.file("extdata/gpcr_list.rds", package = "ligandFinder"))

# --- 5. Parse protein names from directory ------------------------------------

dir_name <- basename(test_dir)
run_name <- basename(dirname(test_dir))

parsed <- parse_dirname(pairing_dir = dir_name,
                        delim_proteins = "_",
                        delim_ranges = "x",
                        delim_start_end = "x")

proteins <- parsed %>%
  unnest(parsed_pair) %>%
  filter(annotation == "name") %>%
  pull(value)

message("Proteins: ", paste(proteins, collapse = ", "))

# --- 6. Load residue data ----------------------------------------------------

pq_path <- paste0(get_db_path(), "/residue_db")

if (!dir.exists(pq_path)) {
  message("Residue DB not found at ", pq_path)
  message("Run download_residue_db() or set_db_path() first.")
  stop("residue_db not available")
}

res_db <- arrow::open_dataset(source = pq_path)

residue_data <- res_db %>%
  filter(uni_gene %in% proteins) %>%
  collect()

message("Loaded residue data: ", nrow(residue_data), " rows for ",
        length(unique(residue_data$uni_gene)), " protein(s)")

# --- 7. Run metric extraction -------------------------------------------------

message("\n--- Running do_metrics() on: ", dir_name, " ---\n")

result <- do_metrics(directory = test_dir,
                     job = "test",
                     res_dat = residue_data,
                     run_name = run_name)

# --- 8. Check outputs ---------------------------------------------------------

metrics_file <- file.path(test_dir, "metrics_v2.csv")
contacts_file <- file.path(test_dir, "metrics_v2c.rds")

if (is.character(result)) {
  message("\ndo_metrics() returned an error:\n  ", result)
} else {
  # Check metrics_v2.csv
  if (file.exists(metrics_file)) {
    metrics <- data.table::fread(metrics_file)
    message("\nmetrics_v2.csv: ", nrow(metrics), " rows x ", ncol(metrics), " cols")
    print(metrics %>% select(any_of(c("code", "run_name", "model_e", "iptm",
                                       "iptm+ptm", "seq_match", "seq_match2",
                                       "depth", "radius", "location"))))
  } else {
    message("\nmetrics_v2.csv was NOT written")
  }

  # Check metrics_v2c.rds
  if (file.exists(contacts_file)) {
    contacts <- readRDS(contacts_file)
    n_non_null <- sum(!sapply(contacts, is.null))
    message("\nmetrics_v2c.rds: ", length(contacts), " entries, ",
            n_non_null, " with contact data")
    if (n_non_null > 0) {
      first_contact <- contacts[!sapply(contacts, is.null)][[1]]
      message("First contact table: ", nrow(first_contact), " rows x ",
              ncol(first_contact), " cols")
      print(head(first_contact))
    }
  } else {
    message("\nmetrics_v2c.rds was NOT written")
  }
}
