
# test_residue_db.R
# Check that the residue_db parquet dataset is reachable and has the expected schema.
# Source on the cluster (login or compute node) before submitting metric jobs.

.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')

pq_paths <- c(
  scratch = "/scratch/groups/ebutcher/deorphan/ligandFinder/residue_db",
  home    = "/home/groups/ebutcher/kevin/ligandFinder/residue_db"
)

required_cols <- c("uni_gene", "cons", "af_missense")
probe_gene    <- "AGTR2"   # any known GPCR; just needs to return > 0 rows

check_one <- function(label, path) {
  cat("\n[", label, "] ", path, "\n", sep = "")

  if(!dir.exists(path)) {
    cat("  FAIL: directory does not exist\n"); return(invisible(FALSE))
  }
  if(dir.exists(file.path(path, "residue_db"))) {
    cat("  WARN: nested residue_db/residue_db/ found - dir_copy(overwrite=TRUE) foot-gun\n")
  }

  pq_files <- list.files(path, pattern = "\\.parquet$", recursive = TRUE)
  cat("  parquet files: ", length(pq_files), "\n", sep = "")
  if(length(pq_files) == 0) {
    cat("  FAIL: no parquet files under this path\n"); return(invisible(FALSE))
  }

  ds <- tryCatch(arrow::open_dataset(path),
                 error = function(e) { cat("  FAIL: open_dataset errored: ", conditionMessage(e), "\n", sep = ""); NULL })
  if(is.null(ds)) return(invisible(FALSE))

  cols <- names(ds$schema)
  missing <- setdiff(required_cols, cols)
  if(length(missing) > 0) {
    cat("  FAIL: missing columns: ", paste(missing, collapse = ", "), "\n", sep = "")
    cat("  schema has: ", paste(cols, collapse = ", "), "\n", sep = "")
    return(invisible(FALSE))
  }

  n <- tryCatch(
    ds |> dplyr::filter(uni_gene == probe_gene) |> dplyr::count() |> dplyr::collect() |> dplyr::pull(n),
    error = function(e) { cat("  FAIL: probe query errored: ", conditionMessage(e), "\n", sep = ""); NA_integer_ }
  )
  if(is.na(n)) return(invisible(FALSE))
  if(n == 0) {
    cat("  WARN: probe gene '", probe_gene, "' returned 0 rows (column present but data may be incomplete)\n", sep = "")
  } else {
    cat("  OK: probe gene '", probe_gene, "' returned ", n, " rows\n", sep = "")
  }
  invisible(TRUE)
}

results <- Map(check_one, names(pq_paths), pq_paths)

cat("\nsummary: ", sum(unlist(results)), " / ", length(results), " paths usable\n", sep = "")
