#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Download AlphaFold-DB structures for every GPCR in ligandFinder's list.
##
## For each unique uniprot_id in inst/extdata/gpcr_list.rds, pulls:
##   - <cache>/pdb/AF-<UNIPROT>-F1-model_v4.pdb           (pLDDT in B-factor)
##   - <cache>/pae/AF-<UNIPROT>-F1-predicted_aligned_error_v4.json
## from the EBI AlphaFold Protein Structure Database (https://alphafold.ebi.ac.uk).
##
## Writes <cache>/index.tsv with one row per uniprot_id (status + paths).
## Idempotent: skips files already present unless --refresh is passed.
##
## Typical Sherlock invocation:
##   Rscript inst/scripts/gpcr_af_cache_download.R \
##     --out /oak/stanford/groups/ebutcher/deorphan-AI-ze/gpcr_af_db
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-o", "--out"), type = "character",
              default = "/oak/stanford/groups/ebutcher/deorphan-AI-ze/gpcr_af_db",
              help = "Output cache directory [default: %default]"),
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "Path to gpcr_list.rds. Defaults to ligandFinder's built-in copy."),
  make_option("--refresh", action = "store_true", default = FALSE,
              help = "Re-download even if files exist locally"),
  make_option("--afver", type = "character", default = "v4",
              help = "EBI AlphaFold DB version tag [default: %default]"),
  make_option("--sleep", type = "double", default = 0.1,
              help = "Seconds to sleep between requests, to be polite to EBI [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

list_path <- opt$list
if (is.null(list_path) || !nzchar(list_path)) {
  list_path <- system.file("extdata", "gpcr_list.rds", package = "ligandFinder")
  if (!nzchar(list_path)) {
    list_path <- "inst/extdata/gpcr_list.rds"
  }
}
if (!file.exists(list_path)) {
  stop("gpcr_list.rds not found. Pass --list /path/to/gpcr_list.rds.")
}

pdb_dir <- file.path(opt$out, "pdb")
pae_dir <- file.path(opt$out, "pae")
dir.create(pdb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pae_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://alphafold.ebi.ac.uk/files"
pdb_name <- function(u) sprintf("AF-%s-F1-model_%s.pdb", u, opt$afver)
pae_name <- function(u) sprintf("AF-%s-F1-predicted_aligned_error_%s.json", u, opt$afver)

gpcr <- readRDS(list_path)
if (!"uniprot_id" %in% names(gpcr)) {
  stop("Expected a 'uniprot_id' column in ", list_path)
}
ids <- unique(stats::na.omit(gpcr$uniprot_id))
message(sprintf("GPCR list: %d rows, %d unique uniprot_ids, %d NA (skipped).",
                nrow(gpcr), length(ids), sum(is.na(gpcr$uniprot_id))))
message("Cache dir: ", opt$out)
message("Mode: ", if (opt$refresh) "refresh (re-download all)" else "incremental (skip existing)")

download_one <- function(url, dest) {
  tryCatch({
    suppressWarnings(
      utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
    )
    file.exists(dest) && file.info(dest)$size > 0
  }, error = function(e) {
    if (file.exists(dest)) file.remove(dest)
    FALSE
  })
}

results <- vector("list", length(ids))
t0 <- Sys.time()
for (i in seq_along(ids)) {
  u <- ids[[i]]
  pdb_path <- file.path(pdb_dir, pdb_name(u))
  pae_path <- file.path(pae_dir, pae_name(u))

  pdb_ok <- if (!opt$refresh && file.exists(pdb_path) && file.info(pdb_path)$size > 0) TRUE
            else download_one(file.path(base_url, pdb_name(u)), pdb_path)
  pae_ok <- if (!opt$refresh && file.exists(pae_path) && file.info(pae_path)$size > 0) TRUE
            else download_one(file.path(base_url, pae_name(u)), pae_path)

  results[[i]] <- data.frame(
    uniprot_id    = u,
    pdb_path      = if (pdb_ok) pdb_path else NA_character_,
    pae_path      = if (pae_ok) pae_path else NA_character_,
    pdb_status    = if (pdb_ok) "ok" else "missing",
    pae_status    = if (pae_ok) "ok" else "missing",
    downloaded_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE
  )

  if (i %% 25 == 0 || i == length(ids)) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf("  %4d / %4d  (%.1fs)  last: %s pdb=%s pae=%s",
                    i, length(ids), elapsed, u,
                    if (pdb_ok) "ok" else "MISS",
                    if (pae_ok) "ok" else "MISS"))
  }
  if (opt$sleep > 0) Sys.sleep(opt$sleep)
}

idx <- do.call(rbind, results)
n_pdb_ok <- sum(idx$pdb_status == "ok")
n_pae_ok <- sum(idx$pae_status == "ok")
message(sprintf("\nDone. PDB: %d / %d  PAE: %d / %d", n_pdb_ok, nrow(idx), n_pae_ok, nrow(idx)))

miss <- idx$uniprot_id[idx$pdb_status != "ok"]
if (length(miss)) {
  message("Missing in EBI (", length(miss), "): ",
          paste(utils::head(miss, 20), collapse = ", "),
          if (length(miss) > 20) ", ..." else "")
}

idx_path <- file.path(opt$out, "index.tsv")
utils::write.table(idx, idx_path, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "")
message("Wrote index: ", idx_path)
