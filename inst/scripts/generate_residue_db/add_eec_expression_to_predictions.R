#!/usr/bin/env Rscript
#
# Append the full enteroendocrine expression block to a predictions CSV in place.
# Rows are never re-serialised: fields are appended textually to the raw lines.
#
# Usage: Rscript add_eec_expression_to_predictions.R <predictions.csv> \
#          <eec_gene_classification.tsv> [proteinatlas.tsv]

suppressPackageStartupMessages(library(data.table))

args      <- commandArgs(trailingOnly = TRUE)
pred_path <- path.expand(args[1])
eec_path  <- path.expand(args[2])
pa_path   <- if (length(args) >= 3) path.expand(args[3]) else NA_character_

lines <- readLines(pred_path, warn = FALSE)
pred  <- data.table::fread(pred_path)
eec   <- data.table::fread(eec_path)
if (length(lines) != nrow(pred) + 1L) stop("line/row mismatch")

add <- c(ntpm_eec = "eec_ntpm", specificity_tier = "eec_tier",
         rank_among_celltypes = "eec_rank_among_celltypes",
         n_celltypes_detected = "eec_n_celltypes_detected",
         top_other_celltype = "eec_top_other_celltype",
         fc_vs_max_other = "eec_fc_vs_max_other",
         fc_vs_mean_other = "eec_fc_vs_mean_other",
         tau = "eec_tau", hpa_sc_specificity = "eec_hpa_specificity",
         is_secreted = "is_secreted", is_gpcr = "is_gpcr")
dup <- intersect(unname(add), names(pred))
if (length(dup)) stop("already present: ", paste(dup, collapse = ", "))

# gene -> row index, with guarded synonym rescue
idx <- match(pred$gene, eec$gene)
if (!is.na(pa_path) && file.exists(pa_path) && anyNA(idx)) {
  miss <- setdiff(unique(pred$gene[is.na(idx)]), eec$gene)
  pa <- data.table::fread(pa_path, quote = "\"",
                          select = c("Gene", "Gene synonym", "Ensembl"))
  setnames(pa, c("gene", "synonym", "ensembl"))
  syn <- pa[nzchar(synonym), .(s = trimws(unlist(strsplit(synonym, ",")))),
            by = .(gene, ensembl)][nzchar(s) & s %in% miss]
  syn <- unique(syn[s %in% syn[, .(n = uniqueN(ensembl)), by = s][n == 1, s]], by = "s")
  syn <- syn[!gene %in% pred$gene]
  if (nrow(syn)) {
    hit <- is.na(idx) & pred$gene %in% syn$s
    idx[hit] <- match(syn$gene[match(pred$gene[hit], syn$s)], eec$gene)
  }
}

fmt <- function(x) {
  if (is.character(x)) ifelse(is.na(x), "NA", paste0('"', x, '"'))
  else ifelse(is.na(x), "NA", as.character(x))
}
block <- do.call(paste, c(lapply(names(add), function(col) fmt(eec[[col]][idx])), sep = ","))

backup <- paste0(pred_path, ".bak-", format(Sys.time(), "%Y%m%d%H%M%S"))
file.copy(pred_path, backup, overwrite = FALSE)
writeLines(c(paste0(lines[1], ",", paste0('"', unname(add), '"', collapse = ",")),
             paste0(lines[-1], ",", block)), pred_path)

check <- data.table::fread(pred_path)
stopifnot(nrow(check) == nrow(pred),
          identical(names(check), c(names(pred), unname(add))),
          identical(check$score_mean, pred$score_mean),
          identical(check$peps, pred$peps))
cat("backup  ", backup, "\nmodified", pred_path, "\nrows", nrow(check),
    "| columns", ncol(pred), "->", ncol(check), "\n")
cat("matched genes:", sum(!is.na(idx)), "/", nrow(pred), "rows\n")
