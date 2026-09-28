## Attaching HPA single-cell (enteroendocrine) expression to a window table.
##
## The join is on GENE SYMBOL, and it has to be: the HPA classification carries
## only `gene` (its own symbol) and `ensembl`, with no UniProt accession and no
## entry name, while `data/id_mapping.rds` carries Entry / Entry Name / Gene
## Names but no Ensembl. There is therefore no accession-level bridge between the
## two locally, and the symbol is the only shared identifier.
##
## Measured on the peptide-end window set (5,573 genes): direct symbol match
## covers 97.1%, the guarded synonym rescue adds 6 genes, and routing through
## UniProt's primary gene name instead is slightly WORSE (96.8%) while rescuing
## nothing -- so the plain symbol, plus the rescue below, is the best available
## key.

#' Reconcile window-table gene symbols against HPA gene symbols
#'
#' Direct symbol matches, plus a rescue through HPA's own gene synonyms -- but
#' only where the synonym is unambiguous (one HPA Ensembl gene) AND the HPA
#' target is not already present in the window table under its own name.
#'
#' Without that second guard, "HLA-H" (the HLA-region pseudogene) is silently
#' handed HFE's expression profile, because HFE was once called HLA-H. The guard
#' is the whole reason this is a function rather than a `left_join`.
#'
#' @param genes the window table's gene symbols (unique or not).
#' @param hpa_genes the HPA classification's `gene` column.
#' @param proteinatlas path to HPA's `proteinatlas.tsv` (for the synonyms); the
#'   rescue is skipped when it is NA or missing.
#' @return data.table `pred_gene`, `hpa_gene`, `match_type` ("direct"/"synonym").
#' @export
lf_eec_gene_lut <- function(genes, hpa_genes, proteinatlas = NA_character_) {
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required", call. = FALSE)
  dt <- data.table::data.table
  genes  <- unique(as.character(genes))
  direct <- intersect(genes, as.character(hpa_genes))
  lut    <- dt(pred_gene = direct, hpa_gene = direct, match_type = "direct")

  missing <- setdiff(genes, as.character(hpa_genes))
  if (is.na(proteinatlas) || !file.exists(proteinatlas) || !length(missing))
    return(lut)

  pa <- data.table::fread(proteinatlas, quote = "\"",
                          select = c("Gene", "Gene synonym", "Ensembl"))
  data.table::setnames(pa, c("gene", "synonym", "ensembl"))
  syn <- pa[nzchar(synonym),
            .(s = trimws(unlist(strsplit(synonym, ",")))), by = .(gene, ensembl)]
  syn <- syn[nzchar(s) & s %in% missing]
  if (!nrow(syn)) return(lut)

  ## one HPA gene per synonym ...
  unambiguous <- syn[, .(n = data.table::uniqueN(ensembl)), by = s][n == 1, s]
  syn <- unique(syn[s %in% unambiguous], by = "s")
  ## ... and the HPA target must not already stand on its own in the table
  syn <- syn[!gene %in% genes]
  if (!nrow(syn)) return(lut)

  data.table::rbindlist(list(
    lut, dt(pred_gene = syn$s, hpa_gene = syn$gene, match_type = "synonym")))
}

#: the EEC columns carried over, and the names they take in the output
LF_EEC_COLS <- c(
  ntpm_eec             = "eec_ntpm",
  class                = "eec_class",
  specificity_tier     = "eec_tier",
  rank_among_celltypes = "eec_rank_among_celltypes",
  n_celltypes_detected = "eec_n_celltypes_detected",
  top_other_celltype   = "eec_top_other_celltype",
  fc_vs_max_other      = "eec_fc_vs_max_other",
  fc_vs_mean_other     = "eec_fc_vs_mean_other",
  tau                  = "eec_tau",
  hpa_sc_specificity   = "eec_hpa_specificity",
  is_secreted          = "is_secreted",
  is_gpcr              = "is_gpcr")

#' Attach the EEC expression block to a window table
#'
#' A gene with no HPA record gets `eec_class` `"no_expression_data"`, kept
#' distinct from `"not_expressed"`: no measurement is not a measured zero.
#'
#' @param pred window table with a `gene` column; one row per window.
#' @param eec the HPA EEC classification (`eec_gene_classification.tsv`).
#' @param proteinatlas path to `proteinatlas.tsv`, for the synonym rescue.
#' @return `pred` with the EEC columns and `eec_match_type` added, same row count
#'   and order.
#' @export
lf_eec_attach <- function(pred, eec, proteinatlas = NA_character_) {
  stopifnot("gene" %in% names(pred), "gene" %in% names(eec))
  pred <- data.table::as.data.table(pred)
  eec  <- data.table::as.data.table(eec)
  have <- intersect(names(LF_EEC_COLS), names(eec))
  clash <- intersect(unname(LF_EEC_COLS[have]), names(pred))
  if (length(clash))
    stop("these columns are already present: ", paste(clash, collapse = ", "),
         call. = FALSE)

  lut <- lf_eec_gene_lut(pred$gene, eec$gene, proteinatlas)
  slim <- eec[, c("gene", have), with = FALSE]
  data.table::setnames(slim, c("hpa_gene", unname(LF_EEC_COLS[have])))

  ## `sort = FALSE` is not a promise about order, so carry an explicit index and
  ## restore it -- callers join these columns back onto other tables by position
  n0 <- nrow(pred)
  pred <- data.table::copy(pred)[, lf_row_idx__ := .I]
  out <- merge(pred, lut, by.x = "gene", by.y = "pred_gene", all.x = TRUE)
  out <- merge(out, slim, by = "hpa_gene", all.x = TRUE)
  ## a left join must never duplicate or drop a window
  if (nrow(out) != n0)
    stop(sprintf("the join changed the row count (%d -> %d); a gene symbol is ",
                 n0, nrow(out)), "not unique in the EEC table", call. = FALSE)
  data.table::setorderv(out, "lf_row_idx__")
  out[, lf_row_idx__ := NULL]

  out[, eec_match_type := ifelse(is.na(match_type), "none", match_type)]
  out[, match_type := NULL]
  if ("eec_class" %in% names(out))
    out[is.na(eec_class), eec_class := "no_expression_data"]
  if ("eec_tier" %in% names(out))
    out[is.na(eec_tier), eec_tier := "no_expression_data"]
  out[]
}
