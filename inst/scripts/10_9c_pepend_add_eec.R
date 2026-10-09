#!/usr/bin/env Rscript
## ---- attach EEC (HPA single-cell) expression to the peptide-end windows --------
## Adds the enteroendocrine-cell expression block from
## inst/scripts/enteroendocrine_expression.R's classification to every window of
## 10_9b's table, one row per window.
##
## THE JOIN IS ON GENE SYMBOL, and it has to be. Verified on this data:
##   * the HPA classification carries `gene` (its own symbol) and `ensembl` --
##     no UniProt accession and no entry name
##   * data/id_mapping.rds carries Entry / Entry Name / Gene Names -- no Ensembl
##   so there is no accession-level bridge locally and the symbol is the only
##   shared identifier. Entry name is not an option: it is absent from the HPA
##   side entirely.
##
## Symbol matching is also good enough here, measured over the 5,573 genes:
##   direct symbol                       5,414 (97.1%)
##   + HPA synonym rescue, guarded          +6 (97.3%)
##   routing via UniProt's primary gene  5,394 (96.8%) -- WORSE, and rescues 0,
##                                       so it is deliberately not used
## The 151 that stay unmatched are overwhelmingly not secreted-peptide genes
## (mitochondrial rRNA, Ig segments, pseudogenes, LOC records).
##
## Reconciliation lives in R/eec_expression.R, including the guard that stops
## HLA-H being handed HFE's profile via a stale synonym.
##
##   Rscript inst/scripts/10_9c_pepend_add_eec.R --term C
##
## Output: ~/AF2_analysis/lf_pepend_all_windows_<term>_eec.csv
## ------------------------------------------------------------------------------

suppressMessages({ library(data.table); library(dplyr) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term    <- toupper(.opt("--term", "C"))
in_csv  <- path.expand(.opt("--input", sprintf("~/AF2_analysis/lf_pepend_all_windows_%s.csv", term)))
eec_p   <- path.expand(.opt("--eec", "~/AF2_analysis/eec_gene_classification.tsv"))
pa_p    <- path.expand(.opt("--proteinatlas", "~/AF2_analysis/hpa_cache/proteinatlas.tsv"))
out_csv <- path.expand(.opt("--out", sub("\\.csv$", "_eec.csv", in_csv)))
score   <- .opt("--score-col", "pepend_mean_allseeds")

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "eec_expression.R"))

pred <- fread(in_csv)
eec  <- fread(eec_p)
stopifnot(score %in% names(pred))
message(sprintf("%s windows, %d genes  <-  %s",
                format(nrow(pred), big.mark = ","), uniqueN(pred$gene), basename(in_csv)))
message(sprintf("%s HPA gene records  <-  %s",
                format(nrow(eec), big.mark = ","), basename(eec_p)))

## ---- the key, stated and checked ---------------------------------------------
message("\nidentifiers available for the join:")
message("  window table : ", paste(intersect(names(pred), c("gene", "accession")), collapse = ", "))
message("  HPA table    : ", paste(intersect(names(eec), c("gene", "ensembl")), collapse = ", "))
message("  -> joining on GENE SYMBOL (no accession or entry name on the HPA side)")

lut <- lf_eec_gene_lut(pred$gene, eec$gene, pa_p)
ng  <- uniqueN(pred$gene)
message(sprintf("\nsymbol match: direct %d | synonym %d | unmatched %d  of %d genes",
                sum(lut$match_type == "direct"), sum(lut$match_type == "synonym"),
                ng - nrow(lut), ng))
if (any(lut$match_type == "synonym"))
  print(as.data.frame(lut[match_type == "synonym"]), row.names = FALSE)

out <- lf_eec_attach(pred, eec, pa_p)
stopifnot(nrow(out) == nrow(pred), identical(out$peps, pred$peps))

## keep the window identity first, then scores, then the expression block
front <- intersect(c("peps", "gene", "accession", "hpa_gene", "eec_match_type"), names(out))
out <- out[, c(front, setdiff(names(out), front)), with = FALSE]
fwrite(out, out_csv)
message(sprintf("\n%s rows x %d columns -> %s", format(nrow(out), big.mark = ","),
                ncol(out), out_csv))

## ---- what the merge produced --------------------------------------------------
message("\nwindows by EEC class:")
out[, .(windows = .N, genes = uniqueN(gene),
        median_score = round(median(get(score), na.rm = TRUE), 4),
        max_score = round(max(get(score), na.rm = TRUE), 4),
        known = sum(known == 1)), by = eec_class][order(-windows)] %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\nEEC-specific genes: highest-scoring windows")
out[eec_class == "specifically_expressed"][order(-get(score))][1:20,
    .(gene, anchor, win_type, score = round(get(score), 4),
      n_gt = get(grep("^pepend_n_seeds_gt", names(out), value = TRUE)[1]),
      eec_ntpm = round(eec_ntpm, 1), eec_tier, is_gpcr,
      known = ifelse(is.na(pep_name), "", pep_name))] %>%
  as.data.frame() %>% print(row.names = FALSE)

message("\nNOVEL (known == 0) windows in EEC-specific genes, top 20")
out[eec_class == "specifically_expressed" & known == 0][order(-get(score))][1:20,
    .(gene, anchor, win_type, score = round(get(score), 4),
      eec_ntpm = round(eec_ntpm, 1), eec_fc_vs_max_other = round(eec_fc_vs_max_other, 1),
      eec_tier, scan_rank_peak)] %>%
  as.data.frame() %>% print(row.names = FALSE)
