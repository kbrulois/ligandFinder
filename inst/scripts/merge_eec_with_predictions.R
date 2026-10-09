#!/usr/bin/env Rscript
#
# Attach the enteroendocrine-cell expression classification to the ensemble
# window predictions, one row per candidate cleavage window.
#
# Gene symbols are matched directly where possible and rescued through HPA gene
# synonyms otherwise -- but only when the synonym is unambiguous AND the HPA
# target is not already standing on its own in the prediction table. Without
# that second guard "HLA-H" (the HLA-region pseudogene) is silently handed HFE's
# expression profile, because HFE was once called HLA-H.
#
# Genes with no HPA record get eec_class "no_expression_data", kept distinct
# from "not_expressed": no measurement is not the same as a measured zero.
#
# Usage:
#   Rscript inst/scripts/merge_eec_with_predictions.R <predictions.csv> \
#           <eec_gene_classification.tsv> <out_dir> [proteinatlas.tsv]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("need <predictions.csv> <eec_classification.tsv> <out_dir>")

pred_path <- args[1]
eec_path  <- args[2]
out_dir   <- args[3]
pa_path   <- if (length(args) >= 4) args[4] else NA_character_

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pred <- data.table::fread(pred_path)
eec  <- data.table::fread(eec_path)
stopifnot("gene" %in% names(pred), "gene" %in% names(eec))

## ---- symbol reconciliation ------------------------------------------------

pred_genes <- unique(pred$gene)
direct     <- intersect(pred_genes, eec$gene)
missing    <- setdiff(pred_genes, eec$gene)

lut <- data.table(pred_gene = direct, hpa_gene = direct, match_type = "direct")

if (!is.na(pa_path) && file.exists(pa_path) && length(missing)) {
  pa <- data.table::fread(pa_path, quote = "\"",
                          select = c("Gene", "Gene synonym", "Ensembl"))
  setnames(pa, c("gene", "synonym", "ensembl"))

  syn <- pa[nzchar(synonym),
            .(s = trimws(unlist(strsplit(synonym, ",")))), by = .(gene, ensembl)]
  syn <- syn[nzchar(s) & s %in% missing]

  # one HPA gene per synonym, and the target must not already be its own row
  unambiguous <- syn[, .(n = uniqueN(ensembl)), by = s][n == 1, s]
  syn <- unique(syn[s %in% unambiguous], by = "s")
  syn <- syn[!gene %in% pred_genes]

  if (nrow(syn)) {
    lut <- rbind(lut, data.table(pred_gene = syn$s, hpa_gene = syn$gene,
                                 match_type = "synonym"))
  }
}

## ---- merge ----------------------------------------------------------------

keep <- c("gene", "ntpm_eec", "class", "specificity_tier", "rank_among_celltypes",
          "n_celltypes_detected", "top_other_celltype", "fc_vs_max_other",
          "fc_vs_mean_other", "tau", "hpa_sc_specificity", "is_secreted", "is_gpcr")
eec_slim <- eec[, ..keep]
setnames(eec_slim,
         c("gene", "ntpm_eec", "class", "specificity_tier", "rank_among_celltypes",
           "n_celltypes_detected", "top_other_celltype", "fc_vs_max_other",
           "fc_vs_mean_other", "tau", "hpa_sc_specificity", "is_secreted", "is_gpcr"),
         c("hpa_gene", "eec_ntpm", "eec_class", "eec_tier", "eec_rank_among_celltypes",
           "eec_n_celltypes_detected", "eec_top_other_celltype", "eec_fc_vs_max_other",
           "eec_fc_vs_mean_other", "eec_tau", "eec_hpa_specificity",
           "is_secreted", "is_gpcr"))

out <- pred %>%
  left_join(lut, by = c("gene" = "pred_gene")) %>%
  left_join(eec_slim, by = "hpa_gene") %>%
  mutate(
    eec_match_type = ifelse(is.na(match_type), "none", match_type),
    eec_class      = ifelse(is.na(eec_class), "no_expression_data", eec_class),
    eec_tier       = ifelse(is.na(eec_tier),  "no_expression_data", eec_tier)
  ) %>%
  select(-match_type) %>%
  as.data.table()

stopifnot(nrow(out) == nrow(pred))   # left join must not duplicate windows

out_file <- file.path(out_dir, "ensemble_global_predictions_with_eec.csv")
data.table::fwrite(out, out_file)

## ---- summary --------------------------------------------------------------

cat("\nwindows:", nrow(out), " genes:", uniqueN(out$gene), "\n")
cat("symbol match: direct", length(direct),
    "| synonym", sum(lut$match_type == "synonym"),
    "| unmatched", length(setdiff(pred_genes, lut$pred_gene)), "\n\n")

cat("windows by EEC class\n")
print(out[, .(windows = .N, genes = uniqueN(gene),
              median_score = round(median(score_mean), 3),
              top_score = round(max(score_mean), 3),
              known_peptides = sum(known == 1)), by = eec_class][order(-windows)])

cat("\nhighest-scoring windows in EEC-specific genes\n")
print(head(out[eec_class == "specifically_expressed"][order(-score_mean),
      .(peps, gene, term, score_mean = round(score_mean, 3), known,
        eec_ntpm, eec_tier, is_gpcr)], 25))

cat("\nnovel (known == 0) EEC-specific windows, top 15\n")
print(head(out[eec_class == "specifically_expressed" & known == 0][order(-score_mean),
      .(peps, gene, term, score_mean = round(score_mean, 3),
        eec_ntpm, eec_fc_vs_max_other, eec_tier)], 15))

cat("\nwrote ", out_file, "\n", sep = "")
