#!/usr/bin/env Rscript
#
# Enteroendocrine-cell gene expression classification.
#
# Pulls the Human Protein Atlas single-cell consensus (nTPM per gene per cell
# type) and classifies every gene with respect to enteroendocrine cells as one
# of:
#
#   not_expressed              below the detection floor in EECs
#   non_specifically_expressed detected, but no better than the rest of the body
#   specifically_expressed     detected and enriched/enhanced over other cell types
#
# HPA v23 is used because it is the last release carrying a standalone
# "Enteroendocrine cells" cluster; v24+ merges it into "neuroendocrine cells",
# which pools gut EECs with lung and other neuroendocrine populations.
#
# Usage:  Rscript inst/scripts/enteroendocrine_expression.R [out_dir] [cache_dir]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

## ---- parameters -----------------------------------------------------------

TARGET     <- "Enteroendocrine cells"
DETECT     <- 1     # nTPM >= this counts as detected (HPA convention)
FOLD       <- 4     # fold-change for enrichment (HPA convention)
GROUP_MAX  <- 10    # largest set allowed for "group enriched"
DOM_TOL    <- 2     # another cell type may not beat EECs by more than this fold
                    # for a gene to count as EEC-specific (see the guard below)
BREADTH    <- 1/3   # a gene detected in more than this share of cell types is too
                    # widely expressed to be called specific on the enhanced rule
                    # alone (HPA's "detected in some" / "in many" boundary)

args      <- commandArgs(trailingOnly = TRUE)
out_dir   <- if (length(args) >= 1) args[1] else "."
cache_dir <- if (length(args) >= 2) args[2] else file.path(out_dir, "hpa_cache")

dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

## ---- fetch ----------------------------------------------------------------

fetch_hpa <- function(url, zip_name, tsv_name) {
  zip_path <- file.path(cache_dir, zip_name)
  tsv_path <- file.path(cache_dir, tsv_name)
  if (!file.exists(tsv_path)) {
    if (!file.exists(zip_path)) {
      message("downloading ", url)
      utils::download.file(url, zip_path, mode = "wb", quiet = TRUE)
    }
    utils::unzip(zip_path, files = tsv_name, exdir = cache_dir)
  }
  tsv_path
}

sc_path <- fetch_hpa("https://v23.proteinatlas.org/download/rna_single_cell_type.tsv.zip",
                     "rna_single_cell_type.tsv.zip", "rna_single_cell_type.tsv")
pa_path <- fetch_hpa("https://v23.proteinatlas.org/download/proteinatlas.tsv.zip",
                     "proteinatlas.tsv.zip", "proteinatlas.tsv")

sc <- data.table::fread(sc_path, showProgress = FALSE)
setnames(sc, c("Gene", "Gene name", "Cell type", "nTPM"),
             c("ensembl", "gene", "cell_type", "ntpm"))

stopifnot(TARGET %in% sc$cell_type)

## ---- gene x cell type matrix ----------------------------------------------

wide <- data.table::dcast(sc, ensembl + gene ~ cell_type, value.var = "ntpm")
key  <- wide[, .(ensembl, gene)]
mat  <- as.matrix(wide[, !c("ensembl", "gene"), with = FALSE])
stopifnot(!anyNA(mat))

n_ct   <- ncol(mat)
target <- mat[, TARGET]
other  <- mat[, setdiff(colnames(mat), TARGET), drop = FALSE]

# Row-wise descending sort of the full profile; ties are handled by comparing
# against the sorted value rather than by rank, so a gene tied at the k-th
# position still counts as inside the top-k group.
sorted <- t(apply(mat, 1, sort, decreasing = TRUE))

max_other  <- apply(other, 1, max)
mean_other <- (rowSums(mat) - target) / (n_ct - 1)
top_other  <- colnames(other)[max.col(other, ties.method = "first")]

detected   <- target >= DETECT
n_detected <- rowSums(mat >= DETECT)
rank_eec   <- rowSums(mat > target) + 1L

## ---- specificity tiers ----------------------------------------------------
# All tiers are anchored on the target cell type: a gene that is enriched in
# some other cell type is not "specific" here, whatever HPA calls it globally.

# The dominance guard matters. HPA's group-enriched and enhanced rules ask only
# whether a gene clears a bar relative to the *rest* of the body, so any gene
# restricted to gut epithelium sweeps EECs in alongside the enterocytes: ZG16 is
# ~100x higher in goblet cells, FABP2 ~50x higher in enterocytes, and EPCAM is
# pan-epithelial, yet all three clear the bar. Requiring that no other cell type
# beats EECs by more than DOM_TOL-fold holds those in the non-specific bucket and
# labels them separately as eec_shared_lineage.

is_top    <- target >= sorted[, 1]
dominated <- target * DOM_TOL < max_other

enriched <- detected & is_top & (sorted[, 1] >= FOLD * sorted[, 2])

# Group enriched: some k in 2..GROUP_MAX splits the profile so that the k-th
# highest value is >= FOLD x the (k+1)-th, and the target sits inside that top k.
in_group <- rep(FALSE, nrow(mat))
for (k in 2:min(GROUP_MAX, n_ct - 1)) {
  splits   <- sorted[, k] >= FOLD * sorted[, k + 1]
  in_group <- in_group | (detected & splits & (target >= sorted[, k]))
}

# The enhanced rule compares EECs against the *average* of the other 80 cell
# types, which a gene expressed everywhere can clear just by peaking in EECs --
# GCLC, KCTD12 and QTRT1 are each detected in 80+ of 81 cell types. Enriched and
# group-enriched genes are restricted by construction (they need a 4x gap to the
# rest), so the breadth guard is only needed here.
elevated <- detected & (target >= FOLD * mean_other) &
  (n_detected <= BREADTH * n_ct)

tier <- data.table::fcase(
  !detected,                         "not_detected",
  enriched,                          "eec_enriched",
  in_group & !dominated,             "eec_group_enriched",
  elevated & !dominated,             "eec_enhanced",
  (in_group | elevated) & dominated, "eec_shared_lineage",
  default =                          "low_specificity"
)

class3 <- data.table::fcase(
  tier == "not_detected", "not_expressed",
  tier %in% c("eec_enriched", "eec_group_enriched", "eec_enhanced"),
                          "specifically_expressed",
  default =               "non_specifically_expressed"
)

## ---- tau specificity index (on log2 scale, HPA/Yanai convention) ----------

lg    <- log2(mat + 1)
lgmax <- apply(lg, 1, max)
tau   <- ifelse(lgmax > 0, rowSums(1 - lg / lgmax) / (n_ct - 1), NA_real_)

## ---- assemble -------------------------------------------------------------

res <- tibble(
  gene                = key$gene,
  ensembl             = key$ensembl,
  ntpm_eec            = round(target, 1),
  class               = class3,
  specificity_tier    = tier,
  rank_among_celltypes = rank_eec,
  n_celltypes_detected = n_detected,
  frac_celltypes_detected = round(n_detected / n_ct, 3),
  top_other_celltype  = top_other,
  max_other_ntpm      = round(max_other, 1),
  mean_other_ntpm     = round(mean_other, 2),
  fc_vs_max_other     = round(target / pmax(max_other, 0.01), 2),
  fc_vs_mean_other    = round(target / pmax(mean_other, 0.01), 2),
  tau                 = round(tau, 3)
)

## ---- HPA's own call + secretome / GPCR annotation -------------------------

pa <- data.table::fread(pa_path, showProgress = FALSE, quote = "\"",
                        select = c("Gene", "Ensembl",
                                   "RNA single cell type specificity",
                                   "RNA single cell type specific nTPM",
                                   "Secretome location", "Protein class"))
setnames(pa, c("gene_pa", "ensembl", "hpa_sc_specificity",
               "hpa_sc_specific_ntpm", "secretome_location", "protein_class"))

res <- res %>%
  left_join(pa[, .(ensembl, hpa_sc_specificity, hpa_sc_specific_ntpm,
                   secretome_location, protein_class)],
            by = "ensembl") %>%
  mutate(is_secreted = !is.na(secretome_location) & nzchar(secretome_location))

gpcr_env <- new.env()
load(file.path("R", "sysdata.rda"), envir = gpcr_env)
res <- res %>%
  mutate(is_gpcr = gene %in% gpcr_env$gpcr_list$gene_name_primary) %>%
  select(-protein_class) %>%
  arrange(factor(class, levels = c("specifically_expressed",
                                   "non_specifically_expressed",
                                   "not_expressed")),
          desc(fc_vs_max_other), desc(ntpm_eec))

out_file <- file.path(out_dir, "eec_gene_classification.tsv")
data.table::fwrite(res, out_file, sep = "\t")

# Two-column lookup for joining onto anything keyed by gene symbol.
slim_file <- file.path(out_dir, "eec_gene_class.csv")
data.table::fwrite(res[, c("gene", "class")], slim_file)

## ---- summary --------------------------------------------------------------

cat("\nHPA v23 single-cell consensus  |  target:", TARGET,
    " |  cell types:", n_ct, " |  genes:", nrow(res), "\n\n")

cat("classification\n")
print(res %>% count(class) %>% mutate(pct = round(100 * n / sum(n), 1)))

cat("\nspecificity tier\n")
print(res %>% count(specificity_tier) %>% mutate(pct = round(100 * n / sum(n), 1)))

cat("\nsecreted genes by class\n")
print(res %>% filter(is_secreted) %>% count(class))

cat("\ntop 25 EEC-specific genes\n")
print(res %>% filter(class == "specifically_expressed") %>%
        select(gene, ntpm_eec, fc_vs_max_other, top_other_celltype,
               specificity_tier, tau, is_secreted) %>%
        head(25), n = 25)

cat("\nagreement with HPA's own single-cell specificity call\n")
print(res %>%
        mutate(hpa_names_eec = grepl("Enteroendocrine", hpa_sc_specific_ntpm, fixed = TRUE)) %>%
        count(class, hpa_sc_specificity, hpa_names_eec) %>%
        arrange(class, desc(n)), n = 40)

cat("\nwrote ", out_file, "\n", sep = "")
cat("wrote ", slim_file, "\n", sep = "")
