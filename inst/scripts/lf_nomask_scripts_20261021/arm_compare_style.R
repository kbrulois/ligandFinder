## One colour per arm, the same in every viz_pipeline.sh figure (sourced by the
## arm_compare_*.R plot scripts). arm_metrics_compare_plot.R keeps its own scheme
## (the arm under test red, its comparator dark grey, references light grey).
ARM_COLS <- c(curated = "#A6A6A6", flatten = "#A6A6A6", curated_AAswap = "#E08E0B",
              curated_AAswap_padg = "#7F8C8D", curated_AAswap_msacons_padg = "#2E86C1",
              msa = "#C0392B", msa_gap = "#8E44AD", curated_wAAid = "#5D6D7E", esmonly36 = "#16A085",
              esmonly36x = "#1ABC9C", motifind = "#B7950B")
## colours for `arms`, named by `labels` (default the arm names); arms without a
## fixed colour take the next free one of a colour-blind-safe set
arm_cols <- function(arms, labels = arms) {
  spare <- setdiff(c("#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#000000"), ARM_COLS)
  col <- ARM_COLS[arms]
  col[is.na(col)] <- spare[seq_len(sum(is.na(col)))]
  setNames(unname(col), labels)
}
