## =============================================================================
## Append model-score column(s) to peptide_umap_full.csv -- NO recompute of the
## params / PCA / UMAP. Just reads the CSV, joins scores on peps+target, rewrites.
##
##   pred_april17  <- April17_peps.rds$pred   (the saved April-17 model scores)
##   pred_nn       <- nn_input_comb$pred       (current session scores, if present)
##
## Join key is peps+target (C / N / loop_C / loop_N); a few duplicate keys are
## collapsed by taking the max score.
## =============================================================================

library(tidyverse)

csv_path  <- path.expand("~/AF2_analysis/peptide_umap_full.csv")
april_rds <- path.expand("~/AF2_analysis/April17_peps.rds")

read_csv_fast  <- function(p) if (requireNamespace("data.table", quietly = TRUE))
                                as_tibble(data.table::fread(p)) else readr::read_csv(p, show_col_types = FALSE)
write_csv_fast <- function(d, p) if (requireNamespace("data.table", quietly = TRUE))
                                   data.table::fwrite(d, p) else readr::write_csv(d, p)

stopifnot(file.exists(csv_path))
df  <- read_csv_fast(csv_path)
key <- intersect(c("peps", "target"), names(df))
stopifnot("peps" %in% key)
message(sprintf("CSV: %d rows, joining on %s", nrow(df), paste(key, collapse = "+")))

## collapse a scored table to one score per key (max), renamed to `new` ---------
score_by_key <- function(tbl, score_col, new) {
  tbl %>%
    dplyr::select(dplyr::all_of(c(key, score_col))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
    dplyr::summarise("{new}" := max(.data[[score_col]], na.rm = TRUE), .groups = "drop")
}

## --- April17 scores from the rds -------------------------------------------
ap  <- readRDS(april_rds)
stopifnot("pred" %in% names(ap))
df  <- dplyr::left_join(df, score_by_key(ap, "pred", "pred_april17"), by = key)
rm(ap); gc()
message(sprintf("pred_april17: matched %d / %d windows",
                sum(!is.na(df$pred_april17)), nrow(df)))

## --- current-session nn_input_comb scores (optional) -----------------------
if (exists("nn_input_comb")) {
  sc <- intersect(c("pred", "pred_raw", "pred_cal"), names(nn_input_comb))[1]
  if (!is.na(sc)) {
    df <- dplyr::left_join(df, score_by_key(nn_input_comb, sc, "pred_nn"), by = key)
    message(sprintf("pred_nn (from nn_input_comb$%s): matched %d / %d windows",
                    sc, sum(!is.na(df$pred_nn)), nrow(df)))
  } else message("nn_input_comb has no pred/pred_raw/pred_cal column -- skipping pred_nn")
} else message("nn_input_comb not in session -- skipping pred_nn")

## put the score columns up front, next to the coords -------------------------
anchor <- intersect(c("target", "terminus"), names(df))[1]
if (!is.na(anchor))
  df <- dplyr::relocate(df, dplyr::any_of(c("pred_april17", "pred_nn")),
                        .after = dplyr::all_of(anchor))

write_csv_fast(df, csv_path)
message("updated ", csv_path, "  (", nrow(df), " x ", ncol(df), ")")
