## For each trunk: the held-out known ends' scan scores (positives) and every
## other anchor (negatives). Train positives are DROPPED, not counted as
## negatives -- the model saw them, and scoring them highly is correct.
suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
np <- reticulate::import("numpy", convert = FALSE)
arms <- tribble(
  ~tag,                   ~term, ~label,                      ~input,
  "C_near",               "C",   "C: current",                "lf_pepend_nn_input_near.rds",
  "C_near_slack8",        "C",   "C: mask fixed",             "lf_pepend_nn_input_near.rds",
  "C_near_ctpad",         "C",   "C: CT-padded, old mask",    "lf_pepend_nn_input_near_ctpad.rds",
  "C_near_ctpad_slack8",  "C",   "C: CT-padded, mask fixed",  "lf_pepend_nn_input_near_ctpad.rds",
  "N_near",               "N",   "N: current",                "lf_pepend_nn_input_near.rds",
  "N_near_ntpad_slack8",  "N",   "N: NT-padded, mask fixed",  "lf_pepend_nn_input_near_ntpad.rds")
out <- list()
for (i in seq_len(nrow(arms))) {
  a <- arms[i, ]
  f <- file.path(A2, sprintf("lf_pepend_scan_%s.rds", a$tag))
  if (!file.exists(f)) { message("skip (no scan yet): ", a$tag); next }
  s <- readRDS(f)
  k <- readRDS(file.path(A2, a$input))$knowns %>% filter(term == a$term) %>%
    mutate(key = paste(accession, anchor))
  val <- k$key[k$split == "val"]
  sc <- s$scan %>% mutate(key = paste(accession, anchor))
  pos <- sc$score[sc$key %in% val]
  neg <- sc$score[!sc$key %in% k$key]
  stopifnot(length(pos) == length(val))
  p <- file.path(S2, sprintf("scanauc_%s.npz", a$tag))
  np$savez_compressed(p, pos = np$asarray(pos), neg = np$asarray(neg))
  message(sprintf("%-26s pos %3d  neg %9s  (dropped %d train positives)", a$label,
                  length(pos), format(length(neg), big.mark = ","), sum(k$split == "train")))
  out[[length(out) + 1L]] <- list(label = a$label, npz = p, term = a$term)
}
jsonlite::write_json(out, file.path(S2, "scan_inputs.json"), auto_unbox = TRUE)
