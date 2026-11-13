## For each trunk: the held-out known ends' scan scores (positives) and every
## other anchor (negatives). Train positives are DROPPED, not counted as
## negatives -- the model saw them, and scoring them highly is correct.
suppressMessages(library(dplyr))
S2 <- commandArgs(TRUE)[1]; A2 <- path.expand("~/AF2_analysis")
## the extracted score vectors live beside the scans, not in a scratch dir: the
## scans cost ~10 min each and these are what the AUCs are computed from
OUT <- file.path(A2, "pepend_auc_inputs"); dir.create(OUT, showWarnings = FALSE)
np <- reticulate::import("numpy", convert = FALSE)
## ONLY the trunks the protein pages draw. The other arms (C mask-fixed,
## C CT-padded with the old mask) exist for the comparison and never appear on
## a page, so they are not in this figure.
arms <- tribble(
  ~tag,                   ~term, ~label,                      ~input,
  "N_near",               "N",   "N terminus",                "lf_pepend_nn_input_near.rds",
  "C_near",               "C",   "C terminus",                "lf_pepend_nn_input_near.rds",
  "N_near_ntpad_slack8",  "N",   "N terminus, NT-padded",     "lf_pepend_nn_input_near_ntpad.rds",
  "C_near_ctpad_slack8",  "C",   "C terminus, CT-padded",     "lf_pepend_nn_input_near_ctpad.rds")
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
  heads <- list(`CNN attention head` = sc$score)
  ## the per-member xgb head over the same anchors, when its full scan exists.
  ## It is written in the residues-npz order (prot_idx, anchor), which is the
  ## order lf_dcnn.scan walks too, so the two line up row for row -- checked.
  xf <- file.path(A2, sprintf("lf_winxgb_fullscan_%s.npz", a$tag))
  if (file.exists(xf)) {
    z <- np$load(xf)
    xa <- as.integer(reticulate::py_to_r(z[["anchor"]]))
    xp <- as.integer(reticulate::py_to_r(z[["prot_idx"]]))
    stopifnot(length(xa) == nrow(sc), identical(xa, as.integer(sc$anchor)))
    heads[["xgb window head"]] <- as.numeric(reticulate::py_to_r(z[["xgb_pm"]]))
  }
  for (hd in names(heads)) {
    v <- heads[[hd]]
    pos <- v[sc$key %in% val]; neg <- v[!sc$key %in% k$key]
    stopifnot(length(pos) == length(val))
    tag2 <- sprintf("%s_%s", a$tag, if (hd == "CNN attention head") "cnn" else "xgb")
    p <- file.path(OUT, sprintf("scanauc_%s.npz", tag2))
    np$savez_compressed(p, pos = np$asarray(pos), neg = np$asarray(neg))
    message(sprintf("%-26s %-20s pos %3d  neg %9s", a$label, hd,
                    length(pos), format(length(neg), big.mark = ",")))
    out[[length(out) + 1L]] <- list(label = a$label, head = hd, npz = p, term = a$term)
  }
}
jsonlite::write_json(out, file.path(OUT, "scan_inputs.json"), auto_unbox = TRUE)
message("inputs -> ", file.path(OUT, "scan_inputs.json"))
