#!/usr/bin/env Rscript
## ---- db-anchored windows scored by BOTH window models -------------------------
## One row per dibasic-anchored window of the production set (9.2's `all` split,
## 29,599 of them for C), with per-seed score summaries from two models:
##
##   db      the production dibasic-anchored model -- read straight out of the
##           20-seed arm's member files, not recomputed
##   pepend  the peptide-end-anchored model of 10_8_pepend_train.R, run over the
##           same windows here
##
## THE TWO MODELS DO NOT AGREE ON WHERE THE ANCHOR SITS. A db-anchored window
## pins the dibasic at positions 30-31, so the canonical cleavage (cut before the
## pair, carboxypeptidase trims it) leaves the peptide ending at position 29. The
## peptide-end model was trained with that last residue at position 28. Feeding
## it a db-anchored window unchanged therefore presents the anchor one residue
## off, so this writes BOTH:
##
##   pepend_asis     the identical 36x26 window, unshifted (what "score these
##                   windows with the other model" literally means)
##   pepend_aligned  the window shifted +1 so the implied peptide end lands on
##                   position 28 -- the same cleavage hypothesis in the geometry
##                   the model was trained on. Needs the precursor's residue
##                   features, so it is NA where the window cannot be
##                   reconstructed exactly (see below).
##
## Reconstruction is self-validating: every window is rebuilt from the residue
## cache at the offset the `peps` label and the padding channel imply, and only
## windows that come back bit-exact get an aligned score. The rest are NA rather
## than silently misaligned -- ~8% of windows are clipped by the precursor
## boundary, and the two datasets do not always agree on the signal-peptide end.
##
##   Rscript inst/scripts/10_9_pepend_vs_db_windows.R
##   ... --arm ~/AF2_analysis/lf_dcnn_bench_none20/unet_none_w1 --thresh 0.2 --top-k 5
##
## Per model, per window: mean over ALL seeds, mean over the top-k seeds for THAT
## window (10_5's `pred_raw_top3` convention -- "what do the seeds that like it
## see?"), and how many seeds score it above --thresh.
##
## Output: ~/AF2_analysis/lf_pepend_vs_db_windows_C.csv
## ------------------------------------------------------------------------------

suppressMessages({ library(dplyr); library(arrow) })

.args <- commandArgs(trailingOnly = TRUE)
.opt  <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
term    <- toupper(.opt("--term", "C"))
arm     <- path.expand(.opt("--arm", "~/AF2_analysis/lf_dcnn_bench_none20/unet_none_w1"))
in_dir  <- path.expand(.opt("--in-dir", file.path(dirname(arm), "in")))
run_p   <- path.expand(.opt("--pepend-run", sprintf("~/AF2_analysis/lf_pepend_run_%s.rds", term)))
cache_p <- path.expand(.opt("--feature-cache", "~/AF2_analysis/lf_pepend_residue_cache.rds"))
out_csv <- path.expand(.opt("--out", sprintf("~/AF2_analysis/lf_pepend_vs_db_windows_%s.csv", term)))
thresh  <- as.numeric(.opt("--thresh", "0.2"))
top_k   <- as.integer(.opt("--top-k", "5"))

.this <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT  <- if (length(.this)) normalizePath(file.path(dirname(.this[[1]]), "..", "..")) else getwd()
source(file.path(ROOT, "R", "pepend_windows.R"))
source(file.path(ROOT, "R", "dcnn_bridge.R"))
mod <- lf_dcnn_python(path = file.path(ROOT, "inst", "python"))
np  <- reticulate::import("numpy", convert = FALSE)
npc <- reticulate::import("numpy", convert = TRUE)

## ---- 1. the db-anchored windows, and the db model's per-seed scores ------------
meta <- as.data.frame(read_parquet(file.path(in_dir, "meta.parquet")))
X    <- npc$load(file.path(in_dir, "arrays.npz"))[[paste0(term, "__all__x")]]
n    <- nrow(meta)
stopifnot(dim(X)[[1]] == n)
message(sprintf("%s db-anchored %s windows, %d x %d channels",
                format(n, big.mark = ","), term, dim(X)[[2]], dim(X)[[3]]))

mem <- sort(Sys.glob(file.path(arm, "members", sprintf("member_*_%s.npz", term))))
if (!length(mem)) stop("no member files under ", arm, call. = FALSE)
db <- vapply(mem, function(f) as.numeric(npc$load(f)[[paste0(term, "__global__x")]]),
             numeric(n))                                  # (n, n_seeds)
message(sprintf("db model: %d seeds read from %s", ncol(db), basename(arm)))

## ---- 2. the peptide-end model over the same windows ---------------------------
cfg_p <- file.path(paste0(tools::file_path_sans_ext(run_p), "_isolated"), "out")
cfg   <- mod$Config$from_json(file.path(dirname(cfg_p), "in", "config.json"))
models <- mod$scan$load_member_models(cfg_p, cfg, term)
message(sprintf("pepend model: %d members, class_names %s", length(models),
                paste(as.character(cfg$class_names), collapse = "/")))
stopifnot(identical(as.integer(cfg$n_channels), dim(X)[[3]]),
          identical(as.integer(cfg$seq_len), dim(X)[[2]]))

score_with <- function(models, arr) {
  x <- np$ascontiguousarray(np$asarray(arr, dtype = "float32"))
  vapply(models, function(m)
    as.numeric(reticulate::py_to_r(mod$pipeline$predict_all(m, x)[["global"]])),
    numeric(dim(arr)[[1]]))
}
message("scoring as-is ...")
pep_asis <- score_with(models, X)

## ---- 3. the +1-aligned window, where it can be rebuilt exactly ----------------
## `peps` is GENE_w<first real residue>-<last real residue>; the window's own
## position 1 is that first residue only when nothing is padded in front.
mm <- regmatches(meta$peps, regexec("^(.*)_w([0-9]+)-([0-9]+)$", meta$peps))
if (any(lengths(mm) != 4)) stop("could not parse every `peps` label", call. = FALSE)
meta$win_from <- vapply(mm, function(v) as.integer(v[[3]]), 1L)
meta$win_to   <- vapply(mm, function(v) as.integer(v[[4]]), 1L)

fc   <- readRDS(cache_p); prec <- fc$prec
stopifnot(identical(colnames(fc$feats[[1]]), fc$all_params3))
pad_col   <- which(fc$all_params3 == "padding")
first_real <- apply(X[, , pad_col], 1, function(v) { w <- which(v == 0); if (length(w)) w[[1]] else NA_integer_ })
meta$w0   <- meta$win_from - (first_real - 1L)
meta$acc  <- prec$accession[match(meta$gene, prec$gene)]

message("rebuilding windows to verify the offset ...")
jj  <- match(meta$acc, prec$accession)
ok  <- rep(FALSE, n)
shifted <- array(0, dim = dim(X))
for (i in seq_len(n)) {
  j <- jj[[i]]
  if (is.na(j) || is.na(meta$w0[[i]])) next
  f <- fc$feats[[j]]; npr <- prec$n_prot[[j]]; cpr <- prec$c_prot[[j]]
  got <- as.matrix(lf_pepend_slice(f, meta$w0[[i]], npr, cpr)[, fc$all_params3])
  if (max(abs(got - X[i, , ])) >= 1e-5) next          # cannot place this window
  ok[[i]] <- TRUE
  shifted[i, , ] <- as.matrix(lf_pepend_slice(f, meta$w0[[i]] + 1L, npr, cpr)[, fc$all_params3])
}
message(sprintf("  rebuilt bit-exact: %s of %s windows (%.1f%%); %s get an aligned score",
                format(sum(ok), big.mark = ","), format(n, big.mark = ","),
                100 * mean(ok), format(sum(ok), big.mark = ",")))

pep_al <- matrix(NA_real_, n, length(models))
if (any(ok)) {
  message("scoring +1-aligned ...")
  pep_al[ok, ] <- score_with(models, shifted[ok, , , drop = FALSE])
}

## ---- 4. summarise each model per window ---------------------------------------
## top_k: the mean of the k members that score THAT window highest, per window
summarise_model <- function(S, tag) {
  mean_all <- rowMeans(S, na.rm = FALSE)
  k <- min(top_k, ncol(S))
  mean_top <- apply(S, 1, function(v)
    if (anyNA(v)) NA_real_ else mean(sort(v, decreasing = TRUE)[seq_len(k)]))
  n_gt <- apply(S, 1, function(v) if (anyNA(v)) NA_integer_ else sum(v > thresh))
  out <- tibble(mean_all, mean_top, n_gt, sd = apply(S, 1, sd))
  names(out) <- sprintf("%s_%s", tag,
                        c("mean_allseeds", sprintf("mean_top%d", k),
                          sprintf("n_seeds_gt_%s", format(thresh)), "sd_seeds"))
  out
}

res <- bind_cols(
  meta %>% transmute(peps, gene, accession = acc, win_type, target, known,
                     win_from, win_to, win_start = w0,
                     db_pos = w0 + 29L,                     # the dibasic, precursor coords
                     pepend_anchor = w0 + 28L,              # implied peptide end
                     aligned_ok = ok),
  summarise_model(db,       "db"),
  summarise_model(pep_asis, "pepend_asis"),
  summarise_model(pep_al,   "pepend_aligned"))

res <- res %>% arrange(desc(db_mean_allseeds))
readr::write_csv(res, out_csv)
message(sprintf("\n%s rows x %d columns -> %s", format(nrow(res), big.mark = ","),
                ncol(res), out_csv))

## ---- 5. what the file says ----------------------------------------------------
message(sprintf("\nseeds: db %d, pepend %d;  threshold %.2f;  top-k %d",
                ncol(db), length(models), thresh, top_k))
num <- res %>% select(where(is.numeric)) %>% select(matches("mean_|n_seeds_"))
message("\ncolumn medians (over all rows):")
tibble(column = names(num),
       median = vapply(num, function(v) median(v, na.rm = TRUE), 1),
       n_NA   = vapply(num, function(v) sum(is.na(v)), 1L)) %>%
  as.data.frame() %>% print(row.names = FALSE, digits = 4)

k1 <- res %>% filter(known == 1)
message(sprintf("\nknown windows (n = %d) vs the rest -- median mean_allseeds:", nrow(k1)))
for (cl in c("db_mean_allseeds", "pepend_asis_mean_allseeds", "pepend_aligned_mean_allseeds"))
  message(sprintf("  %-30s known %.4f   other %.4f", cl,
                  median(k1[[cl]], na.rm = TRUE),
                  median(res[[cl]][res$known == 0], na.rm = TRUE)))

cc <- res %>% filter(aligned_ok) %>%
  summarise(db_vs_asis    = cor(db_mean_allseeds, pepend_asis_mean_allseeds, method = "spearman"),
            db_vs_aligned = cor(db_mean_allseeds, pepend_aligned_mean_allseeds, method = "spearman"),
            asis_vs_align = cor(pepend_asis_mean_allseeds, pepend_aligned_mean_allseeds,
                                method = "spearman"))
message("\nSpearman between models (rows with an aligned score):")
print(as.data.frame(cc), row.names = FALSE, digits = 3)

## how discriminating the seed-count column actually is: most windows have no
## seed over the threshold at all, so the count is only informative at the top
message(sprintf("\nwindows by how many of the 20 seeds exceed %.2f:", thresh))
cnt_cols <- grep("n_seeds_gt_", names(res), value = TRUE)
bind_rows(lapply(cnt_cols, function(cl) {
  v <- res[[cl]]
  tibble(column = cl, `0` = sum(v == 0, na.rm = TRUE), `1-4` = sum(v >= 1 & v <= 4, na.rm = TRUE),
         `5-9` = sum(v >= 5 & v <= 9, na.rm = TRUE), `10-19` = sum(v >= 10 & v <= 19, na.rm = TRUE),
         `20` = sum(v == 20, na.rm = TRUE), NA_ = sum(is.na(v)))
})) %>% as.data.frame() %>% print(row.names = FALSE)

message(sprintf("\nof the %d known windows, all %d seeds exceed %.2f in:", nrow(k1), ncol(db), thresh))
for (cl in cnt_cols)
  message(sprintf("  %-32s %d of %d", cl, sum(k1[[cl]] == ncol(db), na.rm = TRUE), nrow(k1)))
