## R -> Python -> R round trip for the lf_dcnn bridge.
##
##   /usr/local/bin/Rscript inst/python/tests/roundtrip.R
##
## Builds a small synthetic nn_input in the shape 9.2_add_contact_data.R
## produces, pushes it across the bridge, and checks the boundary: array axis
## order, the one-hot column layout, row order, and that the standalone .npz
## path reproduces the in-process one.

suppressMessages({library(tibble); library(dplyr)})   # NOT keras3: see below
source("R/dcnn_bridge.R")

ok <- 0L; bad <- 0L
check <- function(label, expr) {
  res <- tryCatch(isTRUE(expr), error = function(e) {message("  error: ", conditionMessage(e)); FALSE})
  if (res) { ok <<- ok + 1L; cat("ok   ", label, "\n") }
  else     { bad <<- bad + 1L; cat("FAIL ", label, "\n") }
}

WIN_LEN  <- 36L
CHANNELS <- c("cons_rs_n","min_afm","mean_afm","relASA",
              "SS_P","SS_S","SS_E","SS_-","SS_T","SS_G","SS_B","SS_H","SS_I",
              "Phi_cos","Psi_cos","Phi_sin","Psi_sin",
              "NH->O_1_energy","O->NH_1_energy","NH->O_2_energy","O->NH_2_energy",
              "AA_hydro","AA_charge","AA_mw","AA_pI","padding")
PAD_CH <- which(CHANNELS == "padding")

make_split <- function(n, n_pos, seed, tag) {
  set.seed(seed)
  known <- c(rep(1, n_pos), rep(0, n - n_pos))
  known_idx <- lapply(seq_len(n), function(i) {
    v <- rep("none", WIN_LEN)
    npad <- sample(0:3, 1)
    if (npad > 0) v[seq_len(npad)] <- "padding"
    v[31:32] <- "DB"
    if (known[i] == 1) { a <- sample(6:25, 1); v[a:min(a + 4, 30)] <- "pep_pocket" }
    v
  })
  data <- lapply(seq_len(n), function(i) {
    m <- matrix(runif(WIN_LEN * length(CHANNELS)), nrow = WIN_LEN,
                dimnames = list(NULL, CHANNELS))
    m[, PAD_CH] <- as.numeric(known_idx[[i]] == "padding")
    if (known[i] == 1) m[known_idx[[i]] == "pep_pocket", 2] <- m[known_idx[[i]] == "pep_pocket", 2] + 0.6
    m
  })
  ord <- sample(n)
  tibble(data = data[ord], known = known[ord], known_idx = known_idx[ord],
         peps = paste0(tag, "_", seq_len(n)))
}

nn_input <- list(
  N = list(train = make_split(96, 12, 1, "N"), val = make_split(40, 6, 2, "Nv"),
           all = make_split(60, 8, 3, "Na")),
  C = list(train = make_split(96, 12, 4, "C"), val = make_split(40, 6, 5, "Cv"),
           all = make_split(60, 8, 6, "Ca"))
)

cat("-- bridge --\n")
mod <- lf_dcnn_python()
cfg <- lf_dcnn_config(CHANNELS, epochs = 3L, patience = 2L,
                      start_from_epoch = 1L, seed = 42L, mod = mod)

check("python module imports", inherits(mod, "python.builtin.module"))
check("param count is 2438", mod$build_model(cfg)$count_params() == 2438L)

## ---- the array contract ----------------------------------------------------
arrays <- lf_dcnn_arrays(nn_input, cfg)
xa <- arrays$N$all$x
check("x is (n, win_len, n_channels)",
      identical(dim(xa), c(60L, WIN_LEN, length(CHANNELS))))
check("y_global is (n, 1)", identical(dim(arrays$N$all$y_global), c(60L, 1L)))
check("y_per_index_cat is (n, win_len, K_pi)",
      identical(dim(arrays$N$all$y_per_index_cat), c(60L, WIN_LEN, 8L)))

## axis order: R's data[[i]][pos, ch] must land at x[i, pos, ch], not transposed
probe <- vapply(1:5, function(i)
  max(abs(xa[i, , ] - nn_input$N$all$data[[i]])), numeric(1))
check("axis order survives (aperm is correct)", max(probe) == 0)

## one-hot layout: cols 1:K_cat are [6 real, none], col K_pi is the padding flag
bk  <- lf_dcnn_bookkeeping(cfg)
yc  <- arrays$N$all$y_per_index_cat
lbl <- nn_input$N$all$known_idx
check("every position is one-hot over K_pi", all(apply(yc, c(1, 2), sum) == 1))
check("padding flag sits in the LAST column",
      all((yc[, , bk$K_pi] == 1) == t(vapply(lbl, function(v) v == "padding", logical(WIN_LEN)))))
check("`none` sits at column K_cat",
      all((yc[, , bk$K_cat] == 1) == t(vapply(lbl, function(v) v == "none", logical(WIN_LEN)))))
check("pi_names order is [6 real, none, padding]",
      identical(bk$pi_names, c("CT_cleavage_context","DB","gap","NT_cleavage_context",
                               "pep_other","pep_pocket","none","padding")))
check("validate() rejects a pre-aperm array", {
  bad_arrays <- arrays
  bad_arrays$N$all$x <- aperm(xa, c(2, 3, 1))
  isTRUE(tryCatch({mod$data$as_term_data(bad_arrays, cfg); FALSE},
                  error = function(e) grepl("aperm", conditionMessage(e))))
})

## ---- in-process run --------------------------------------------------------
cat("-- training (3 epochs) --\n")
res <- lf_dcnn_run(nn_input, CHANNELS, epochs = 3L, patience = 2L,
                   start_from_epoch = 1L, seed = 42L, verbose = 0L)

n_all <- sum(vapply(nn_input, function(x) nrow(x$all), integer(1)))
check("term order follows names(nn_input)", identical(res$term_order, c("N","C")))
check("pred is (n,)", length(res$pred) == n_all)
check("pred_raw is (n,)", length(res$pred_raw) == n_all)
check("per_index is (n, win_len, K_pi)",
      identical(dim(res$per_index), c(n_all, WIN_LEN, 8L)))
check("emb is (n, 16)", identical(dim(res$emb), c(n_all, 16L)))
check("emb rows are L2-normalised", max(abs(rowSums(res$emb^2) - 1)) < 1e-6)
check("pred is a probability", all(res$pred >= 0 & res$pred <= 1))
check("pred correlates ~1.0 with rank(pred_raw)",
      abs(cor(res$pred, rank(res$pred_raw), method = "spearman") - 1) < 1e-9)
check("calibration is monotone (ranks identical)",
      identical(rank(res$pred, ties.method = "first"),
                rank(res$pred_raw, ties.method = "first")))
check("one pooled calibrator over both models' val sets",
      res$calibrator$n_obs == sum(vapply(nn_input, function(x) nrow(x$val), integer(1))))
check("per_index tibbles: one per window, win_len rows, named columns",
      length(res$per_index_tbl) == n_all &&
        nrow(res$per_index_tbl[[1]]) == WIN_LEN &&
        identical(names(res$per_index_tbl[[1]]), c("index", bk$pi_names)))

## the session contract 10_1dcnn_new6.R rebuilds for the downstream plotting
## scripts (10_2_per_ind_profiles.R, 10_2_model_importance.R, 10_3d_embed_umap.R)
check("run() returns the Config it used",
      as.integer(res$config$seq_len) == WIN_LEN &&
        as.integer(res$config$n_channels) == length(CHANNELS))
check("run() returns the built arrays as nn_in_all",
      identical(dim(res$arrays$N$all$x), c(60L, WIN_LEN, length(CHANNELS))))
## `predict()` here is deliberately un-namespaced and keras3 is NOT attached:
## this is exactly how the downstream plotting scripts call it, so it checks that
## lf_dcnn_run() left keras3's S3 methods registered.
check("keras3 generics work on the returned models without library(keras3)", {
  p <- predict(res$models[["N"]], res$arrays$N$val$x, verbose = 0)
  em <- keras3::keras_model(res$models[["N"]]$input,
                            keras3::get_layer(res$models[["N"]], "embed")$output)
  identical(names(p), c("global","per_index_cat")) &&
    identical(dim(predict(em, res$arrays$N$val$x, verbose = 0)), c(40L, 16L))
})

## row order must match bind_rows(lapply(nn_input, function(x) x$all))
comb <- bind_rows(lapply(nn_input, function(x) x$all))
check("row count matches bind_rows of the `all` splits", nrow(comb) == n_all)
check("N rows come first, then C",
      identical(comb$peps[1:3], nn_input$N$all$peps[1:3]) &&
        identical(comb$peps[n_all], nn_input$C$all$peps[nrow(nn_input$C$all)]))

## ---- standalone .npz path --------------------------------------------------
cat("-- standalone .npz path --\n")
tmp <- file.path(tempdir(), "lf_dcnn_rt")
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
meta <- comb %>% select(peps, known)
lf_dcnn_export(nn_input, CHANNELS, file.path(tmp, "in"), meta = meta,
               epochs = 3L, patience = 2L, start_from_epoch = 1L, seed = 42L)
check("export wrote arrays.npz + config.json + meta.parquet",
      all(file.exists(file.path(tmp, "in", c("arrays.npz","config.json","meta.parquet")))))

py <- file.path(dirname(reticulate::py_config()$python), "python")
rc <- system2(py, c("-m","lf_dcnn","train","--input-dir", file.path(tmp,"in"),
                    "--output-dir", file.path(tmp,"out"), "--verbose","0"),
              env = c(paste0("PYTHONPATH=", lf_dcnn_path())),
              stdout = FALSE, stderr = FALSE)
check("standalone CLI exits 0", rc == 0)

back <- lf_dcnn_import(file.path(tmp, "out"))
check("standalone pred matches the in-process pred",
      max(abs(back$pred - res$pred)) < 1e-9)
check("standalone pred_raw matches the in-process pred_raw",
      max(abs(back$pred_raw - res$pred_raw)) < 1e-9)
check("standalone emb matches the in-process emb",
      max(abs(back$emb - res$emb)) < 1e-9)
check("predictions.parquet carries the meta columns and the term",
      all(c("peps","known","term","pred","pred_raw") %in% names(back$predictions)) &&
        nrow(back$predictions) == n_all)

cat(sprintf("\n%d/%d checks passed\n", ok, ok + bad))
if (bad > 0) quit(status = 1)
