#!/usr/bin/env Rscript
## Does the Python port reproduce the R model?
##
##   Rscript inst/python/tests/compare_r_python.R                     # 1 seed, full run
##   Rscript inst/python/tests/compare_r_python.R --seeds 1,2,3
##   Rscript inst/python/tests/compare_r_python.R --seeds 42 --epochs 300 --arms r,flat
##
## Trains the reference R model and the Python port on the SAME arrays with the
## same seed, and reports validation PR-AUC / ROC-AUC side by side.
##
## Read the spread, not a single number. With only 6-7 validation positives per
## terminus, PR-AUC swings ~0.4 across seeds within EITHER implementation, so a
## single-seed difference means nothing. Run >=3 seeds and compare the columns.
##
## Arms:
##   r      the frozen R reference model below (keras3 via reticulate)
##   flat   the port, cfg$trunk = "flat"   -- should match `r` within noise
##   unet   the port, cfg$trunk = "unet"   -- the pooling encoder/decoder trunk

suppressMessages({library(keras3); library(tensorflow); library(tfdatasets); library(dplyr)})

## ---- args ------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) {
  i <- match(flag, .args); if (is.na(i) || i == length(.args)) default else .args[[i + 1L]]
}
SEEDS  <- as.integer(strsplit(.opt("--seeds", "42"), ",")[[1]])
EPOCHS <- as.integer(.opt("--epochs", "2000"))
ARMS   <- strsplit(.opt("--arms", "r,flat,unet"), ",")[[1]]
CACHE  <- .opt("--cache", "~/AF2_analysis/lf_dcnn_compare_nn_input.rds")
ROOT   <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."))

message(sprintf("seeds=%s epochs=%d arms=%s", paste(SEEDS, collapse=","), EPOCHS, paste(ARMS, collapse=",")))

## ---- 1. the window data ----------------------------------------------------
## Cached because 9.2 takes minutes even off its own cache.
if (file.exists(path.expand(CACHE))) {
  message("loading cached nn_input from ", CACHE)
  .cc <- readRDS(path.expand(CACHE)); nn_input <- .cc$nn_input; all_params3 <- .cc$all_params3
} else {
  message("building nn_input via 9.2_add_contact_data.R (needs ~/AF2_analysis/nn_dat_cache.rds) ...")
  suppressMessages(source(file.path(ROOT, "inst/scripts/putative_peptide/generate_residue_db_untested/9.2_add_contact_data.R")))
  saveRDS(list(nn_input = nn_input, all_params3 = all_params3), path.expand(CACHE))
  message("cached to ", CACHE)
}
for (t in names(nn_input)) for (s in c("train","val"))
  message(sprintf("  %s/%-5s n=%4d  positives=%3d", t, s,
                  nrow(nn_input[[t]][[s]]), sum(nn_input[[t]][[s]]$known)))

## ---- 2. one set of arrays, shared by every arm ------------------------------
source(file.path(ROOT, "R/dcnn_bridge.R"))
mod  <- lf_dcnn_python(path = file.path(ROOT, "inst/python"))
## r_exact: reproduce the pre-port R model, including its DB-mask overwrite and
## its non-functional per-epoch resample. Without this the arms are not comparable.
cfg0 <- lf_dcnn_config(all_params3, r_exact = TRUE, seed = SEEDS[[1]], mod = mod)
arrays <- lf_dcnn_arrays(nn_input, cfg0, splits = c("train", "val"))

## ---- 3. the frozen R reference model ---------------------------------------
## Lifted verbatim from 10_1dcnn_new6.R before the port. Do not "improve" it --
## its only job is to be the thing the port is compared against.
r_reference_fit <- function(nn_in, term, all_params3, seed, epochs) {
  keras3::clear_session(); keras3::set_random_seed(seed); set.seed(seed)

  classes <- setNames(0:7, c("CT_cleavage_context","DB","gap","NT_cleavage_context",
                             "pep_other","pep_pocket","padding","none"))
  seq_len <- 36L; n_channels <- length(all_params3); K <- 8L
  mid_mask <- CT_mask <- NT_mask <- all_mask <- rep(0, seq_len)
  NT_mask[1:5] <- 1; CT_mask[31:36] <- 1; mid_mask[6:30] <- 1
  mask_matrix <- matrix(0, nrow = seq_len, ncol = K)
  mask_matrix[, which(names(classes) %in% c("DB","NT_cleavage_context"))] <- NT_mask
  mask_matrix[, which(names(classes) %in% c("DB","CT_cleavage_context"))] <- CT_mask
  mask_matrix[, which(names(classes) %in% c("pep_other","pep_pocket","gap"))] <- mid_mask
  mask_matrix[, which(names(classes) %in% c("padding"))] <- all_mask
  real_cols <- which(!names(classes) %in% c("none","padding"))
  cat_cols  <- c(real_cols, which(names(classes) == "none")); K_cat <- length(cat_cols)
  mask_matrix_cat <- cbind(mask_matrix[, real_cols, drop = FALSE], none = 1)
  pi_cols  <- c(cat_cols, which(names(classes) == "padding")); K_pi <- length(pi_cols)
  pi_names <- names(classes)[pi_cols]
  pi_w <- setNames(rep(1, K_pi), pi_names); pi_w["DB"] <- 2
  pi_w[if (term == "C") "CT_cleavage_context" else "NT_cleavage_context"] <- 3
  pi_w <- pi_w / mean(pi_w)
  pad_channel_id <- which(all_params3 == "padding")
  cont_channels  <- which(all_params3 %in% c("cons_rs","cons_rs_n","min_afm","mean_afm","relASA"))

  make_ds <- function(noise_frac = 0.1, n_neg = 3, batch_size = 32) {
    y <- nn_in$train$y_global[, 1]
    pos_idx <- which(y == 1); neg_idx <- which(y == 0)
    n_neg_target <- length(pos_idx) * n_neg
    sampled_neg <- sample(neg_idx, n_neg_target, replace = length(neg_idx) < n_neg_target)
    idx <- sample(c(pos_idx, sampled_neg))
    x_sel <- nn_in$train$x[idx, , , drop = FALSE]
    if (noise_frac > 0) {
      real_mask <- x_sel[, , pad_channel_id] == 0
      for (ch in cont_channels) {
        sl <- x_sel[, , ch]; sd_ch <- stats::sd(sl[real_mask], na.rm = TRUE)
        if (is.finite(sd_ch) && sd_ch > 0)
          x_sel[, , ch] <- sl + matrix(rnorm(length(sl), 0, noise_frac * sd_ch), nrow = nrow(sl)) * real_mask
      }
    }
    tensor_slices_dataset(list(
      tensorflow::as_tensor(x_sel, dtype = "float32"),
      tensorflow::as_tensor(nn_in$train$y_global[idx, , drop = FALSE], dtype = "float32"),
      tensorflow::as_tensor(nn_in$train$y_per_index_cat[idx, , ], dtype = "float32"))) |>
      dataset_shuffle(buffer_size = length(idx)) |> dataset_batch(batch_size) |>
      dataset_map(function(x, yg, yc) list(x, list(global = yg, per_index_cat = yc))) |>
      dataset_repeat()
  }

  inputs <- layer_input(shape = c(seq_len, n_channels))
  regularizer <- regularizer_l2(1e-3)
  pos_channel <- inputs |> layer_lambda(function(x) {
    ones <- op_sum(x, axis = -1L, keepdims = TRUE) * 0 + 1
    (op_cumsum(ones, axis = 2L) - 1) / (seq_len - 1)
  }, output_shape = c(seq_len, 1L))
  shared <- layer_concatenate(list(inputs, pos_channel)) |>
    layer_conv_1d(16, kernel_size = 3, activation = "gelu",
                  kernel_regularizer = regularizer, padding = "same") |>
    layer_dropout(0.2) |>
    layer_conv_1d(8, kernel_size = 3, activation = "gelu",
                  kernel_regularizer = regularizer, padding = "same") |>
    layer_dropout(0.3)
  per_index_logits <- shared |>
    layer_conv_1d(filters = K_pi, kernel_size = 1, kernel_regularizer = regularizer)
  per_index_cat <- per_index_logits |> layer_activation("softmax", name = "per_index_cat")
  masked_sum <- per_index_logits |>
    layer_lambda(function(x) op_take(x, indices = as.integer(0:(K_cat - 1L)), axis = -1L),
                 output_shape = c(seq_len, K_cat)) |>
    layer_lambda(function(x) x * op_convert_to_tensor(mask_matrix_cat, dtype = "float32"),
                 output_shape = c(seq_len, K_cat)) |>
    layer_activation("softmax")
  global_output <- layer_multi_head_attention(num_heads = 2, key_dim = 8)(masked_sum, masked_sum) |>
    layer_layer_normalization() |> layer_global_average_pooling_1d() |>
    layer_dense(16, activation = "gelu", name = "embed") |>
    layer_dropout(0.3) |> layer_dense(1, activation = "sigmoid", name = "global")
  model <- keras_model(inputs = inputs,
                       outputs = list(global = global_output, per_index_cat = per_index_cat))

  per_index_cat_loss <- function(y_true, y_pred, gamma = 2, smoothness_weight = 0.01) {
    eps <- 1e-7
    labels <- y_true[, , 1:K_pi]; none_flag <- y_true[, , K_cat]; keep <- 1 - none_flag
    p  <- op_clip(y_pred, eps, 1 - eps)
    cw <- op_reshape(op_convert_to_tensor(as.numeric(pi_w), dtype = "float32"), c(1L,1L,K_pi))
    ce   <- -op_sum(cw * labels * op_power(1 - p, gamma) * op_log(p), axis = -1L)
    main <- op_sum(ce * keep, axis = 2L) / (op_sum(keep, axis = 2L) + eps)
    sm_mask <- keep[, 2:seq_len] * keep[, 1:(seq_len - 1)]
    sq_pos  <- op_mean(op_square(y_pred[, 2:seq_len, ] - y_pred[, 1:(seq_len - 1), ]), axis = 3L)
    main + smoothness_weight * (op_sum(sq_pos * sm_mask, axis = 2L) / (op_sum(sm_mask, axis = 2L) + eps))
  }
  masked_cat_accuracy <- custom_metric("masked_cat_accuracy", function(y_true, y_pred) {
    labels <- y_true[, , 1:K_pi]; keep <- (1 - y_true[, , K_pi]) * (1 - y_true[, , K_cat])
    correct <- op_cast(op_equal(op_argmax(y_pred, axis = -1L),
                                op_argmax(labels, axis = -1L)), "float32") * keep
    op_sum(correct) / (op_sum(keep) + 1e-7)
  })
  wbce <- function(w0, w1) function(y_true, y_pred)
    -(w1 * y_true * log(y_pred) + w0 * (1 - y_true) * log(1 - y_pred))

  model |> compile(
    optimizer = optimizer_adam(learning_rate = 3e-4, clipnorm = 1.0),
    loss = list(global = wbce(1, 1), per_index_cat = per_index_cat_loss),
    loss_weights = list(global = 0.05, per_index_cat = 1),
    metrics = list(global = list(metric_auc(name = "auc"), metric_auc(name = "pr_auc", curve = "PR")),
                   per_index_cat = list(masked_cat_accuracy)))

  train_ds <- make_ds()
  n_pos <- sum(nn_in$train$y_global)
  h <- model |> fit(train_ds,
    validation_data = list(nn_in$val$x, list(global = nn_in$val$y_global,
                                             per_index_cat = nn_in$val$y_per_index_cat)),
    steps_per_epoch = ceiling(n_pos * 4 / 32), epochs = epochs,
    callbacks = list(callback_early_stopping(monitor = "val_global_pr_auc", patience = 300,
                                             mode = "max", start_from_epoch = 300,
                                             restore_best_weights = TRUE),
                     callback_lambda(on_epoch_begin = function(epoch, logs) {
                       train_ds <<- make_ds() })),
    batch_size = NULL, verbose = 0)
  list(pr = max(h$metrics$val_global_pr_auc), auc = max(h$metrics$val_global_auc),
       epochs = length(h$metrics$val_global_pr_auc), params = model$count_params())
}

## ---- 4. the python arms -----------------------------------------------------
py_fit <- function(term, seed, trunk, epochs) {
  cfg <- lf_dcnn_config(all_params3, r_exact = TRUE, seed = seed,
                        epochs = epochs, trunk = trunk, mod = mod)
  cfg <- lf_dcnn_align_terms(cfg, names(nn_input))
  td  <- mod$data$as_term_data(arrays, cfg)
  mod$pipeline$set_seed(cfg$seed)
  res <- mod$train(term, td[[term]], cfg, verbose = 0L,
                   rng = reticulate::import("numpy")$random$default_rng(seed))
  h <- res[[2]]
  list(pr = max(unlist(h$val_global_pr_auc)), auc = max(unlist(h$val_global_auc)),
       epochs = length(unlist(h$val_global_pr_auc)), params = res[[1]]$count_params())
}

## ---- 5. run every arm on identical arrays -----------------------------------
rows <- list()
for (seed in SEEDS) for (term in names(nn_input)) for (arm in ARMS) {
  t0 <- Sys.time()
  r <- if (arm == "r") r_reference_fit(arrays[[term]], term, all_params3, seed, EPOCHS)
       else py_fit(term, seed, arm, EPOCHS)
  rows[[length(rows) + 1L]] <- data.frame(
    term = term, seed = seed, arm = arm, pr_auc = r$pr, roc_auc = r$auc,
    epochs = r$epochs, params = r$params,
    mins = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  message(sprintf("  [%s] %-4s seed=%-3d pr_auc=%.4f roc_auc=%.4f (%d epochs, %.1f min)",
                  term, arm, seed, r$pr, r$auc, r$epochs, tail(rows,1)[[1]]$mins))
}
res <- dplyr::bind_rows(rows)

cat("\n================ per run ================\n"); print(res, row.names = FALSE)
cat("\n================ summary ================\n")
print(res %>% group_by(term, arm) %>%
        summarise(n = n(), params = first(params),
                  pr_mean = round(mean(pr_auc), 3), pr_sd = round(sd(pr_auc), 3),
                  pr_range = round(diff(range(pr_auc)), 3),
                  roc_mean = round(mean(roc_auc), 4), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\nInterpretation: compare pr_mean BETWEEN arms against pr_sd/pr_range WITHIN\n",
    "an arm. If the between-arm gap is smaller than the within-arm spread, the\n",
    "implementations agree. roc_auc is the stabler metric at these positive counts.\n")

out <- file.path(path.expand("~/AF2_analysis"), sprintf("lf_dcnn_compare_%s.csv", format(Sys.Date(), "%Y%m%d")))
write.csv(res, out, row.names = FALSE); cat("\nwrote", out, "\n")
