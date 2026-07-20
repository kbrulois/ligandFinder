
masked_focal_loss <- function(y_true, y_pred, gamma = 2.0, alpha = 0.25) {
  # 1. Use standard CE as base (shape: batch, 36)
  ce <- loss_categorical_crossentropy(y_true, y_pred, from_logits = FALSE, label_smoothing = 0.1)

  # 2. Compute focal factor: (1 - p)^gamma
  # This makes 'easy' predictions contribute almost zero to the loss
  p_t <- op_sum(y_true * y_pred, axis = -1L)
  focal_factor <- op_power(1.0 - p_t, gamma)

  # 3. Apply the mask (ignore class 7)
  class_ids <- op_argmax(y_true, axis = -1L)
  mask <- op_cast(op_not_equal(class_ids, 8L), "float32")

  # 4. Multiply and return
  return(op_sum(ce * focal_factor * mask, axis = -1L))
}

masked_focal_loss <- function(y_true, y_pred, gamma = 2.0) {

  ce <- loss_categorical_crossentropy(
    y_true, y_pred,
    from_logits = FALSE,
    label_smoothing = 0.1
  )

  p_t <- op_sum(y_true * y_pred, axis = -1L)
  focal_factor <- op_power(1.0 - p_t, gamma)

  loss <- ce * focal_factor

  # Ensure shape = (batch,)
  return(loss)
}

masked_focal_loss <- function(y_true, y_pred, gamma = 2.0) {
  ce <- loss_categorical_crossentropy(y_true, y_pred, from_logits = FALSE, label_smoothing = 0.05)
  p_t <- op_sum(y_true * y_pred, axis = -1L)
  focal_factor <- op_power(1.0 - p_t, gamma)

  # mask out positions where true class is the "none" class (index 0 or whichever you removed)
  class_ids <- op_argmax(y_true, axis = -1L)
  mask <- op_cast(op_not_equal(class_ids, 0L), "float32")  # adjust index to match your "none" class position

  loss <- ce * focal_factor * mask
  op_sum(loss, axis = 2L) / (op_sum(mask, axis = 2L) + 1e-7)  # masked mean not sum
}

masked_focal_loss <- function(y_true, y_pred, gamma = 2.0, alpha = 0.25) {
  # 1. Use standard CE as base (shape: batch, 36)
  ce <- loss_categorical_crossentropy(y_true, y_pred, from_logits = FALSE, label_smoothing = 0.1)

  # 2. Compute focal factor: (1 - p)^gamma
  # This makes 'easy' predictions contribute almost zero to the loss
  p_t <- op_sum(y_true * y_pred, axis = -1L)
  focal_factor <- op_power(1.0 - p_t, gamma)

  # 3. Apply the mask (ignore class 7)
  class_ids <- op_argmax(y_true, axis = -1L)
  mask <- op_cast(op_not_equal(class_ids, 4L), "float32")

  # 4. Multiply and return
  return(op_sum(ce * focal_factor * mask, axis = 2L))
}

masked_focal_loss <- function(y_true, y_pred, gamma = 2.0) {
  ce <- loss_categorical_crossentropy(y_true, y_pred, from_logits = FALSE, label_smoothing = 0.05)
  p_t <- op_sum(y_true * y_pred, axis = -1L)
  focal_factor <- op_power(1.0 - p_t, gamma)

  # mask out positions where true class is the "none" class (index 0 or whichever you removed)
  class_ids <- op_argmax(y_true, axis = -1L)
  mask <- op_cast(op_not_equal(class_ids, 8L), "float32")  # adjust index to match your "none" class position

  loss <- ce * focal_factor * mask
  op_sum(loss, axis = 2L) / (op_sum(mask, axis = 2L) + 1e-7)  # masked mean not sum
}

masked_focal_loss <- function(y_true, y_pred, gamma = 2.0) {

  ce <- loss_categorical_crossentropy(
    y_true, y_pred,
    from_logits = FALSE,
    label_smoothing = 0.15
  )                        # (batch, seq)

  p_t <- op_sum(y_true * y_pred, axis = -1L)   # (batch, seq)

  focal_factor <- op_power(1.0 - p_t, gamma)

  class_ids <- op_argmax(y_true, axis = -1L)

  mask <- op_cast(op_not_equal(class_ids, 8L), "float32")

  loss <- ce * focal_factor * mask

  op_sum(loss, axis = 2L) / (op_sum(mask, axis = 2L) + 1e-7)
}



per_index_loss_fn <- function(y_true, y_pred) {

  pos_weight <- 10

  eps <- 1e-7

  y_pred <- tf$clip_by_value(y_pred, eps, 1 - eps)

  loss <- -(pos_weight * y_true * tf$math$log(y_pred) +
              (1 - y_true) * tf$math$log(1 - y_pred))

  tf$reduce_mean(loss)
}


masked_categorical_loss <- function(y_true, y_pred) {
  # 1. FALSE because your model has a softmax layer at the end
  ce <- loss_categorical_crossentropy(y_true, y_pred, from_logits = FALSE, label_smoothing = 0.1)

  class_ids <- op_argmax(y_true, axis = -1L)
  mask <- op_cast(op_not_equal(class_ids, 8L), "float32")

  # 2. INCREASE THE SIGNAL: multiply the loss of non-7 classes by 5 or 10
  # This forces the model to treat mistakes on real classes as "expensive"
  boosted_loss <- ce * mask * 100.0

  return(op_sum(boosted_loss, axis = -1L))
}

masked_accuracy <- function(y_true, y_pred) {
  # 1. Get predicted class and actual class indices
  y_pred_idx <- op_argmax(y_pred, axis = -1L)
  y_true_idx <- op_argmax(y_true, axis = -1L)

  # 2. Check where they match (Boolean)
  correct_preds <- op_equal(y_true_idx, y_pred_idx)

  # 3. Create mask for everything EXCEPT class 7
  mask <- op_cast(op_not_equal(y_true_idx, 7L), "float32")

  # 4. Multiply matches by mask (zeros out matches on class 7)
  weighted_correct <- op_cast(correct_preds, "float32") * mask

  # 5. Return average: (Correct non-7) / (Total non-7)
  return(op_sum(weighted_correct) / (op_sum(mask) + 1e-7))
}



weighted_binary_crossentropy <- function(weight_0, weight_1) {
  function(y_true, y_pred) {
    - (weight_1 * y_true * log(y_pred) + weight_0 * (1 - y_true) * log(1 - y_pred))
  }
}



library(keras3)
library(tensorflow)
library(yardstick)



class_cols <- setNames(c("#FED439FF", "#370335FF", "#8A9197FF", "#D2AF81FF",
                         "#D5E4A2FF", "#197EC0FF", "grey85", "#075149FF"),
                       c("CT_cleavage_context", "DB", "gap", "NT_cleavage_context", "pep_other", "pep_pocket", "padding", "none")
)


classes <- lapply(nn_input, function(x) x$train$known_idx) %>% do.call(`c`, .) %>% do.call(`c`, .) %>% unique
classes <- names(class_cols)
classes <- setNames(0:(length(classes) - 1), classes)


seq_len <- nn_input[[1]]$train$data %>% `[[`(1) %>% nrow
n_channels <- nn_input[[1]]$train$data %>% `[[`(1) %>% ncol
K <- lapply(nn_input, function(x) x$train$known_idx) %>% do.call(`c`, .) %>% do.call(`c`, .) %>% unique %>% length
kernel_size <- 8

mid_mask <- CT_mask <- NT_mask <- all_mask <- rep(0, seq_len)
NT_mask[1:5] <- 1
CT_mask[31:36] <- 1

mid_mask[6:30] <- 1

mask_matrix <- matrix(0, nrow = seq_len, ncol = K)

mask_matrix[, which(names(classes) %in% c("DB", "NT_cleavage_context"))] <- NT_mask
mask_matrix[, which(names(classes) %in% c("DB", "CT_cleavage_context"))] <- CT_mask
mask_matrix[, which(names(classes) %in% c("pep_other", "pep_pocket", "gap"))] <- mid_mask
mask_matrix[, which(names(classes) %in% c("padding"))] <- all_mask

# --- per-index multi-label head: real classes only ------------------------
# The per-index head predicts an INDEPENDENT sigmoid per real class. `none` is
# not a class (background = all sigmoids low), and `padding` is not a class
# either: padding positions are masked out of the loss and zeroed in the output
# using the known input padding channel. The global head uses a softmax over the
# same real classes, so `none`/`padding` never take softmax mass there either.
real_cols        <- which(!names(classes) %in% c("none", "padding"))   # 1-based real-class columns
K_real           <- length(real_cols)

## The per-index head is a single-label SOFTMAX over the real classes PLUS an
## explicit "none" (background) class; padding is still excluded (masked in loss).
none_col <- which(names(classes) == "none")
cat_cols <- c(real_cols, none_col)                                     # 6 real + none (none last)
K_cat    <- length(cat_cols)

## per-class weights for the (softmax categorical) per-index loss. Up-weight the
## classes that matter (DB, pep_pocket, pep_other), down-weight positionally-trivial
## contexts and background `gap`, and down-weight the abundant `none` so it doesn't
## dominate the softmax. Reordered to the output-column order names(classes)[cat_cols].
class_w <- c(CT_cleavage_context = 1, DB = 5, gap = 0.3, NT_cleavage_context = 1,
             pep_other = 2, pep_pocket = 5, none = 0.3)[names(classes)[cat_cols]]
stopifnot(!anyNA(class_w))   # every class (incl none) must have a weight
class_w <- class_w / mean(class_w)   # mean 1: relative emphasis, no magnitude inflation
mask_matrix_real <- mask_matrix[, real_cols, drop = FALSE]             # (seq_len, K_real) position mask
mask_matrix_cat  <- cbind(mask_matrix_real, none = 1)                  # (seq_len, K_cat): none allowed at all positions
padding_col      <- which(names(classes) == "padding")                # 1-based padding column (target flag)
pad_channel_id   <- if (exists("all_params3")) which(all_params3 == "padding") else n_channels
stopifnot(length(pad_channel_id) == 1L)                               # padding must be a single input channel
channel_onehot   <- as.numeric(seq_len(n_channels) == pad_channel_id) # (C,) selects the pad input channel

## continuous channels that receive train-time noise augmentation. The one-hot /
## flag channels (AA_*, SS_*, padding, end_type_ch) are deliberately left intact
## so augmented inputs stay on the manifold.
cont_channels <- if (exists("all_params3"))
  which(all_params3 %in% c("cons_rs", "cons_rs_n", "min_afm", "mean_afm", "relASA")) else integer(0)

## noise-augmentation strength (fraction of each continuous channel's sd). Set to
## 0 to fall back to plain exact-repeat oversampling for an A/B comparison.
aug_noise_frac <- 0.1

library(tfdatasets)

make_oversampled_dataset <- function(nn_in, n_negatives_per_positive = 3, batch_size = 32,
                                     noise_frac = 0.1) {

  y <- nn_in[["train"]][["y_global"]][, 1]
  pos_idx <- which(y == 1)
  neg_idx <- which(y == 0)

  x_train <- nn_in[["train"]][["x"]]
  y_global <- nn_in[["train"]][["y_global"]]
  y_cat <- nn_in[["train"]][["y_per_index_cat"]]

  # repeat positives to match desired ratio
  n_pos <- length(pos_idx)
  n_neg_target <- n_pos * n_negatives_per_positive

  sampled_neg <- sample(neg_idx, n_neg_target, replace = length(neg_idx) < n_neg_target)
  idx <- sample(c(rep(pos_idx, 1), sampled_neg))  # shuffle

  ## --- input-noise augmentation ---------------------------------------------
  ## Jitter the continuous channels at REAL (non-padding) positions so the tiny,
  ## repeated positive set can't be memorized. One-hot/flag channels are left
  ## untouched, padding positions stay 0, and val is never touched (train only).
  ## Noise is per-channel (scaled to that channel's sd) and regenerated every
  ## epoch (resample_callback rebuilds the dataset). Set noise_frac = 0 to disable.
  x_sel <- x_train[idx, , , drop = FALSE]                       # (n, seq, n_channels)
  if (noise_frac > 0 && length(cont_channels) > 0) {
    real_mask <- x_sel[, , pad_channel_id] == 0                 # (n, seq): TRUE at real positions
    for (ch in cont_channels) {
      sl    <- x_sel[, , ch]                                    # (n, seq)
      sd_ch <- stats::sd(sl[real_mask], na.rm = TRUE)
      if (is.finite(sd_ch) && sd_ch > 0) {
        noise <- matrix(rnorm(length(sl), 0, noise_frac * sd_ch), nrow = nrow(sl))
        x_sel[, , ch] <- sl + noise * real_mask                # add only at real positions
      }
    }
  }

  x_t <- tensorflow::as_tensor(x_sel, dtype = "float32")
  yg_t <- tensorflow::as_tensor(y_global[idx,, drop = FALSE], dtype = "float32")

  tensor_slices_dataset(list(x_t, yg_t)) |>       # single-output model: global target only
    dataset_shuffle(buffer_size = length(idx)) |>
    dataset_batch(batch_size) |>
    dataset_map(function(x, yg) {
      list(x, list(global = yg))
    }) |>
    dataset_repeat()  # repeat so keras can run multiple epochs
}
models <- list()

val_pred <- list()

for(term in names(nn_input)) {

  rm(model,per_index_logits, per_index_cat, shared,  per_index_logits, per_index_output, per_index_output2, global_output, early_stopping)
  keras3::clear_session()

  generate_keras_input <- function(nn_in = nn_input[[term]]) {

    Map(function(input) {

      n_samples <- length(input[["data"]])

      x <- array(
        unlist(input[["data"]]),
        dim = c(seq_len, n_channels, length(input[["data"]]))
      )

      x <- aperm(x, c(3, 1, 2))

      y_global <- matrix(input[["known"]], ncol = 1)

      onehot <- array(input[["known_idx"]] %>%
                        do.call(`c`, .) %>%
                        classes[.] %>%
                        to_categorical(., num_classes = K),
                      dim = c(seq_len, n_samples, K)) %>%
        aperm(., c(2, 1, 3))                                    # (n, seq, K) one-hot over all classes

      # Single-label target: one-hot over the K_cat classes (6 real + explicit
      # none); the last column carries a padding flag so the loss/metric can mask
      # padding positions.
      y_per_index_cat <- array(0, dim = c(n_samples, seq_len, K_cat + 1L))
      y_per_index_cat[, , 1:K_cat]    <- onehot[, , cat_cols]
      y_per_index_cat[, , K_cat + 1L] <- onehot[, , padding_col]


      list(x = x,
           y_global = y_global,
           y_per_index_cat = y_per_index_cat)

    }, nn_in)

  }

  nn_in <- generate_keras_input(nn_in = nn_input[[term]])


  inputs <- layer_input(shape = c(seq_len, n_channels))

  regularizer <- regularizer_l2(1e-3)

  # Positional encoding: the windows are anchored/aligned (fixed cleavage/DB
  # layout), so absolute position is meaningful, but conv is translation-
  # equivariant and the attention/pooling are order-agnostic. Append a normalized
  # position ramp [0,1] as an extra input channel so the trunk (and everything
  # above it) can learn where along the window it is. Built in-graph so the raw
  # `inputs` (and its padding-channel logic below) stays untouched.
  pos_channel <- inputs |>
    layer_lambda(function(x) {
      ones <- op_sum(x, axis = -1L, keepdims = TRUE) * 0 + 1  # (batch, seq, 1) ones
      idx  <- op_cumsum(ones, axis = 2L) - 1                  # (batch, seq, 1): 0..seq-1
      idx / (seq_len - 1)                                     # normalized ramp [0,1]
    }, output_shape = c(seq_len, 1L))

  conv_in <- layer_concatenate(list(inputs, pos_channel))               # (batch, seq, n_channels + 1)

  # Slim stack for the tiny (<100 example) training sets: two small conv layers
  # instead of the old four wide ones. The old k=16 / k=8, 32-filter layers held
  # ~26k of the ~33k params; this brings the conv trunk to ~4k, far more
  # defensible given the 4-way split leaves few positives per model.
  shared <- conv_in |>
    layer_conv_1d(16, kernel_size = 3, activation = "gelu",
                  kernel_regularizer = regularizer, padding = "same") |>
    layer_dropout(0.2) |>
    layer_conv_1d(8, kernel_size = 3, activation = "gelu",
                  kernel_regularizer = regularizer, padding = "same") |>
    layer_dropout(0.3)

  per_index_logits <- shared |>
    layer_conv_1d(filters = K_cat, kernel_size = 1,
                  kernel_regularizer = regularizer)

  # Per-index output kept as an UNSUPERVISED head: softmax over the K_cat classes
  # (6 real + none), emitted for visualization only. It has NO loss in compile, so
  # it doesn't compete with / shape training -- per_index_logits is shaped solely
  # by the global head via masked_sum below.
  per_index_cat <- per_index_logits |>
    layer_activation("softmax", name = "per_index_cat")

  # Global head input: masked softmax over ALL K_cat classes (6 real + none), so
  # the ranking head's class-structure view includes the background/none level.
  masked_sum <- per_index_logits |>
    layer_lambda(function(x) {
      mask <- op_convert_to_tensor(mask_matrix_cat, dtype = "float32")
      x * mask
    }, output_shape = c(seq_len, K_cat)) |>
    layer_activation("softmax")

  if(FALSE) {
    global_output <- masked_sum |>
      layer_global_average_pooling_1d() |>
      layer_dense(16, activation = "gelu") |>
      layer_dense(1, activation = "sigmoid", name = "global")
  }

  # Global head: a SMALL multi-head attention (2 heads, key_dim 8, ~0.5k params)
  # over the masked class softmax, then pool. A plain GAP of the softmax washed
  # out all positional signal and tanked global AUC; the attention restores it at
  # ~1/6 the cost of the old 4-head/key_dim-32 version. The dense named "embed"
  # is the 16-d representation reused (no new loss) for nearest-known retrieval.
  global_output <- layer_multi_head_attention(num_heads = 2, key_dim = 8)(masked_sum, masked_sum) |>
    layer_layer_normalization() |>
    layer_global_average_pooling_1d() |>
    layer_dense(16, activation = "gelu", name = "embed") |>
    layer_dropout(0.3) |>
    layer_dense(1, activation = "sigmoid", name = "global")



  lr_reduce <- callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 20, min_lr = 1e-6)

  early_stopping <- callback_early_stopping(
    monitor = "val_global_pr_auc",
    patience = 200,
    mode = 'max',
    start_from_epoch = 30,
    restore_best_weights = TRUE
  )


  model <- keras_model(
    inputs = inputs,
    outputs = list(global = global_output,
                   per_index_cat = per_index_cat)   # per_index_cat is UNSUPERVISED (no loss in compile)
  )

  message(sprintf("[%s] trainable params: %s", term,
                  format(model$count_params(), big.mark = ",")))



  # Categorical focal cross-entropy over the K_cat classes (6 real + none).
  # y_true is packed: columns 1:K_cat are the one-hot class (incl none), and the
  # last column is a padding flag (1 = padding) used to mask those positions.
  per_index_cat_loss <- function(y_true, y_pred, gamma = 2, smoothness_weight = 0.01) {

    eps      <- 1e-7
    labels   <- y_true[, , 1:K_cat]         # (batch, seq, K_cat) one-hot incl none
    pad_flag <- y_true[, , K_cat + 1L]      # (batch, seq): 1 at padding
    keep     <- 1 - pad_flag                # (batch, seq): 0 at padding, 1 elsewhere

    p  <- op_clip(y_pred, eps, 1 - eps)     # softmax probs (batch, seq, K_cat)
    cw <- op_reshape(op_convert_to_tensor(as.numeric(class_w), dtype = "float32"),
                     c(1L, 1L, K_cat))      # (1, 1, K_cat) per-class weights

    # focal categorical CE per position: -sum_c w_c * y_c * (1-p_c)^gamma * log(p_c)
    ce   <- -op_sum(cw * labels * op_power(1 - p, gamma) * op_log(p), axis = -1L)  # (batch, seq)
    main <- op_sum(ce * keep, axis = 2L) / (op_sum(keep, axis = 2L) + eps)          # (batch,)

    # smoothness - penalize large changes between adjacent positions (softmax probs),
    # skipping any adjacency touching a padding position.
    pred_curr <- y_pred[, 2:seq_len,]       # positions 2:36
    pred_prev <- y_pred[, 1:(seq_len-1),]     # positions 1:35
    sm_mask   <- keep[, 2:seq_len] * keep[, 1:(seq_len-1)]              # (batch, seq-1)
    sq_pos    <- op_mean(op_square(pred_curr - pred_prev), axis = 3L)  # (batch, seq-1)
    smoothness <- op_sum(sq_pos * sm_mask, axis = 2L) /
                  (op_sum(sm_mask, axis = 2L) + eps)                    # (batch,)

    main + smoothness_weight * smoothness
  }


  masked_cat_accuracy <- custom_metric("masked_cat_accuracy", function(y_true, y_pred) {
    labels   <- y_true[, , 1:K_cat]                        # (batch, seq, K_cat) one-hot incl none
    pad_flag <- y_true[, , K_cat + 1L]                     # (batch, seq)
    keep     <- 1 - pad_flag                               # (batch, seq)
    correct  <- op_cast(op_equal(op_argmax(y_pred,  axis = -1L),
                                 op_argmax(labels, axis = -1L)), "float32") * keep
    op_sum(correct) / (op_sum(keep) + 1e-7)                # single-label accuracy, non-pad only
  })

  class_counts <- nn_in$train$y_global %>% table

  class_weights <- list(
    "0" = class_counts[2] / class_counts[1],  # Weight for majority class
    "1" = class_counts[2] / class_counts[2]   # Weight for minority class
  )

  # The oversampler (make_oversampled_dataset, n_neg = 3) already rebalances each
  # batch to ~1:3, so ALSO upweighting positives here double-corrects the
  # imbalance and overfits the knowns. Keep the global class weights ~1:1.
  # was: list("0" = 1, "1" = class_counts["0"] / class_counts["1"])  # ~20x on positives
  class_weights <- list(
    "0" = 1,
    "1" = 1
  )




  model |> compile(
    optimizer = optimizer_adam(
      # was: learning_rate_schedule_cosine_decay(4e-4, decay_steps = 10000).
      # With only ~600-1200 total training steps here vs decay_steps = 10000, the
      # cosine barely moved (LR stayed ~0.98x initial), so it was effectively a
      # constant LR anyway. Use a plain, slightly lower constant for the slim model.
      learning_rate = 3e-4,
      clipnorm = 1.0),
    loss = list(
      global = weighted_binary_crossentropy(weight_1 = class_weights[["1"]], weight_0 = class_weights[["0"]])
    ),
    metrics = list(
      global = list(
        metric_auc(name = "auc"),
        metric_auc(name = "pr_auc", curve = "PR")
      )
    )
  )


  resample_callback <- callback_lambda(
    on_epoch_begin = function(epoch, logs) {
      # rebuild dataset with fresh negative sample each epoch
      train_ds <<- make_oversampled_dataset(nn_in, n_negatives_per_positive = 3, noise_frac = aug_noise_frac)
    }
  )

  train_ds <- make_oversampled_dataset(nn_in, n_negatives_per_positive = 3, noise_frac = aug_noise_frac)

  n_pos <- sum(nn_in[["train"]][["y_global"]])
  steps_per_epoch <- ceiling(n_pos * (1 + 3) / 32)

  model |> fit(
    train_ds,
    validation_data = list(
      nn_in[["val"]][["x"]],
      list(global = nn_in[["val"]][["y_global"]])
    ),
    steps_per_epoch = steps_per_epoch,
    epochs = 1000,
    callbacks = list(early_stopping, resample_callback),
    batch_size = NULL
  )

  models[[term]] <- model

}








# Per-model Keras inputs, built once and reused for scoring (all), calibration
# (val) and the retrieval embeddings below. `all` is predicted ONCE per model.
# All concatenation is in names(nn_input) order so it stays aligned with
# bind_rows() of the $all sets.
nn_in_all   <- lapply(nn_input, generate_keras_input)

# predict the all-set ONCE per model -> global score + the unsupervised per_index_cat (for plots)
all_preds <- setNames(
  lapply(names(nn_input), function(term)
    predict(models[[term]], nn_in_all[[term]][["all"]][["x"]], verbose = 0)),
  names(nn_input))

# robust global-score extractor (predict returns list(global, per_index_cat))
.global_pred <- function(m, x) {
  p <- predict(m, x, verbose = 0)
  as.numeric((if (is.list(p)) p[["global"]] else p)[, 1])
}

# Raw (uncalibrated) global scores, kept alongside the calibrated ones.
raw_pred_comb <- do.call(c, lapply(all_preds, function(p) as.numeric(p[["global"]][, 1]))) %>% unname

# --- single POOLED Platt calibrator --------------------------------------
# Per-model calibration is unstable here (~few val positives -> near-perfect
# separation), and the models feed one shared global ranking, so fit ONE
# logistic calibrator on the val scores + labels POOLED across all models.
# This maps every model's raw score onto a common probability scale. Because the
# mapping is monotone it leaves each model's ROC/PR-AUC unchanged; it only makes
# the score scales comparable for the pooled ranking.
val_scores_pooled <- do.call(c, lapply(names(nn_input), function(term)
  .global_pred(models[[term]], nn_in_all[[term]][["val"]][["x"]])))
val_labels_pooled <- do.call(c, lapply(nn_in_all, function(nn) nn[["val"]][["y_global"]][, 1]))

cal_model <- glm(label ~ score,
                 data   = data.frame(score = val_scores_pooled, label = val_labels_pooled),
                 family = binomial)

cal_pred_comb <- unname(predict(cal_model,
                                newdata = data.frame(score = raw_pred_comb),
                                type = "response"))

# Rank on the calibrated (comparable) scores; raw kept for comparison (below).
val_pred_comb <- cal_pred_comb

metrics <- metric_set(
  roc_auc,
  pr_auc,
  accuracy,
  mcc,
  f_meas,
  precision,
  recall,
)

df <- tibble(
  truth = factor(do.call(c, lapply(nn_input, function(x) x$all$known))),
  .pred_class = ifelse(val_pred_comb > 0.6, 1, 0) %>% factor(., levels = c(0,1)),
  .pred_1 = val_pred_comb
)


qc_mets <- metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second")
qc_mets

extract_per_index <- function(arr, class_names) {
  d     <- dim(arr)
  n     <- d[1L]; seq_n <- d[2L]; K <- d[3L]
  idx   <- seq_len(seq_n)
  if (length(class_names) != K) class_names <- class_names[seq_len(K)]  # array has K real-class cols
  col_names <- c("index", class_names)

  lapply(seq_len(n), function(i) {
    m    <- arr[i, , ]                                          # (seq, K)
    cols <- c(list(idx), lapply(seq_len(K), function(k) m[, k]))
    names(cols) <- col_names
    tibble::new_tibble(cols, nrow = seq_n)                      # skip as_tibble/dplyr validation
  })
}

nn_input_comb <- bind_rows(lapply(nn_input, function(x) x$all))

nn_input_comb$pred     <- val_pred_comb   # calibrated (pooled Platt, common scale)
nn_input_comb$pred_raw <- raw_pred_comb   # uncalibrated global score, for comparison

# per_index_cat is an UNSUPERVISED output (no loss) kept only for visualization.
# NOTE: it is NOT trained against known_idx, so per-index *accuracy* is meaningless
# now -- it reflects whatever the global-ranking objective shaped, not the labels.
nn_input_comb$per_index <- do.call(c, lapply(all_preds, function(p) {
  extract_per_index(p[["per_index_cat"]], names(classes)[cat_cols])
}))

# --- nearest known-peptide retrieval (reuse the trained "embed" layer) ------
# Tap the 32-d penultimate representation, L2-normalise, and for each window
# find the most similar KNOWN peptide by cosine similarity. No new loss: the
# embedding is whatever the global classifier already learned to sit on.
l2norm <- function(m) m / sqrt(pmax(rowSums(m^2), 1e-12))

emb_all <- do.call(rbind, lapply(names(nn_input), function(term) {
  embed_model <- keras_model(models[[term]]$input, get_layer(models[[term]], "embed")$output)
  l2norm(predict(embed_model, nn_in_all[[term]][["all"]][["x"]]))
}))                                                          # rows aligned with nn_input_comb

# reference bank = known peptides only
ref_i     <- which(nn_input_comb$known == 1)
ref_emb   <- emb_all[ref_i, , drop = FALSE]
ref_names <- nn_input_comb$peps[ref_i]

sim_mat <- emb_all %*% t(ref_emb)                            # (n_windows, n_ref) cosine sims
sim_mat[outer(nn_input_comb$peps, ref_names, `==`)] <- -Inf # never match a window to itself
best <- max.col(sim_mat, ties.method = "first")

nn_input_comb$nn_closest_peptide <- ref_names[best]
nn_input_comb$nn_closest_sim     <- sim_mat[cbind(seq_along(best), best)]

nn_input_comb <- nn_input_comb %>%
  # per-model (terminus) category, and its combination with win_type
  mutate(model    = stringr::str_remove(as.character(target), "^loop_"),  # "N" or "C"
         category = paste0(model, "_", win_type)) %>%                     # e.g. "N_db", "C_chym", "C_pep_end"
  arrange(desc(pred)) %>%
  #filter(win_type == "db") %>%
  mutate(rank = row_number(), .by = known) %>%                            # existing global rank (within known/unknown)
  # rank the score WITHIN each combined (model x win_type) category, kept
  # separate for known vs candidate. Drop `known` from .by to rank the two together.
  mutate(rank_cat = row_number(), .by = c(category, known))

uniprot_peps <- data.table::fread("~/Desktop/Peptides/uniprot_peptides.csv") %>% as_tibble()

nn_input_comb <- nn_input_comb %>%
  mutate(gene = stringr::str_extract(peps, "^[^_]+")) %>%
  {left_join(., secretome %>% distinct(gene, .keep_all = T) %>% select(gene, location), by = "gene")} %>%
  mutate(uni_pep = if_else(gene %in% uniprot_peps$gene, 1, 0))

## peptide length from the pep_id "<start>x<end>" range. Vectorized (one str_match
## over the whole column) instead of a per-row map_int -- the latter ran the regex
## ~40K times and returned NA for every non-known window. NA where no range (controls).
nn_input_comb <- nn_input_comb %>%
  mutate(length = {
    rng <- stringr::str_match(pep_id, "(\\d+)x(\\d+)")
    as.integer(rng[, 3]) - as.integer(rng[, 2])
  })

## --- annotate windows that ANCHOR a uniprot-peptide terminus -----------------
## Windows are built as anchor + [-5,30] (N) / [-30,5] (C), so the putative peptide
## boundary sits at a fixed anchor residue: wN+5 for N windows, wC-5 for C windows.
## A window "hits" a uniprot peptide only if that peptide's matching terminus --
## start (N-terminus) for an N window, end (C-terminus) for a C window -- lands at
## the window's anchor residue (+/- anchor_tol), i.e. at the correct position, not
## merely somewhere inside the window span.
anchor_tol <- 2L

## per-model terminus (N/C); derive here so this block is self-contained even if
## the ranking mutate above hasn't been run on this nn_input_comb
if (!"model" %in% names(nn_input_comb))
  nn_input_comb$model <- stringr::str_remove(as.character(nn_input_comb$target), "^loop_")

wm <- stringr::str_match(nn_input_comb$peps, "_w(\\d+)-(\\d+)$")
nn_input_comb$wN         <- as.integer(wm[, 2])
nn_input_comb$wC         <- as.integer(wm[, 3])
nn_input_comb$anchor_res <- ifelse(nn_input_comb$model == "N",
                                   nn_input_comb$wN + 5L,     # expected peptide N-terminus
                                   nn_input_comb$wC - 5L)     # expected peptide C-terminus

starts_by_gene <- split(as.integer(uniprot_peps$start), uniprot_peps$gene)  # peptide N-ends
ends_by_gene   <- split(as.integer(uniprot_peps$end),   uniprot_peps$gene)  # peptide C-ends

anchor_hits <- function(gene, model, anchor, tol) {
  ter <- if (model == "N") starts_by_gene[[gene]] else ends_by_gene[[gene]]
  if (is.null(ter) || is.na(anchor)) return(FALSE)
  any(abs(ter - anchor) <= tol)
}

nn_input_comb$pep_terminus_hit <- FALSE
idx <- which(nn_input_comb$gene %in% names(starts_by_gene))   # only genes with uniprot peptides
if (length(idx) > 0) {
  nn_input_comb$pep_terminus_hit[idx] <- mapply(
    anchor_hits,
    nn_input_comb$gene[idx], nn_input_comb$model[idx], nn_input_comb$anchor_res[idx],
    MoreArgs = list(tol = anchor_tol))
}

message(sprintf("uniprot-terminus hits: %d windows (N: %d, C: %d)",
                sum(nn_input_comb$pep_terminus_hit),
                sum(nn_input_comb$pep_terminus_hit & nn_input_comb$model == "N"),
                sum(nn_input_comb$pep_terminus_hit & nn_input_comb$model == "C")))

## --- violin: score of hits vs non-hits, segregated by window terminus (N/C) ---
## swap `y = pred` for `y = rank_cat` (or `rank`) to plot rank instead of score.
pep_hit_violin <- nn_input_comb %>%
  mutate(hit = factor(if_else(pep_terminus_hit, "anchors uniprot pep terminus", "no"),
                      levels = c("no", "anchors uniprot pep terminus"))) %>%
  ggplot(aes(x = hit, y = pred, fill = hit)) +
  geom_violin(scale = "width", alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(data = ~ dplyr::filter(.x, pep_terminus_hit),
              width = 0.15, size = 0.9, alpha = 0.7) +
  facet_grid(rows = vars(win_type), cols = vars(model)) +
  labs(x = NULL, y = "score (pred)",
       title = "Windows anchoring a uniprot-peptide terminus vs not, by model x win_type") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("~/AF2_analysis/uniprot_terminus_score_violin.svg", pep_hit_violin, width = 9, height = 5)

nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(peps %in% c("ANO8_w12-47", "ANO8_w36-71", "ANO8_w14-49")) %>%
  View()

nn_input_comb %>% filter(known == 0) %>%
  filter(gene == "ANO8") %>%
  View()


nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(grepl("GDF", gene)) %>%
  View()

nn_input_comb %>% filter(known == 0) %>%
  filter(grepl("ASIP", gene)) %>%
  View()


View(nn_input_comb %>%
       filter(known == 0 & location %in% c("2t", "3t", "2l", "3l", "4l")) %>%
       filter(win_type == "chym") %>%
       filter(target %in% c("C", "C_loop")))

View(nn_input_comb)

nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(grepl("BRINP3", gene) & target == "C") %>%
  View()






secretome <- dplyr::left_join(
  secretome,
  secretome_aa %>% dplyr::select(!matches("_lead|_lag")) %>%
    dplyr::group_by(accession) %>% tidyr::nest(.key = "aa_scores"),
  by = "accession"
)

secretome[["aa_scores"]] <- purrr::map2(
  secretome[["aa_scores"]], secretome[["sequence_uni"]],
  \(x, y) {
    if (is.null(x)) return(x)
    tmp <- tibble(index = 1:nchar(y), AA = stringr::str_split(y, "", simplify = TRUE) %>% c)
    dplyr::left_join(tmp, x, by = c("index", "AA"))
  }
)


data.table::fwrite(nn_input_comb %>% select(!where(is.list)), "~/Desktop/scores3.csv")




genes <- c("CPXM1", "MDK", "FGF5", "ARSI", "LACRT", "PDGFB", "PCSK1N",
           "KRTDAP", "NDFIP1", "RTBDN", "HMCN2", "IL21", "GNPTG", "ANGPTL8",
           "MANF", "GRP", "CCDC3", "DMKN", "CASP4", "SERPINA9", "CXCL3",
           "CXCL17", "DHRS4L2", "IBSP", "PROS1", "GAS6", "CCDC126", "GNPTG",
           "NMB", "SEMG2", "LPO", "PCSK1N", "FGF6", "BRINP3", "ENAM", "INSL5",
           "SOSTDC1")


plot_dir <- "~/AF2_analysis/new_meth_plot2"


#secretome <- readRDS("~/AF2_analysis/secretome_latest.rds")

#peps_tp <- readRDS("~/AF2_analysis/peps_tp.rds")


mets <- list(cons = c("blos_wt_all_n", "cons_rs", "blos_wt_mam", "blos_wt_all", "gran_wt_all"),
             af_missense = c("mean_afm", "min_afm"),
             dssp = c("relASA"),
             aa_scores = c("pep_xgb4c", "chem_xgb3c", "pep_nn4c", "chem_nn4c")
)

mets$aa_scores <- c(mets$aa_scores, paste0(mets$aa_scores, "_s6"))

all_mets <- do.call(`c`, mets) %>% unname

names(all_mets) <- rep(names(mets), sapply(mets, length))





dir.create(plot_dir)
unlink(plot_dir)

the_input <- secretome %>%
  filter(!accession %in% c("A0AAG2TCD0", "A0AAG2UXZ5")) %>%
  filter(gene %in% !!genes) %>%
  mutate(aa_scores = map(aa_scores, \(x) x[, colnames(x) %in% mets[["aa_scores"]]])) %>%
  group_split(gene)

pep_input <- peps_tp %>%
  filter(gene %in% genes) %>%
  group_split(gene)

input_gene <- map_chr(the_input, \(x) x$gene)

pep_input_gene <- map_chr(pep_input, \(x) x$gene)

identical(input_gene, pep_input_gene)

pred_to_plot <- nn_input_comb


species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))


devtools::load_all("/Users/kbrulois/R_projects/ligandFinder")

# baseline reference
#make_protein_plot(the_input[[1]], pep_input[[1]])
# new variant

for(i in seq_along(the_input)) {
make_protein_plot_win(the_input[[i]], pred_to_plot, plot_dir, pep_input[[i]])
}





devtools::load_all("/Users/kbrulois/R_projects/ligandFinder")
make_protein_plot_win(the_input[[2]], nn_input_comb, plot_dir, pep_input[[2]])

# verify nnwin_ data-id is now in SVG
html <- readLines(file.path(plot_dir, paste0(the_input[[1]]$gene, ".html")), warn = FALSE)
sum(stringr::str_count(html, "data-id='nnwin_"))
sum(stringr::str_count(html, "onclick='nn_show_panel"))
sum(stringr::str_count(html, "nn_show_panel"))
sum(stringr::str_count(html, "nnwin_"))
sum(stringr::str_count(html, '<svg'))
sum(stringr::str_count(html, "nn-detail-panel.active"))
sum(stringr::str_count(html, "#nn-detail-container"))
sum(stringr::str_count(html, "nn-status-badge"))
sum(stringr::str_count(html, "nn_status"))
file.info(file.path(plot_dir, paste0(the_input[[1]]$gene, ".html")))$mtime




nn_input$val$pred <- val_pred[,1]

p <- yardstick::roc_curve(data = df, truth = "truth", ".pred_1",
                          event_level = "second") %>%
  autoplot() +
  ggtitle("Mixed training set; Mixed test",
          subtitle = paste(paste0("roc auc: ", round(qc_mets %>% filter(.metric == "roc_auc") %>% pull(.estimate), 2)),
                           paste0("pr auc: ", round(qc_mets %>% filter(.metric == "pr_auc") %>% pull(.estimate), 2)), collapse = "\n"))


ggsave("~/AF2_analysis/all_peps_roc_AUC_mixed-train-mixed_test.svg", p)

model_stats <- metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second")

write.csv(model_stats, "~/AF2_analysis/model_stats_new.csv")




predict_combined <- function(x, category) {
  ifelse(
    category == 1,
    predict(model_A, x),
    predict(model_B, x)
  )
}














