
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


classes <- c(nn_input$N$train$known_idx, nn_input$C$train$known_idx) %>% do.call(`c`, .) %>% unique
classes <- names(class_cols)
classes <- setNames(0:(length(classes) - 1), classes)


seq_len <- nn_input$N$train$data %>% `[[`(1) %>% nrow
n_channels <- nn_input$N$train$data %>% `[[`(1) %>% ncol
K <- c(nn_input$N$train$known_idx, nn_input$C$train$known_idx) %>% do.call(`c`, .) %>% unique %>% length
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

library(tfdatasets)

make_oversampled_dataset <- function(nn_in, n_negatives_per_positive = 3, batch_size = 32) {

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

  x_t <- tensorflow::as_tensor(x_train[idx,,], dtype = "float32")
  yg_t <- tensorflow::as_tensor(y_global[idx,, drop = FALSE], dtype = "float32")
  yc_t <- tensorflow::as_tensor(y_cat[idx,,], dtype = "float32")

  tensor_slices_dataset(list(x_t, yg_t, yc_t)) |>
    dataset_shuffle(buffer_size = length(idx)) |>
    dataset_batch(batch_size) |>
    dataset_map(function(x, yg, yc) {
      list(x, list(global = yg, per_index_cat = yc))
    }) |>
    dataset_repeat()  # repeat so keras can run multiple epochs
}
models <- list()

val_pred <- list()

for(term in names(nn_input)[1:2]) {

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

      y_per_index_cat <- array(input[["known_idx"]] %>%
                                 do.call(`c`, .) %>%
                                 classes[.] %>%
                                 to_categorical(.),
                               dim = c(seq_len, n_samples,  K)) %>%
        aperm(., c(2,1,3))


      list(x = x,
           y_global = y_global,
           y_per_index_cat = y_per_index_cat)

    }, nn_in)

  }

  nn_in <- generate_keras_input(nn_in = nn_input[[term]])


  inputs <- layer_input(shape = c(seq_len, n_channels))

  regularizer <- regularizer_l2(1e-3)

  shared <- inputs |>
    layer_conv_1d(32, kernel_size = 16, activation = "gelu",
                  kernel_regularizer = regularizer, padding = "same") |>
    layer_conv_1d(32, kernel_size = 8, activation = "gelu",
                  kernel_regularizer = regularizer, padding = "same") |>
    layer_conv_1d(16, kernel_size = 4, padding = "same",
                  activation = "gelu", dilation_rate = 1, kernel_regularizer = regularizer) |>
    #layer_batch_normalization() |>
    layer_dropout(0.2) |>
    layer_conv_1d(16, kernel_size = 2, padding = "same",
                  activation = "gelu", dilation_rate = 1, kernel_regularizer = regularizer) |>
    #layer_batch_normalization() |>
    layer_dropout(0.4)

  per_index_logits <- shared |>
    layer_conv_1d(filters = K, kernel_size = 1,
                  kernel_regularizer = regularizer)

  per_index_cat <- per_index_logits |>
    layer_activation("softmax", name = "per_index_cat")


  masked_sum <- per_index_logits |>
    layer_lambda(function(x) {
      mask <- keras$ops$convert_to_tensor(mask_matrix, dtype = "float32")
      keras$ops$multiply(x, mask)
    }, output_shape = c(seq_len, K)) |>
    layer_activation("softmax")

  if(FALSE) {
    global_output <- masked_sum |>
      layer_global_average_pooling_1d() |>
      layer_dense(16, activation = "gelu") |>
      layer_dense(1, activation = "sigmoid", name = "global")
  }

  global_output <- layer_multi_head_attention(num_heads = 4, key_dim = 32)(masked_sum, masked_sum) |>
    layer_dropout(0.4) |>
    layer_layer_normalization() |>
    layer_global_average_pooling_1d() |>
    layer_dense(8, activation = "gelu") |>
    layer_dropout(0.3) |>
    layer_dense(1, activation = "sigmoid", name = "global")



  lr_reduce <- callback_reduce_lr_on_plateau(monitor = "val_loss", factor = 0.5, patience = 20, min_lr = 1e-6)

  early_stopping <- callback_early_stopping(
    monitor = "val_global_auc",
    patience = 50,
    mode = 'max',
    start_from_epoch = 30,
    restore_best_weights = TRUE
  )


  model <- keras_model(
    inputs = inputs,
    outputs = list(global = global_output,
                   per_index_cat = per_index_cat)
  )



  masked_focal_loss <- function(y_true, y_pred, gamma = 2.0, smoothness_weight = 0.1) {

    # focal loss component
    ce <- loss_categorical_crossentropy(y_true, y_pred, from_logits = FALSE, label_smoothing = 0.05)
    p_t <- op_sum(y_true * y_pred, axis = -1L)
    focal_factor <- op_power(1.0 - p_t, gamma)
    class_ids <- op_argmax(y_true, axis = -1L)
    mask <- op_cast(op_not_equal(class_ids, 8L), "float32")
    loss <- ce * focal_factor * mask
    focal <- op_sum(loss, axis = 2L) / (op_sum(mask, axis = 2L) + 1e-7)  # (batch,)

    # smoothness component - penalize large changes between adjacent positions
    # y_pred shape: (batch, seq_len, 8)
    pred_curr <- y_pred[, 2:seq_len,]    # positions 2:36
    pred_prev <- y_pred[, 1:(seq_len-1),]  # positions 1:35
    smoothness <- op_mean(op_square(pred_curr - pred_prev), axis = c(2L, 3L))  # (batch,)

    focal + smoothness_weight * smoothness
  }


  masked_cat_accuracy <- custom_metric("masked_cat_accuracy", function(y_true, y_pred) {
    class_ids <- op_argmax(y_true, axis = -1L)
    mask <- op_cast(op_not_equal(class_ids, 8L), "float32")
    correct <- op_cast(
      op_equal(op_argmax(y_true, axis = -1L), op_argmax(y_pred, axis = -1L)),
      "float32"
    )
    op_sum(correct * mask) / (op_sum(mask) + 1e-7)
  })

  class_counts <- nn_in$train$y_global %>% table

  class_weights <- list(
    "0" = class_counts[2] / class_counts[1],  # Weight for majority class
    "1" = class_counts[2] / class_counts[2]   # Weight for minority class
  )

  class_weights <- list(
    "0" = 1,
    "1" = class_counts["0"] / class_counts["1"]  # upweight minority
  )




  model |> compile(
    optimizer = optimizer_adam(
      learning_rate = learning_rate_schedule_cosine_decay(
        initial_learning_rate = 4e-4,
        decay_steps = 10000
      ),
      clipnorm = 1.0),
    loss = list(
      global = weighted_binary_crossentropy(weight_1 = class_weights[["1"]], weight_0 = class_weights[["0"]]),
      per_index_cat = masked_focal_loss
    ),
    loss_weights = list(
      global = 0.2,
      per_index_cat = 1),
    metrics = list(
      global = list(
        metric_auc(name = "auc"),
        metric_auc(name = "pr_auc", curve = "PR")
      ),
      per_index_cat = list(
        masked_cat_accuracy,
        "categorical_accuracy"
      )
    )
  )


  resample_callback <- callback_lambda(
    on_epoch_begin = function(epoch, logs) {
      # rebuild dataset with fresh negative sample each epoch
      train_ds <<- make_oversampled_dataset(nn_in, n_negatives_per_positive = 3)
    }
  )

  train_ds <- make_oversampled_dataset(nn_in, n_negatives_per_positive = 3)

  n_pos <- sum(nn_in[["train"]][["y_global"]])
  steps_per_epoch <- ceiling(n_pos * (1 + 3) / 32)

  model |> fit(
    train_ds,
    validation_data = list(
      nn_in[["val"]][["x"]],
      list(
        global = nn_in[["val"]][["y_global"]],
        per_index_cat = nn_in[["val"]][["y_per_index_cat"]]
      )
    ),
    steps_per_epoch = steps_per_epoch,
    epochs = 300,
    callbacks = list(early_stopping, resample_callback),
    batch_size = NULL
  )

  models[[term]] <- model

  val_pred[[term]] <- models[[term]] %>%
    predict(nn_in[["all"]][["x"]]) %>%
    `[[`("global")

}








calibrate_model <- function(model, nn_in) {
  # Get raw predictions on validation set
  val_preds <- model |> predict(nn_in[["val"]][["x"]])
  val_scores <- val_preds[["global"]][, 1]  # positive class score
  val_labels <- nn_in[["val"]][["y_global"]][, 1]  # positive class label

  # Fit logistic regression to map raw scores -> calibrated probabilities
  cal_data <- data.frame(score = val_scores, label = val_labels)
  cal_model <- glm(label ~ score, data = cal_data, family = binomial)
  cal_model
}

apply_calibration <- function(raw_scores, cal_model) {
  predict(cal_model, newdata = data.frame(score = raw_scores), type = "response")
}

nn_in_N <- generate_keras_input(nn_input[["N"]])
nn_in_C <- generate_keras_input(nn_input[["C"]])

cal_N <- calibrate_model(models$N, nn_in_N)
cal_C <- calibrate_model(models$C, nn_in_C)

# Apply to get comparable scores
scores_N <- models$N |> predict(nn_in_N[["all"]][["x"]])
scores_C <- models$C |> predict(nn_in_C[["all"]][["x"]])

calibrated_N <- apply_calibration(scores_N[["global"]][, 1], cal_N)
calibrated_C <- apply_calibration(scores_C[["global"]][, 1], cal_C)




#val_pred_comb <- map(val_pred, ~.[,1]) %>% do.call(`c`, .) %>% unname
val_pred_comb <- c(calibrated_N, calibrated_C)

val_pred_comb <- c(scores_N[["global"]][, 1], scores_C[["global"]][, 1])

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
  truth = factor(c(nn_input$N$all$known,
                   nn_input$C$all$known)),
  .pred_class = ifelse(val_pred_comb > 0.6, 1, 0) %>% factor(., levels = c(0,1)),
  .pred_1 = val_pred_comb
)


qc_mets <- metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second")
qc_mets

nn_input_comb <- bind_rows(nn_input$N$all, nn_input$C$all)

nn_input_comb$pred <- val_pred_comb

nn_input_comb <- nn_input_comb %>%
  arrange(desc(pred)) %>%
  #filter(win_type == "db") %>%
  mutate(rank = row_number(), .by = known)

uniprot_peps <- data.table::fread("~/Desktop/Peptides/uniprot_peptides.csv") %>% as_tibble()

nn_input_comb <- nn_input_comb %>%
  mutate(gene = stringr::str_extract(peps, "^[^_]+")) %>%
  {left_join(., secretome %>% distinct(gene, .keep_all = T) %>% select(gene, location), by = "gene")} %>%
  mutate(uni_pep = if_else(gene %in% uniprot_peps$gene, 1, 0))

nn_input_comb <- nn_input_comb %>%
  mutate(length = map_int(pep_id, \(x) {

    inds <- stringr::str_extract(x, "\\d+x\\d+") %>% stringr::str_split(., "x", simplify = TRUE) %>% `c` %>% as.numeric
    if(!is.na(inds[1])) {
      inds[2] - inds[1]
    } else {
      NA
    }}))

nn_input_comb %>% filter(known == 0 & win_type == "db") %>%
  filter(peps %in% c("ANO8_w12-47", "ANO8_w36-71", "ANO8_w14-49") | gene %in% c("ANO8", "CXCL14", "GDF5", "GDF6") | grepl("BRINP", gene)) %>%
  View()

View(nn_input_comb %>%
       filter(known == 0 & location %in% c("2t", "3t", "2l", "3l", "4l")) %>%
       filter(win_type == "db") %>%
       filter(target %in% c("C", "C_loop")))

View(nn_input_comb)





test <- nn_input_comb %>%
            filter(known == 0 & location %in% c("2t", "3t", "2l", "3l", "4l")) %>%
            filter(win_type == "db") %>%
            filter(target %in% c("C", "C_loop")) %>%
            slice(1:50)

test2 <- nn_input_comb %>%
  filter(known == 0 & location %in% c("2t", "3t", "2l", "3l", "4l")) %>%
  filter(win_type == "chym") %>%
  filter(target %in% c("C", "C_loop")) %>%
  slice(1:50)


saveRDS(nn_input_comb, "~/AF2_analysis/April17_peps.rds")

saveRDS(test, "~/AF2_analysis/db_peps.rds")

saveRDS(test2, "~/AF2_analysis/chym_peps.rds")

test <- readRDS("~/AF2_analysis/chym_peps.rds")

id_map <- readRDS(system.file("/data/id_mapping.rds", package = "ligandFinder"))

test2 <- ligandFinder:::window_to_mature_peptide(meta_data = test$meta_data, target = rep("C", 50), gene = test$gene)

gene_to_uni <- setNames(id_map$`Entry Name`, id_map$`Gene Names (primary)`)




test2 <- test2 %>%
  mutate(
    uniprot_name = gene_to_uni[gene],
    start = ifelse(is.na(mature_peptide),
                   purrr::map_int(post_inserting_index, ~ if (length(.x)) min(.x) else NA_integer_),
                   purrr::map_int(mature_index,         ~ if (length(.x)) min(.x) else NA_integer_)),
    end   = ifelse(is.na(mature_peptide),
                   purrr::map_int(post_inserting_index, ~ if (length(.x)) max(.x) else NA_integer_),
                   purrr::map_int(mature_index,         ~ if (length(.x)) max(.x) else NA_integer_))
  )

test3 <- test2 %>% select(uniprot_name, start, end) %>% filter(!is.na(uniprot_name))



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







