
vocab_size <- 20000   # total number of tokens
d_model <- 64
num_heads <- 4
ff_dim <- 128

inputs <- layer_input(shape = c(seq_len))

x <- inputs |>
  layer_embedding(
    input_dim = vocab_size,
    output_dim = d_model
  )

pos_encoding <- positional_encoding(seq_len, d_model)

pos_encoding <- tf$expand_dims(pos_encoding, axis = 0)

x <- x + pos_encoding
# --------------------------------------------------
# 3️⃣ Transformer Block
# --------------------------------------------------

# Multi-head self attention
attn_output <- layer_multi_head_attention(
  num_heads = num_heads,
  key_dim = d_model
)(x, x)

attn_output <- attn_output |>
  layer_dropout(rate = 0.3)
# Residual + LayerNorm
x1 <- layer_layer_normalization()(x + attn_output)

# Feed-forward network
ffn_output <- x1 |>
  layer_dense(ff_dim, activation = "relu",
              kernel_regularizer = regularizer) |>
  layer_dense(d_model)

ffn_output <- ffn_output |>
  layer_dropout(rate = 0.3)
# Second residual + norm
x2 <- layer_layer_normalization()(x1 + ffn_output)

# --------------------------------------------------
# 4️⃣ Per-index Head
# --------------------------------------------------

per_index_logits <- x2 |>
  layer_dense(K)

per_index_output <- per_index_logits |>
  layer_activation("softmax", name = "per_index")

# --------------------------------------------------
# 5️⃣ Global Head
# --------------------------------------------------

global_output <- per_index_logits |>
  layer_global_max_pooling_1d() |>
  layer_dense(32, activation = "relu") |>
  layer_dense(1, activation = "sigmoid", name = "global")

# --------------------------------------------------
# 6️⃣ Build Model
# --------------------------------------------------

model <- keras_model(
  inputs = inputs,
  outputs = list(
    global = global_output,
    per_index = per_index_output
  )
)












library(keras3)
library(tensorflow)
library(yardstick)



classes <- c(nn_input$N$train$known_idx, nn_input$C$train$known_idx) %>% do.call(`c`, .) %>% unique
classes <- classes[classes != "unknown"]
classes <- setNames(0:(length(classes) - 1), classes)
seq_len <- nn_input$N$train$data %>% `[[`(1) %>% nrow
n_channels <- nn_input$N$train$data %>% `[[`(1) %>% ncol
K <- nn_input$N$train$known_idx %>% do.call(`c`, .) %>% unique %>% length
kernel_size <- 8


models <- list()

val_pred <- list()

for(term in c("N", "C")) {

  rm(model,per_index_logits, shared,  per_index_logits, per_index_output, per_index_output2, global_output, early_stopping)
  keras3::clear_session()

  nn_in <- nn_input[[term]]

  n_samples <- length(nn_in$train$data)

  x_train <- array(
    unlist(nn_in$train$data),
    dim = c(seq_len, n_channels, length(nn_in$train$data))
  )

  x_train <- aperm(x_train, c(3, 1, 2))


  y_train_global <- matrix(nn_in$train$known,
                           ncol = 1)


  y_train_per_index <- array(nn_in$train$known_idx %>% do.call(`c`, .) %>% classes[.] %>% to_categorical(.),
                             dim = c(seq_len, n_samples,  K)) %>% aperm(., c(2,1,3))
  #y_train_per_index[,,8] <- 0



  n_samples <- length(nn_in$val$data)

  x_val <- array(
    unlist(nn_in$val$data),
    dim = c(seq_len, n_channels, length(nn_in$val$data))
  )

  x_val <- aperm(x_val, c(3, 1, 2))


  y_val_global <- matrix(nn_in$val$known,
                         ncol = 1)

  y_val_per_index <- array(nn_in$val$known_idx %>% do.call(`c`, .) %>% classes[.] %>% to_categorical(.),
                           dim = c(seq_len, n_samples,  K)) %>% aperm(., c(2,1,3))


  inputs <- layer_input(shape = c(seq_len, n_channels))

  regularizer <- regularizer_l2(1e-4)

  shared <- inputs |>
    layer_conv_1d(32, kernel_size = 3, padding = "same",
                  activation = "gelu", dilation_rate = 1, kernel_regularizer = regularizer) |>
    #layer_batch_normalization() |>
    layer_conv_1d(64, kernel_size = 6, padding = "same",
                  activation = "gelu", dilation_rate = 2, kernel_regularizer = regularizer)
    #layer_batch_normalization() |>
    #layer_dropout(0.3)

  # Per-index head
  per_index_logits <- shared |>
    layer_conv_1d(filters = K, kernel_size = 1, activation = NULL)

  per_index_output <- per_index_logits |>
    layer_activation("softmax", name = "per_index")

  per_index_output2 <- shared |>
    layer_conv_1d(filters = K - 1, kernel_size = 1, activation = NULL) |>
    layer_activation("softmax", name = "per_index2")


  # Global head (separate branch)
  global_output <- shared |>
    layer_conv_1d(32, kernel_size = 8, activation = "gelu") |>
    layer_global_max_pooling_1d() |>
    layer_dense(1, activation = "sigmoid", name = "global")


  early_stopping <- callback_early_stopping(
    monitor = "per_index_accuracy",
    patience = 40,
    mode = 'max',
    start_from_epoch = 20,
    restore_best_weights = TRUE
  )


  model <- keras_model(
    inputs = inputs,
    outputs = list(global = global_output,
                   per_index = per_index_output,
                   per_index2 = per_index_output2,
                   per_index_logits = per_index_logits)
  )

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
    mask <- op_cast(op_not_equal(class_ids, 7L), "float32")

    # 2. INCREASE THE SIGNAL: multiply the loss of non-7 classes by 5 or 10
    # This forces the model to treat mistakes on real classes as "expensive"
    boosted_loss <- ce * mask * 100.0

    return(op_sum(boosted_loss, axis = -1L))
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
    mask <- op_cast(op_not_equal(class_ids, 7L), "float32")

    # 4. Multiply and return
    return(op_sum(ce * focal_factor * mask, axis = -1L))
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



  model |> compile(
    optimizer = optimizer_adam(
      learning_rate = learning_rate_schedule_cosine_decay(
        initial_learning_rate = 1e-3,
        decay_steps = 10000
      ),
      clipnorm = 1.0),
    loss = list(
      global = "binary_crossentropy",
      per_index = masked_focal_loss
    ),
    loss_weights = list(
      global = 0.2,
      per_index = 1
    ),
    metrics = list(
      global = list(
        "accuracy",
        metric_auc(name = "auc")
      ),
      per_index = list(
        "accuracy"
      )
    )
  )


  model |> fit(
    x = x_train,
    y = list(
      global = y_train_global,
      per_index = y_train_per_index
    ),
    epochs = 300,
    callbacks = list(early_stopping),
    batch_size = 32
  )

  models[[term]] <- model

  val_pred[[term]] <- models[[term]] %>%
    predict(x_val) %>%
    `[[`("global")

}

models[["N"]] |> save_model("/Users/kbrulois/AF2_analysis/model_N_Feb24.keras")
models[["C"]] |> save_model("/Users/kbrulois/AF2_analysis/model_C_Feb24.keras")


val_pred_comb <- map(val_pred, ~.[,1]) %>% do.call(`c`, .) %>% unname


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
  truth = factor(c(nn_input$N$val$known,
                   nn_input$C$val$known)),
  .pred_class = ifelse(val_pred_comb > 0.6, 1, 0) %>% factor(., levels = c(0,1)),
  .pred_1 = val_pred_comb
)

#df <- bind_rows(df_c, df_n)

metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second")


nn_input_comb <- bind_rows(nn_input$N$val, nn_input$C$val)

nn_input_comb$pred <- val_pred_comb

nn_input_comb <- nn_input_comb %>%
  arrange(desc(pred)) %>%
  mutate(rank = row_number(), .by = known)

nn_input_comb <- nn_input_comb %>%
  mutate(gene = stringr::str_extract(pep_name, "^[^_]+")) %>%
  {left_join(., secretome %>% distinct(gene, .keep_all = T) %>% select(gene, location), by = "gene")} %>%
  mutate(uni_pep = if_else(gene %in% uniprot_peps$gene, 1, 0))


nn_input_comb %>% filter(known == 0) %>% filter(pep_name %in% c("ANO8_w12-47", "ANO8_w36-71"))

View(nn_input_comb %>% filter(known == 0))

uniprot_peps <- data.table::fread("~/Desktop/uniprot_peptides.csv") %>% as_tibble()


nn_input$val$pred <- val_pred[,1]

p <- yardstick::roc_curve(data = df, truth = "truth", ".pred_1",
                          event_level = "second") %>%
  autoplot() + ggtitle("All peptide prediction")


ggsave("~/AF2_analysis/all_peps_roc_AUC.svg", p)

model_stats <- bind_cols(
  metrics(df_n, truth = truth, estimate = .pred_class, .pred_1, event_level = "second") %>%
    dplyr::rename(N_term = .estimate) %>%
    select(.metric, N_term),
  metrics(df_c, truth = truth, estimate = .pred_class, .pred_1, event_level = "second") %>%
    dplyr::rename(C_term = .estimate) %>%
    select(C_term),
  metrics(df, truth = truth, estimate = .pred_class, .pred_1, event_level = "second") %>%
    dplyr::rename(all = .estimate) %>%
    select(all)
)

write.csv(model_stats, "~/AF2_analysis/model_stats_new.csv")




predict_combined <- function(x, category) {
  ifelse(
    category == 1,
    predict(model_A, x),
    predict(model_B, x)
  )
}










saliency_fn <- function(model, x_batch) {
  x_tensor <- tf$convert_to_tensor(x_batch)

  with(tf$GradientTape() %as% tape, {
    tape$watch(x_tensor)
    preds <- model(x_tensor)[[1]]  # global output
  })

  grads <- tape$gradient(preds, x_tensor)
  as.array(grads)
}


targs <- list(
  c("C", "loop_C"),
  c("N", "loop_N")
)


met_it <- list(core = all_params3[1:4],
               motif1 = c("AA_Y", "AA_W", "AA_F"),
               motif2 = c("AA_Q", "AA_E", "AA_D"),
               motif3 = c("AA_G", "AA_H", "AA_V"),
               SS = c("SS_G", "SS_B", "SS_H"))

for(y in seq_along(met_it)) {

  mets <- met_it[[y]]
  p_final <- list()


  for(x in seq_along(targs)) {

    targ <- targs[[x]]

    shap_vals <- saliency_fn(models[[targ[1]]], x_val)


    dat_toplot <- lapply(mets, \(x) shap_vals[y_val_global == 1,,which(all_params3 == x)] %>% as_tibble %>% mutate(metric = x)) %>%
      bind_rows() %>%
      pivot_longer(cols = -metric) %>%
      mutate(index = setNames(1:36, paste0("V", 1:36))[name]) %>%
      select(-name)

    non_zero_mean <- function(x) {
      mean(x[!x == 0], na.rm = TRUE)
    }

    center_func <- if(names(met_it)[y] == "core") {center_func <- non_zero_mean} else {center_func <- mean}

    dat_toplot2 <- map(known_dat$data[known_dat$target %in% targ], \(x) {

      x[, mets] %>%
        mutate(index = row_number())

    }) %>%
      bind_rows %>%
      group_by(index) %>%
      summarise(across(everything(), center_func),
                index = first(index)) %>%
      pivot_longer(cols = -index) %>%
      mutate(type = "average of residues within or adjacent to known GPCR ligands")

    dat_toplot3 <- map(c_dat$data[c_dat$target %in% targ], \(x) {

      x[, mets] %>%
        mutate(index = row_number())

    }) %>%
      bind_rows %>%
      pivot_longer(cols = -index) %>%
      mutate(type = "average of residues from ~40K candidate regions")



    p <- ggplot2::ggplot(dat_toplot %>% group_by(index, metric) %>%
                           summarise(across(everything(), mean, na.rm = T))) +
      #ggplot2::geom_violin(aes(x = name, y = value), trim = TRUE) +
      ggplot2::geom_point(aes(x = index, y = value)) +
      ylab("model sensitivity") +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5)) +
      ggplot2::facet_grid(rows = vars(metric), scales = "free_y") +
      theme_bw() +
      theme(
        strip.background = element_blank(),     # removes grey box
        strip.text = element_text(
          color = "black",
          size = 10
        )  )


    p2 <- ggplot2::ggplot(data = dat_toplot3, aes(x = index, y = value)) +
      stat_summary(
        fun = mean,
        geom = "point",
        size = 1, mapping = aes(shape = type)
      )

    if(names(met_it)[y] == "core") {
      p2 <- p2 +
        stat_summary(
          fun.min = ~ quantile(.x, 0.25),
          fun.max = ~ quantile(.x, 0.75),
          geom = "errorbar",
          width = 0.2, linewidth = 0.2
        )
    }

    p2 <- p2 +
      ylab("input data") +
      scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5)) +
      ggplot2::geom_point(data = dat_toplot2, aes(x = index, y = value, shape = type), size = 2) +
      ggplot2::facet_grid(rows = vars(name), scales = "free_y") +
      scale_shape_manual(
        values = c(
          "average of residues within or adjacent to known GPCR ligands" = 21,
          "average of residues from ~40K candidate regions" = 15
        )
      ) +
      theme_bw() +
      theme(
        strip.background = element_blank(),     # removes grey box
        strip.text = element_text(
          color = "black",
          size = 10
        )  )


    seq_dat <- nn_input_comb %>%
      filter(target %in% targ & known == 1) %>%
      mutate(meta_data = map2(meta_data, pep_id, \(x,y) x %>% mutate(metric = y) %>% mutate(index_og = index) %>% mutate(index = row_number()))) %>%
      pull(meta_data) %>%
      bind_rows


    seq_dat$metric <- factor(seq_dat$metric,
                             levels = rev(sort(unique(seq_dat$metric))))

    p3 <- ggplot2::ggplot(seq_dat) +
      ggplot2::geom_tile(mapping = aes(x = index, y = metric, fill = known_idx)) +
      ggiraph::geom_point_interactive(data = seq_dat %>%
                                        mutate(tt_value = "test"),
                                      aes(x = index, y = metric, tooltip = tt_value, data_id = index), pch = 15, size = 2.5, color = "grey85") +

      scale_fill_discrete(palette = ggsci::pal_simpsons()) +
      ylab("") +
      xlab("") +

      scale_x_continuous(expand = expansion(mult = c(0, 0)),
                         breaks = seq(0, max_index, by = 10),
                         minor_breaks = seq(0, max_index, by = 5)) +
      ggplot2::geom_text(aes(x = index, y = metric, label = AA), size = 1.8, fontface = "bold", color = "black") +
      theme_bw()





    p_final[[targ[1]]] <- list(p3, p2, p)



  }

  p_final <- unlist(p_final, recursive = FALSE)

  p_final2 <- Reduce(`+`, p_final) +  patchwork::plot_layout(ncol = 2, nrow = 3, byrow = FALSE, heights = c(0.6,1,1),
                                                             guides = "collect") &
    theme(legend.position = "top", legend.title = element_blank())

  ggsave(filename = paste0("~/AF2_analysis/shap_1dcnnddfd_", names(met_it)[y], ".svg"), p_final2, width = 16, height = 14)


}

















library(fastshap)

pred_fun <- function(object, newdata) {

  n <- nrow(newdata)

  reshaped <- array(
    newdata,
    dim = c(n, 36, 110)
  )

  preds <- predict(object, reshaped)

  as.numeric(preds$global[,1])
}


x_val_flat <- matrix(
  x_val,
  nrow = dim(x_val)[1],
  ncol = prod(dim(x_val)[-1])
)

baseline <- mean(pred_fun(models[["C"]], x_val_flat))

shap_values <- fastshap::explain(
  object = models[["C"]],
  X = x_val_flat,
  pred_wrapper = pred_fun,
  nsim = 100,
  baseline = baseline
)































