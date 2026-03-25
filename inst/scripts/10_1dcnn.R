



library(keras3)
library(tensorflow)
library(yardstick)




seq_len <- nn_input$N$train$data %>% `[[`(1) %>% nrow
n_channels <- nn_input$N$train$data %>% `[[`(1) %>% ncol

inputs <- layer_input(shape = c(seq_len, n_channels))

# Shared convolutional backbone
regularizer <- regularizer_l2(1e-4)

x <- inputs |>
  layer_conv_1d(filters = 32, kernel_size = 5, padding = "same", activation = "relu", kernel_regularizer = regularizer) |>
  layer_dropout(0.3) |>
  layer_conv_1d(filters = 64, kernel_size = 5, padding = "same", activation = "relu", kernel_regularizer = regularizer) |>
  layer_dropout(0.3)

# -------------------------
# Output 1: Per-index prediction
# -------------------------
per_index_output <- x |>
  layer_conv_1d(filters = 1, kernel_size = 1, activation = "sigmoid",
                name = "per_index")

# Shape: (batch, seq_len, 1)

# -------------------------
# Output 2: Global prediction
# -------------------------
global_output <- x |>
  layer_global_average_pooling_1d() |>
  layer_dense(1, activation = "sigmoid", name = "global")

# Shape: (batch, 1)

model <- keras_model(
  inputs = inputs,
  outputs = list(global = global_output,
                 per_index = per_index_output)
)

model


model |> compile(
  optimizer = optimizer_adam(learning_rate = 1e-4),
  loss = list(
    global = "binary_crossentropy",
    per_index = "binary_crossentropy"
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

models <- list(N = model,
               C = model)

val_pred <- list()

for(term in c("N", "C")) {

nn_in <- nn_input[[term]]

n_samples <- length(nn_in$train$data)

x_train <- array(
  unlist(nn_in$train$data),
  dim = c(seq_len, n_channels, length(nn_in$train$data))
  )

x_train <- aperm(x_train, c(3, 1, 2))


y_train_global <- matrix(nn_in$train$known,
                   ncol = 1)

y_train_per_index <- array(rep(nn_in$train$known, seq_len),
                     dim = c(n_samples, seq_len, 1))




n_samples <- length(nn_in$val$data)

x_val <- array(
  unlist(nn_in$val$data),
  dim = c(seq_len, n_channels, length(nn_in$val$data))
)

x_val <- aperm(x_val, c(3, 1, 2))


y_val_global <- matrix(nn_in$val$known,
                   ncol = 1)

y_val_per_index <- array(rep(nn_in$val$known, seq_len),
                     dim = c(n_samples, seq_len, 1))






models[[term]] |> fit(
  x = x_train,
  y = list(
    global = y_train_global,
    per_index = y_train_per_index
  ),

  epochs = 500,
  batch_size = 32
)


val_pred[[term]] <- models[[term]] %>%
              predict(x_val) %>%
              `[[`("global")

}



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

View(nn_input_comb %>% arrange(desc(pred)) %>% mutate(rank = row_number()))


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


targ <- c("C", "loop_C")
targ <- c("N", "loop_N")

shap_vals <- saliency_fn(models[[targ[1]]], x_val)

mets <- all_params3



dat_toplot <- lapply(mets, \(x) shap_vals[y_val_global == 1,,which(all_params3 == x)] %>% as_tibble %>% mutate(metric = x)) %>%
                bind_rows() %>%
                  pivot_longer(cols = -metric) %>%
                  mutate(index = setNames(1:36, paste0("V", 1:36))[name]) %>%
                  select(-name)

dat_toplot2 <- map(known_dat$data[known_dat$target %in% targ], \(x) {

  x[, mets] %>%
    mutate(index = row_number())

}) %>%
  bind_rows %>%
  group_by(index) %>%
  summarise(across(everything(), mean, na.rm = TRUE),
            index = first(index)) %>%
  pivot_longer(cols = -index)

dat_toplot3 <- map(c_dat$data[c_dat$target %in% targ], \(x) {

  x[, mets] %>%
    mutate(index = row_number())

}) %>%
  bind_rows %>%
  pivot_longer(cols = -index)



p <- ggplot2::ggplot(dat_toplot %>% group_by(index, metric) %>%
                       summarise(across(everything(), mean, na.rm = T))) +
  #ggplot2::geom_violin(aes(x = name, y = value), trim = TRUE) +
  ggplot2::geom_point(aes(x = index, y = value)) +
  ggplot2::facet_grid(rows = vars(metric), scales = "free_y") +
  theme_bw()


p2 <- ggplot2::ggplot(data = dat_toplot3, aes(x = index, y = value)) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 1, pch = 15,
  ) +
  stat_summary(
    fun.min = ~ quantile(.x, 0.25),
    fun.max = ~ quantile(.x, 0.75),
    geom = "errorbar",
    width = 0.2
  ) +
  ggplot2::geom_point(data = dat_toplot2, aes(x = index, y = value), size = 2) +
  ggplot2::facet_grid(rows = vars(name), scales = "free_y") +
  theme_bw()

p_final <- p / p2


ggsave(filename = paste0("~/AF2_analysis/shap_1dcnn_", targ[1], "_test2.svg"), p_final, width = 6, height = 10)















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































