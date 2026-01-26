

###run 9_expand_by_residue.R


offlimits_features <- c("signal peptide", "E")

black_balled <- c("thrombin light chain", "epiregulin", "INSL5 (A chain)")


one_hot_encode <- function(x, tag = "") {
  x[is.na(x)] <- "unknown"
  unique_aa <- unique(x)
  to_return <- t(sapply(x, function(y) as.integer(unique_aa == y)))
  colnames(to_return) <- paste0(tag, "_", unique_aa)
  return(as_tibble(to_return))
}



AAs <- c("AA_Q", "AA_P", "AA_V", "AA_L", "AA_T", "AA_S", "AA_A",
         "AA_G", "AA_R", "AA_F", "AA_C", "AA_I", "AA_N", "AA_Y", "AA_W",
         "AA_K", "AA_D", "AA_E", "AA_M", "AA_H", "AA_U", "minAA_E", "minAA_A", "minAA_I",
         "minAA_V", "minAA_L", "minAA_P", "minAA_T", "minAA_S", "minAA_Q",
         "minAA_G", "minAA_D", "minAA_F", "minAA_N", "minAA_R", "minAA_H",
         "minAA_Y", "minAA_K", "minAA_M", "minAA_C", "minAA_W", "maxAA_C",
         "maxAA_D", "maxAA_K", "maxAA_W", "maxAA_P", "maxAA_F", "maxAA_M",
         "maxAA_E", "maxAA_H", "maxAA_L", "maxAA_I", "maxAA_G", "maxAA_Q",
         "maxAA_N", "maxAA_T", "maxAA_Y", "maxAA_R", "maxAA_S", "maxAA_A",
         "maxAA_V")

cons_metrics <- c("cons_rs", "cons_lrs", "blos_wt_all",
                  "blos_uw_all", "blos_wt_mam", "blos_uw_mam", "gran_wt_all", "gran_uw_all",
                  "gran_wt_mam", "gran_uw_mam", "blos_nr_all", "blos_nr_mam", "gran_nr_all",
                  "gran_nr_mam")

basic <- c(cons_metrics, paste0(cons_metrics, "_n"), c("max_afm", "min_afm", "mean_afm"), "relASA")

SS <- c("SS_P", "SS_S",
        "SS_E", "SS_-", "SS_T", "SS_G", "SS_B", "SS_H", "SS_I")

angles <- c("Phi_cos", "Psi_cos", "Phi_sin", "Psi_sin")

energy <- c("NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy")

sites <- c("DB_Dibasic",  "SV_sequence variant")






nn <- list()
nn[["nn1"]] <- c(basic, SS)
nn[["nn2"]] <- c(nn[["nn1"]], angles, energy)
nn[["nn3"]] <- c(nn[["nn2"]], sites)
nn[["nn4"]] <- c(nn[["nn3"]], AAs)


nn <- tibble(neural_net = names(nn),
             parameters = nn)

all_params <- unique(do.call(c, nn[["parameters"]]))





context <- c("_lag1", "_lag2", "_lead1", "_lead2")

nn2 <- nn %>%
  mutate(parameters = map(parameters, .f = ~c(., do.call(c, map(., ~paste0(., context)))))) %>%
  mutate(neural_net = paste0(neural_net, "c"))

names(nn2[["parameters"]]) <- paste0(names(nn2[["parameters"]]), "c")

nn <- bind_rows(nn, nn2)

nn <- nn %>%
  mutate(parameters = map(parameters, .f = ~c(., "known")))


if(FALSE) {
to_remove <- "^SS_.*_(lead|lag)\\d+$"

nn3 <- nn %>%
  filter(neural_net %in% c("nn5c", "nn6c", "nn7c")) %>%
  mutate(parameters = map(parameters, .f = ~.[!str_detect(., to_remove)]))

names(nn3$parameters) <- nn3$neural_net <- paste0(nn3$neural_net, "noSS")

nn <- bind_rows(nn, nn3)

  nn_og <- readRDS("~/peptide_alg/nn_models_256_64_4_leaky_relu_lion_opt_large_val_better_angles_w_context.rds")
  nn <- bind_cols(nn, nn_og[,3:5])
}


secretome_aa <- secretome_aa %>%
  mutate(SS = replace_na(SS, "-")) %>%
  mutate(one_hot_encode(AA, tag = "AA")) %>%
  mutate(one_hot_encode(SS, tag = "SS")) %>%
  mutate(one_hot_encode(sites_Dibasic_type, tag = "DB")) %>%
  mutate(one_hot_encode(`uniprot_sequence variant_type`, tag = "SV")) %>%
  mutate(one_hot_encode(gpcrdb_gtp_type, tag = "known")) %>%
  mutate(one_hot_encode(min_AA, tag = "minAA")) %>%
  mutate(one_hot_encode(max_AA, tag = "maxAA")) %>%
  mutate(across(c("Phi", "Psi"), .fns = ~tibble(cos = cos(. * pi/180),
                                                sin = sin(. * pi/180)), .unpack = TRUE)) %>%
  dplyr::group_by(accession) %>%
  mutate(across(ends_with("_energy"), .fns = \(x) {
    x[is.na(x)] <- 0
    scales::rescale(x, to = c(0, 1))})) %>%
  ungroup() %>%
  select(!ends_with("_unknown")) %>%
  mutate(known = rowSums(across(starts_with("known_gpcr_pep"))) > 0) %>%
  mutate(across(all_of(cons_metrics), ~ .x * species_limit, .names = "{.col}_n"))


secretome_aa <- secretome_aa %>%
  group_by(accession) %>%
  mutate(across(all_of(all_params),
                .fns = list(lag1 = ~lag(., n = 1),
                            lag2 = ~lag(., n = 2),
                            lag3 = ~lag(., n = 3),
                            lag4 = ~lag(., n = 4),
                            lead1 = ~lead(., n = 1),
                            lead2 = ~lead(., n = 2),
                            lead3 = ~lead(., n = 3),
                            lead4 = ~lead(., n = 4)),
                .unpack = TRUE)) %>%
  ungroup()





known_genes <- secretome_aa %>%
  filter(1 %in% known, .by = gene) %>%
  mutate(train = if_else(grepl("^(CCL|CXCL|XCL|CX3CL)", gene), "chem", "pep")) %>%
  select(gene, train) %>%
  split(.[["train"]]) %>%
  lapply(., \(x) unique(x[["gene"]]))

ctrl_genes <- secretome_aa %>%
  filter(!1 %in% known, .by = gene) %>%
  pull(gene) %>%
  unique(.)

ctrl_genes_sub <- sample(ctrl_genes, size = 100)

known_samp <- sample(c(TRUE, FALSE),
                          length(known_genes),
                          replace=TRUE,
                          prob=c(0.5,0.5))

ctrl_samp <- sample(c(TRUE, FALSE),
                    length(ctrl_genes_sub),
                    replace=TRUE,
                    prob=c(0.5,0.5))


train <- secretome_aa[secretome_aa[["accession"]] %in% c(gtp_accession[gtp_samp], ctrl_accession[ctrl_samp]), ]

validate <- secretome_aa[secretome_aa[["accession"]] %in% c(gtp_accession[!gtp_samp], ctrl_accession[!ctrl_samp]), ]

train <- secretome_aa[secretome_aa[["accession"]] %in% c(gtp_accession[gtp_samp]), ]

validate <- secretome_aa[secretome_aa[["accession"]] %in% c(gtp_accession[!gtp_samp]), ]





library(keras3)

weighted_binary_crossentropy <- function(weight_0, weight_1) {
  function(y_true, y_pred) {
    - (weight_1 * y_true * log(y_pred) + weight_0 * (1 - y_true) * log(1 - y_pred))
  }
}



for(x in nn[["neural_net"]]) {
  for(tag in names(known_genes)) {

    message("computing ", tag, " for ", x)

    start <- Sys.time()

    x_train <- secretome_aa %>%
      filter(gene %in% c(known_genes[[tag]], ctrl_genes_sub)) %>%
      select(nn[["parameters"]][[x]]) %>%
      drop_na %>%
      select(-known) %>%
      as.matrix

    y_train <- secretome_aa %>%
      filter(gene %in% c(known_genes[[tag]], ctrl_genes_sub)) %>%
      select(nn[["parameters"]][[x]]) %>%
      drop_na %>%
      mutate(known = as.integer(known)) %>%
      pull(known)

    # x_val <- validate %>%
    #   select(nn[["parameters"]][[x]]) %>%
    #   drop_na %>%
    #   select(-known) %>%
    #   as.matrix
    #
    # y_val <- validate %>%
    #   select(nn[["parameters"]][[x]]) %>%
    #   drop_na %>%
    #   mutate(known = as.integer(known)) %>%
    #   pull(known)

    class_counts <- table(y_train)

    class_weights <- list(
      "0" = class_counts[2] / class_counts[1],  # Weight for majority class
      "1" = class_counts[2] / class_counts[2]   # Weight for minority class
    )
    rm(model)

    early_stopping <- callback_early_stopping(
      monitor = "val_binary_accuracy",
      patience = 40,
      mode = 'max',
      start_from_epoch = 20,
      restore_best_weights = TRUE
    )


    model <- keras_model_sequential() %>%
      layer_batch_normalization() %>%
      layer_dense(units = 256,
                  activation = "leaky_relu",
                  kernel_regularizer = regularizer_l2(0.01)) %>%
      layer_batch_normalization() %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 64,
                  activation = "leaky_relu",
                  kernel_regularizer = regularizer_l2(0.01)) %>%
      layer_batch_normalization() %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 4,
                  activation = "leaky_relu",
                  kernel_regularizer = regularizer_l2(0.01)) %>%
      layer_batch_normalization() %>%
      layer_dropout(rate = 0.5) %>%
      layer_dense(units = 1,
                  activation = "sigmoid",
                  kernel_regularizer = regularizer_l2(0.01))

    # Compile the model
    model %>% compile(
      optimizer = optimizer_adam(learning_rate = 0.0001),
      loss = loss_binary_crossentropy(),
      #loss = weighted_binary_crossentropy(weight_0 = class_weights[["0"]], weight_1 = class_weights[["1"]]),
      metrics = c(metric_binary_accuracy(),
                  metric_specificity_at_sensitivity(sensitivity = 0.8),
                  metric_sensitivity_at_specificity(specificity = 0.8),
                  metric_auc())
    )


    # Fit the model
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 100,
      class_weight = class_weights,
      validation_split = 0.4,
      #validation_data = list(x_val, y_val),
      callbacks = list(early_stopping),
      verbose = 1  # Set to 1 to view progress
    )

    end <- Sys.time()
    comp_time <- end - start

    message("compute time ", comp_time)

    nn[nn$neural_net == x, "comp_time"] <- comp_time

    nn[nn$neural_net == x, "model"] <- list(list(model))

    nn[nn$neural_net == x, "history"] <- list(list(plot(history)))


    secretome_aa[[paste0(tag, "_",x)]] <- model %>%
      predict(secretome_aa %>%
                select(nn[["parameters"]][[x]]) %>%
                select(-known) %>%
                as.matrix)
}
}

message(paste0(paste(nn$neural_net, "comp time: ", round(nn$comp_time, 1)), collapse = "\n"))



params <- c("blos_wt_mam", "blos_wt_all_n", "relASA", "max_afm", "min_afm", "mean_afm",
            "NH->O_1_energy", "O->NH_1_energy", "NH->O_2_energy", "O->NH_2_energy",
            "Phi_cos", "Phi_sin", "Psi_cos", "Psi_sin")
#smooth scores


secretome_aa <- secretome_aa %>%
  mutate(across(starts_with(c("pep_nn", "pep_xgb", "chem_nn", "chem_xgb")), .fns = ~scales::rescale(., to = c(0,1)), .unpack = TRUE)) %>%
  group_by(accession) %>%
  mutate(across(starts_with(c("pep_nn", "pep_xgb", "chem_nn", "chem_xgb")), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE)) %>%
  mutate(across(all_of(params), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE)) %>%
  ungroup()


saveRDS(secretome_aa, paste0(s_localDir, "/processed/secretome_aa.rds"))











library(pROC)


cutoff <- 0.5

roc_res <- Map(\(x) {

  roc(secretome_aa$known_peptide, ifelse(secretome_aa[[x]] > cutoff, 1, 0), ci = TRUE)

}, grep("^score_nn", colnames(secretome_aa), value = TRUE))

names(roc_res) <- paste(names(roc_res), "\n[AUC:", round(sapply(roc_res, auc), 4), "]")


extra_color <- c("#FD7446FF","#FD8CC1FF","#FED439FF", "#197EC0FF", "#46732EFF", "#C80813FF",
                 "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF",
                 "#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF")


ggroc(roc_res, linewidth = 1) + ggplot2::theme_bw() +
  scale_color_discrete(name = "", type = extra_color) +
  ggtitle(paste("ROC Analysis"), subtitle = cutoff) +
  theme(
    legend.key.height = unit(3, "lines")
  )

roc_res <- Map(\(x) {

  roc(secretome_aa %>%
        select(nn[["parameters"]][[x]]) %>%
        pull(known_peptide),
      predict(nn[["model_known_peptide"]][nn[["neural_net"]] == x][[1]],
              secretome_aa %>%
                select(nn[["parameters"]][[x]]) %>%
                select(-all_of(paste0("known_peptide", c("", "_n", "_c")))) %>%
                as.matrix), ci = TRUE)

}, nn$neural_net[!sapply(nn$model_known_peptide, is.null)])

names(roc_res) <- paste(names(roc_res), "\n[AUC:", round(sapply(roc_res, auc), 4), "]")

nn[["AUC"]] <- round(sapply(roc_res, auc), 4)


ggroc(roc_res, linewidth = 1) + ggplot2::theme_bw() +
  scale_color_discrete(name = "", type = extra_color) +
  ggtitle(paste("ROC Analysis"), subtitle = cutoff) +
  theme(
    legend.key.height = unit(3, "lines")
  )
names(roc_res)


saveRDS(nn, "~/peptide_alg/nn_models_256_64_4_leaky_relu_lion_opt_large_val_better_angles_w_context_NC.rds")


to_export <- nn %>%
  select(-model, -training) %>%
  mutate(parameters = map_chr(parameters, .f = ~paste0(., collapse = "\n")))

openxlsx::write.xlsx(to_export, "~/peptide_alg/nn_parameters.xlsx")


nn_tp <- tibble(score_name = secretome_aa %>% select(starts_with("score_nn")) %>% colnames(.),
                nn_type = str_extract(score_name, "nn\\d{1,2}"),
                tmp = str_remove(score_name, paste0("score_", nn_type)))

nn_tp <- nn_tp %>%
  mutate(nn_type = ifelse(str_detect(tmp, "^c"), paste0(nn_type, "c"), nn_type)) %>%
  mutate(nn_type2 = ifelse(str_detect(tmp, "^c"), "context", "isolated")) %>%
  mutate(tmp = str_remove(tmp, "^c")) %>%
  mutate(model_type = str_extract(tmp, "^.")) %>%
  mutate(model_type = ifelse(model_type == "_", "FL", model_type))


nns <- paste0("nn", c(1,3,9,10))

nn_tp <- nn_tp %>%
  filter(nn_type %in% c(nns, paste0(nns, "c"))) %>%
  dplyr::rename(params = score_name)



params_tp <- c("AF_Phi_cos", "AF_Psi_cos", "AF_Phi_sin", "AF_Psi_sin", "AF_NH->O_1_energy", "AF_O->NH_1_energy", "AF_NH->O_2_energy", "AF_O->NH_2_energy")

