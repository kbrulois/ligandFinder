
secretome <- readRDS("~/peptide_alg/secretome_5.rds")


combine_scores <- function(scoreA, scoreB) {
  penalty <- ifelse(scoreB < 0.75, (1.6*(0.75 - scoreB))^2,
                    ifelse(scoreB > 0.9, (1.6*(scoreB - 0.9))^2, 0))

  combined_score <- scoreA * (1 - penalty)
  return(combined_score)
}


ss_features <- c(setNames("Alpha helix (4-12)", "H"),
                 setNames("Isolated beta-bridge residue", "B"),
                 setNames("Strand", "E"),
                 setNames("3-10 helix", "G"),
                 setNames("Pi helix", "I"),
                 setNames("Turn", "T"),
                 setNames("Bend", "S"),
                 setNames("Kappa helix", "P"))

ss_features_coded <- 1:length(ss_features)
names(ss_features_coded) <- names(ss_features)
ss_features_coded <- c(setNames(0, "-"), ss_features_coded)


offlimits_features <- c("signal peptide", "E")

black_balled <- c("thrombin light chain", "epiregulin", "INSL5 (A chain)")

to_expand <- c("uniprot", "sites")

factorize_intervals <- function(type, start, end, dat_size) {
  to_return <- rep(NA, dat_size)
  for(i in seq_along(type)) {
    start_na <- is.na(start[i]) | start[i] > dat_size
    end_na <- is.na(end[i]) | end[i] > dat_size
    if(!start_na & end_na) {
      to_return[start[i]] <- type[i]
    } else if(!start_na & !end_na) {
      to_return[start[i]:end[i]] <- type[i]
    }
  }
  return(to_return)
}

expand_by_AA <- \(sequence_uni, accession, gene, af_mapped, features, cons_mapped, af_missense_mapped, af_xyz_mapped) {

  aa_sequence <- strsplit(sequence_uni, "")[[1]]

  dat_size <- length(aa_sequence)

  if("frequency" %in% names(cons_mapped$ms)) {
    conservation <- cons_mapped$ms$frequency
  } else {
    conservation <- rep(NA, dat_size)
  }

  if("relASA" %in% names(af_mapped$ms)) {
    asa <- af_mapped$ms %>%
      select(-AA, -index, -relASA_s, -relASA_ss, -relASA_sss) %>%
      dplyr::rename_with(\(x) {
        ifelse(!x %in% c("SS", "relASA"), paste0("AF_", x), x)
      })
  } else {
    asa <- Map(\(x) {rep(NA, dat_size)},
               c("SS", "relASA", "AF_Phi", "AF_Psi", "AF_NH->O_1_relidx", "AF_NH->O_1_energy",
                 "AF_O->NH_1_relidx", "AF_O->NH_1_energy", "AF_NH->O_2_relidx",
                 "AF_NH->O_2_energy", "AF_O->NH_2_relidx", "AF_O->NH_2_energy")) %>%
      as_tibble
  }

  if("mean_af_missense" %in% names(af_missense_mapped[["ms"]])) {
    afm <- af_missense_mapped$ms$mean_af_missense
  } else {
    afm <- rep(NA, dat_size)
  }

  if("AA" %in% names(af_xyz_mapped[["ms"]])) {
    afxyz <- af_xyz_mapped[["ms"]] %>% select(-AA)
  } else {
    afxyz <- rep(NA, dat_size)
  }

  dat <- tibble(AA = aa_sequence,
                conservation = as.numeric(conservation),
                pathogenicity = as.numeric(afm),
                accession = accession,
                gene = gene,
                known_peptide = 0,
                known_peptide_n = 0,
                known_peptide_c = 0,
                signal_peptide_or_Strand = 0,
                dibasic_Cys_W_T = as.character(NA))

  dat <- bind_cols(dat, asa, afxyz)

  gtp <- features %>%
    filter(source == "gtp") %>%
    filter(!type %in% black_balled)

  if(nrow(gtp) > 0) {
    gtp_pep <- gtp %>%
      rowwise() %>%
      mutate(locations = map2(start, end, .f = \(x, y) x:y)) %>%
      pull(locations) %>%
      do.call(c, .) %>%
      unique(.)

  dat[["known_peptide"]][gtp_pep] <- 1
  dat[["known_peptide_n"]][gtp[["start"]]] <- 1
  dat[["known_peptide_c"]][gtp[["end"]]] <- 1
  }

  offlimits <- features %>%
    filter(type == "signal peptide" | type == "E") %>%
    rowwise() %>%
    mutate(locations = map2(start, end, .f = \(x, y) {if(!is.na(x) & !is.na(y)) x:y})) %>%
    pull(locations) %>%
    do.call(c, .) %>%
    unique(.)

  offlimits <- offlimits[offlimits %in% 1:nrow(dat)]

  dat[["signal_peptide_or_Strand"]][offlimits] <- 1

  for(feat in to_expand) {
    dat[[feat]] <- features %>%
      filter(source == feat) %>%
      {factorize_intervals(.[["type"]], .[["start"]], .[["end"]], dat_size)}
  }

  dat

}



start <- Sys.time()

secretome_aa <- secretome %>%
  rowwise %>%
  reframe(expand_by_AA(sequence_uni, accession, gene, af_mapped, features, cons_mapped, af_missense_mapped, af_xyz_mapped))

end <- Sys.time()
end - start




one_hot_encode <- function(x, tag = "") {
  x[is.na(x)] <- "unknown"
  unique_aa <- unique(x)
  to_return <- t(sapply(x, function(y) as.integer(unique_aa == y)))
  colnames(to_return) <- paste0(tag, "_", unique_aa)
  return(as_tibble(to_return))
}



secretome_aa <- secretome_aa %>%
  mutate(conservation = scales::rescale(log(conservation + 0.1), c(1,0))) %>%
  mutate(score = combine_scores(conservation, relASA)) %>%
  mutate(alpha_helix = ifelse(SS == "H", 1, 0)) %>%
  mutate(SS = replace_na(SS, "-")) %>%
  mutate(integerSS = ss_features_coded[SS]) %>%
  mutate(one_hot_encode("AA", tag = "AA")) %>%
  mutate(one_hot_encode(SS, tag = "SS")) %>%
  mutate(one_hot_encode(sites, tag = "sites")) %>%
  mutate(across(c("AF_Phi", "AF_Psi"), .fns = ~tibble(cos = cos(. * pi/180),
                                                      sin = sin(. * pi/180)), .unpack = TRUE)) %>%
  dplyr::group_by(accession) %>%
  mutate(conservation_norm = conservation/mean(conservation, na.rm = TRUE)) %>%
  dplyr::rename(conservation_og = conservation) %>%
  mutate(across(ends_with("_energy"), .fns = ~scales::rescale(., to = c(0, 1)))) %>%
  ungroup()


secretome_aa_og <- readRDS("~/peptide_alg/secretome_aa_Dec11.rds")
secretome_aa[["nn_score"]] <- secretome_aa_og[["score_nn7c_noSS_s8"]]

derivative <- function(x) {
  subsetter <- !is.na(x) & !is.nan(x) & !is.infinite(x)
  to_return <- as.numeric(rep(NA, length(x)))
  to_return[subsetter] <- tryCatch({predict(smooth.spline(x[subsetter]), deriv = 1)[["y"]]},
                                   error = function(e) {as.numeric(rep(NA, sum(subsetter)))})
  return(to_return)
}

secretome_aa <- secretome_aa %>%
  group_by(accession) %>%
  mutate(rate_of_change_score = derivative(nn_score)) %>%
  mutate(rate_of_change_score = scales::rescale(rate_of_change_score, to = c(-1, 1))) %>%
  ungroup


ggplot2::ggplot(secretome_aa %>% mutate(known_peptide2 = ifelse(known_peptide == 1, "known", "unknown"))) +
  ggplot2::geom_boxplot(aes(x = known_peptide2, y = `AF_Psi_sin`))



basic <- c("conservation_norm", "pathogenicity", "relASA", "signal_peptide_or_Strand")
angles <- c("AF_Phi_cos", "AF_Psi_cos", "AF_Phi_sin", "AF_Psi_sin")

nn <- list()
nn[["nn1"]] <- basic
nn[["nn2"]] <- c("conservation_norm", "pathogenicity", "relASA", "signal_peptide_or_Strand")
nn[["nn3"]] <- c(nn[["nn1"]], "^SS_")
nn[["nn4"]] <- c(nn[["nn1"]], angles)
nn[["nn5"]] <- c(nn[["nn3"]], angles)
nn[["nn6"]] <- c(nn[["nn3"]], "_energy$")
nn[["nn7"]] <- c(nn[["nn5"]], "_energy$")
nn[["nn8"]] <- c(nn[["nn7"]], "nn_score")
nn[["nn9"]] <- c(nn[["nn8"]], "^sites_")
nn[["nn10"]] <- c(nn[["nn9"]], "rate_of_change_score")


nn <- Map(\(x) do.call(c, lapply(x, \(y) grep(y, colnames(secretome_aa), value = TRUE))), nn)

nn <- tibble(neural_net = names(nn),
             parameters = nn)

all_params <- unique(do.call(c, nn[["parameters"]]))

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


context <- c("_lag1", "_lag2", "_lead1", "_lead2")

nn2 <- nn %>%
      mutate(parameters = map(parameters, .f = ~c(., do.call(c, map(., ~paste0(., context)))))) %>%
      mutate(neural_net = paste0(neural_net, "c"))

names(nn2[["parameters"]]) <- paste0(names(nn2[["parameters"]]), "c")

nn <- bind_rows(nn, nn2)

nn <- nn %>%
      mutate(parameters = map(parameters, .f = ~c(., paste0("known_peptide", c("", "_n", "_c")))))


to_remove <- "^SS_.*_(lead|lag)\\d+$"

nn3 <- nn %>%
        filter(neural_net %in% c("nn5c", "nn6c", "nn7c")) %>%
        mutate(parameters = map(parameters, .f = ~.[!str_detect(., to_remove)]))

names(nn3$parameters) <- nn3$neural_net <- paste0(nn3$neural_net, "noSS")

nn <- bind_rows(nn, nn3)

if(FALSE) {
  nn_og <- readRDS("~/peptide_alg/nn_models_256_64_4_leaky_relu_lion_opt_large_val_better_angles_w_context.rds")
  nn <- bind_cols(nn, nn_og[,3:5])
}

gtp_accession <- secretome_aa %>%
  filter(1 %in% known_peptide, .by = accession) %>%
  pull(accession) %>%
  unique(.)

ctrl_accession <- secretome_aa %>%
  filter(!accession %in% secretome[["accession"]][secretome[["non_secreted_goi"]]]) %>%
  filter(!1 %in% known_peptide, .by = accession) %>%
  pull(accession) %>%
  unique(.)

gtp_samp <- sample(c(TRUE, FALSE),
                   length(gtp_accession),
                   replace=TRUE,
                   prob=c(0.5,0.5))

ctrl_samp <- sample(c(TRUE, FALSE),
                    length(ctrl_accession),
                    replace=TRUE,
                    prob=c(0.05,0.95))


train <- secretome_aa[secretome_aa[["accession"]] %in% c(gtp_accession[gtp_samp], ctrl_accession[ctrl_samp]), ]

validate <- secretome_aa[secretome_aa[["accession"]] %in% c(gtp_accession[!gtp_samp], ctrl_accession[!ctrl_samp]), ]


#
# models <- list(conservation = "conservation",
#                relASA = "relASA",
#                pathogenicity = "pathogenicity",
#                "conservation +\nrelASA +\npathogenicity" = c("conservation", "relASA", "pathogenicity"))
#
# formula <- as.formula(paste("known_peptide", "~", paste(models[[4]], collapse = " + ")))
#
# subsetter <- complete.cases(secretome_aa %>% select(conservation, relASA, pathogenicity, known_peptide))
#
# model_glm <- glm(formula, family="binomial", data=secretome_aa[subsetter, ])
#
# secretome_aa$score_glm <- as.numeric(NA)
#
# secretome_aa$score_glm[subsetter] <- predict(model_glm, type="response")
#



library(keras3)

weighted_binary_crossentropy <- function(weight_0, weight_1) {
  function(y_true, y_pred) {
    - (weight_1 * y_true * log(y_pred) + weight_0 * (1 - y_true) * log(1 - y_pred))
  }
}

targets <- paste0("known_peptide", c("", "_n", "_c"))

for(x in nn[["neural_net"]]) {
  for(y in targets) {

  message("computing ", x, " ", y)

  start <- Sys.time()

  x_train <- train %>%
    select(nn[["parameters"]][[x]]) %>%
    drop_na %>%
    select(-starts_with("known_peptide")) %>%
    as.matrix

  y_train <- train %>%
    select(nn[["parameters"]][[x]]) %>%
    drop_na %>%
    pull(!!sym(y))

  x_val <- validate %>%
    select(nn[["parameters"]][[x]]) %>%
    drop_na %>%
    select(-starts_with("known_peptide")) %>%
    as.matrix

  y_val <- validate %>%
    select(nn[["parameters"]][[x]]) %>%
    drop_na %>%
    pull(!!sym(y))

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
                activation = "leaky_relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 64,
                activation = "leaky_relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 8,
                activation = "leaky_relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 1,
                activation = "softmax")

  # Compile the model
  model %>% compile(
    optimizer = optimizer_lion(learning_rate = 0.00001),
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
    #validation_split = 0.4,
    validation_data = list(x_val, y_val),
    callbacks = list(early_stopping),
    verbose = 1  # Set to 1 to view progress
  )

  end <- Sys.time()
  comp_time <- end - start

  message("compute time ", comp_time)

  nn[nn$neural_net == x, paste0("comp_time", "_", y)] <- comp_time

  nn[nn$neural_net == x, paste0("model", "_", y)] <- list(list(model))

  nn[nn$neural_net == x, paste0("history", "_", y)] <- list(list(plot(history)))


  secretome_aa[[paste0("score_", x, "_", sub("known_peptide", "", y))]] <- model %>%
    predict(secretome_aa %>%
              select(nn[["parameters"]][[x]]) %>%
              select(-all_of(paste0("known_peptide", c("", "_n", "_c")))) %>%
              as.matrix)
  }
}

message(paste0(paste(nn$neural_net, "comp time: ", round(nn$comp_time, 1)), collapse = "\n"))



colnames(secretome_aa) <- sub("__n", "N", colnames(secretome_aa))
colnames(secretome_aa) <- sub("__c", "C", colnames(secretome_aa))


#for(x in nn[["neural_net"]]) {

  secretome_aa[[paste0("score_", x)]] <- nn[["model"]][nn[["neural_net"]] == x][[1]] %>%
    predict(secretome_aa %>%
              select(nn[["parameters"]][[x]]) %>%
              select(-all_of(paste0("known_peptide", c("", "_n", "_c")))) %>%
              as.matrix)

}


params <- c("conservation_og", "conservation_norm", "relASA", "pathogenicity",
            "AF_NH->O_1_energy", "AF_O->NH_1_energy", "AF_NH->O_2_energy", "AF_O->NH_2_energy",
            "AF_Phi_cos", "AF_Phi_sin", "AF_Psi_cos", "AF_Psi_sin")
#smooth scores

win_sizes <- c(4,6,8)

smoother_func <- function(x, window_size = win_sizes, append_name = "a") {
  bind_cols(
    lapply(window_size, \(y) {
      tibble(!!paste0(append_name, y) := slider::slide_dbl(x, ~mean(., na.rm = TRUE), .before = y %/% 2, .after = y %/% 2))
    }
    ))
}

secretome_aa <- secretome_aa %>%
  mutate(across(starts_with("score_nn"), .fns = ~scales::rescale(., to = c(0,1)), .unpack = TRUE)) %>%
  group_by(accession) %>%
  mutate(across(score, .fns = smoother_func, .unpack = TRUE)) %>%
  mutate(across(all_of(params), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE))


for(win_size in win_sizes) {
  secretome_aa <- secretome_aa %>%
    group_by(accession) %>%
    mutate(!!paste0("score_b", win_size) := combine_scores(!!sym(paste0("conservation_og_s", win_size)),
                                                           !!sym(paste0("relASA_s", win_size))))
}

secretome_aa <- secretome_aa %>%
  group_by(accession) %>%
  mutate(across(starts_with("score_nn"), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE))



secretome_aa <- secretome_aa %>%
  ungroup

secretome_aa_og <- readRDS("~/peptide_alg/secretome_aa_DEC19.rds") %>%
                      select(starts_with("score_nn7c_noSS"))

colnames(secretome_aa_og) <- sub("_noSS", "noSS", colnames(secretome_aa_og))

secretome_aa <- bind_cols(secretome_aa, secretome_aa_og)


saveRDS(secretome_aa, "~/peptide_alg/secretome_aa.rds")



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




