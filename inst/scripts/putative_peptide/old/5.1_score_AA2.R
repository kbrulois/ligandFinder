
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


expand_by_AA <- \(sequence_uni, accession, gene, af_mapped, features, cons_mapped, af_missense_mapped) {
  
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
  
  if("mean_af_missense" %in% names(af_missense_mapped$ms)) {
    afm <- af_missense_mapped$ms$mean_af_missense
  } else {
    afm <- rep(NA, dat_size)
  }
  
  dat <- tibble(AA = aa_sequence,
                conservation = as.numeric(conservation),
                pathogenicity = as.numeric(afm),
                accession = accession,
                gene = gene,
                known_peptide = 0,
                signal_peptide_or_Strand = 0,
                dibasic_Cys_W_T = as.character(NA))
  
  dat <- bind_cols(dat, asa)
  
  gtp <- features %>%
    filter(source == "gtp") %>%
    filter(!type %in% black_balled)
  
  if(nrow(gtp) > 0) {
    gtp <- gtp %>%
      rowwise() %>%
      mutate(locations = map2(start, end, .f = \(x, y) x:y)) %>%
      pull(locations) %>%
      do.call(c, .) %>%
      unique(.)
  } else {
    gtp <- NULL
  }
  
  dat[["known_peptide"]][gtp] <- 1
  
  
  offlimits <- features %>%
    filter(type == "signal peptide" | type == "E") %>%
    rowwise() %>%
    mutate(locations = map2(start, end, .f = \(x, y) {if(!is.na(x) & !is.na(y)) x:y})) %>%
    pull(locations) %>%
    do.call(c, .) %>%
    unique(.)
  
  offlimits <- offlimits[offlimits %in% 1:nrow(dat)]
  
  dat[["signal_peptide_or_Strand"]][offlimits] <- 1
  
  factorize_intervals <- function(type, start, end) {
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
  
  for(feat in to_expand) {
    dat[[feat]] <- features %>%
      filter(source == feat) %>%
      {factorize_intervals(.[["type"]], .[["start"]], .[["end"]])}
  }
  
  dat
  
}



start <- Sys.time()

secretome_aa <- secretome %>%
  rowwise %>%
  reframe(expand_by_AA(sequence_uni, accession, gene, af_mapped, features, cons_mapped, af_missense_mapped))

end <- Sys.time()
end - start



#data.table::fwrite(secretome_aa, "~/peptide_alg/per_aa_secretome.csv")


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
  mutate(one_hot_encode(`AA`, tag = "AA")) %>%
  mutate(one_hot_encode(SS, tag = "SS")) %>%
  mutate(AF_Phi = scales::rescale(AF_Phi, from = c(-360, 360), to = c(-1, 1))) %>%
  mutate(AF_Psi = scales::rescale(AF_Psi, from = c(-360, 360), to = c(-1, 1))) %>%
  dplyr::group_by(accession) %>%
  mutate(conservation_norm = conservation/mean(conservation, na.rm = TRUE)) %>%
  ungroup()



ggplot2::ggplot(secretome_aa %>% mutate(known_peptide2 = ifelse(known_peptide == 1, "known", "unknown"))) +
  ggplot2::geom_boxplot(aes(x = known_peptide2, y = `AF_Psi`)) 





gtp_accession <- secretome_aa %>%
  filter(1 %in% known_peptide, .by = accession) %>%
  pull(accession) %>%
  unique(.)

ctrl_accession <- secretome_aa %>%
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
#                `conservation +\nrelASA +\npathogenicity` = c("conservation", "relASA", "pathogenicity"))
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

nn <- list()
nn[["nn1"]] <- c("conservation", "pathogenicity", "relASA", "signal_peptide_or_Strand")
nn[["nn2"]] <- c("conservation_norm", "pathogenicity", "relASA", "signal_peptide_or_Strand")
nn[["nn3"]] <- c(nn[["nn2"]], "alpha_helix")
nn[["nn4"]] <- c(nn[["nn2"]], "integerSS")
nn[["nn5"]] <- c(nn[["nn2"]], "^SS_")
nn[["nn6"]] <- c(nn[["nn2"]], "AF_Phi", "AF_Psi")
nn[["nn7"]] <- c(nn[["nn5"]], "AF_Phi", "AF_Psi")
nn[["nn8"]] <- c(nn[["nn5"]], "^AF_")

nn <- Map(\(x) do.call(c, lapply(x, \(y) grep(y, colnames(secretome_aa), value = TRUE))), nn)

nn <- Map(\(x) c(x, "known_peptide"), nn)

nn <- tibble(neural_net = names(nn),
             parameters = nn)

weighted_binary_crossentropy <- function(weight_0, weight_1) {
  function(y_true, y_pred) {
    - (weight_1 * y_true * log(y_pred) + weight_0 * (1 - y_true) * log(1 - y_pred))
  }
}

for(x in nn$neural_net) {
  
  message("computing ", x)
  
  start <- Sys.time()
  
  x_train <- train %>% 
    select(nn[["parameters"]][[x]]) %>% 
    drop_na %>% 
    select(-known_peptide) %>%
    as.matrix
  
  y_train <- train %>% 
    select(nn[["parameters"]][[x]]) %>% 
    drop_na %>% 
    pull(known_peptide) 
  
  x_val <- validate %>% 
    select(nn[["parameters"]][[x]]) %>% 
    drop_na %>% 
    select(-known_peptide) %>%
    as.matrix
  
  y_val <- validate %>% 
    select(nn[["parameters"]][[x]]) %>% 
    drop_na %>% 
    pull(known_peptide) 
  
  class_counts <- table(y_train)
  
  class_weights <- list(
    "0" = class_counts[2] / class_counts[1],  # Weight for majority class
    "1" = class_counts[2] / class_counts[2]   # Weight for minority class
  )
  rm(model)
  
  early_stopping <- callback_early_stopping(
    monitor = "val_loss",   
    patience = 200,  
    mode = 'min',
    start_from_epoch = 50,
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
    optimizer = optimizer_adam(learning_rate = 0.00001),
    loss = loss_binary_crossentropy(),
    #loss = weighted_binary_crossentropy(weight_0 = class_weights[["0"]], weight_1 = class_weights[["1"]]),
    metrics = c(metric_binary_accuracy(),
                metric_specificity_at_sensitivity(sensitivity = 0.8),
                metric_sensitivity_at_specificity(specificity = 0.8))
  )
  
  
  # Fit the model
  history <- model %>% fit(
    x_train, y_train,
    epochs = 400,
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
  
  
  secretome_aa[[paste0("score_", x)]] <- model %>% 
    predict(secretome_aa %>% 
              select(nn[["parameters"]][[x]]) %>%
              select(-known_peptide) %>%
              as.matrix)
  
}

message(paste0(paste(nn$neural_net, "comp time: ", round(nn$comp_time, 1)), collapse = "\n"))





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
  group_by(accession) %>%
  mutate(across(starts_with("score_nn"), .fns = ~scales::rescale(., to = c(0,1)), .unpack = TRUE)) %>%
  mutate(across(score, .fns = smoother_func, .unpack = TRUE)) %>%
  mutate(across(c("conservation", "relASA", "pathogenicity"), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE))


for(win_size in win_sizes) {
  secretome_aa <- secretome_aa %>%
    group_by(accession) %>%
    mutate(!!paste0("score_b", win_size) := combine_scores(!!sym(paste0("conservation_s", win_size)),
                                                           !!sym(paste0("relASA_s", win_size))))
}

secretome_aa <- secretome_aa %>%
  group_by(accession) %>%
  mutate(across(starts_with("score_nn"), .fns = ~smoother_func(x = ., append_name = "s"), .unpack = TRUE))



secretome_aa <- secretome_aa %>%
  ungroup



saveRDS(secretome_aa, "~/peptide_alg/secretome_aa.rds")



cutoff <- 0.5

roc_res <- Map(\(x) { 
  
  roc(secretome_aa$known_peptide, ifelse(secretome_aa[[x]] > cutoff, 1, 0), ci = TRUE)
  
}, grep("^score_nn", colnames(secretome_aa), value = TRUE))

names(roc_res) <- paste(names(roc_res), "\n[AUC:", round(sapply(roc_res, auc), 4), "]")




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
      predict(nn[["model"]][nn[["neural_net"]] == x][[1]], 
              secretome_aa %>% 
                select(nn[["parameters"]][[x]]) %>% 
                select(-known_peptide) %>%
                as.matrix), ci = TRUE)
  
}, nn$neural_net)

names(roc_res) <- paste(names(roc_res), "\n[AUC:", round(sapply(roc_res, auc), 4), "]")




ggroc(roc_res, linewidth = 1) + ggplot2::theme_bw() + 
  scale_color_discrete(name = "", type = extra_color) +
  ggtitle(paste("ROC Analysis"), subtitle = cutoff) +  
  theme(
    legend.key.height = unit(3, "lines")  
  )



saveRDS(nn, "~/peptide_alg/nn_models_256_64_4_leaky_relu.rds")






