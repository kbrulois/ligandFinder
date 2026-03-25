



#install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))

library(xgboost)



start <- Sys.time()

    x_train <- models$C$ %>%
      mutate(data = map(data, \(x) {x %>%
          mutate(index = row_number()) %>%
          pivot_wider(names_from = index, values_from = -index)})) %>%
      pull(data) %>%
      bind_rows %>%
      as.matrix

    x_train <- xgb_input %>% as.matrix

    y_train <- nn_input$C$val$known

    x_val <- secretome_aa %>%
      filter(gene %in% c(known_genes[[tag]][!known_samp], ctrl_genes_sub[!ctrl_samp])) %>%
      select(nn[["parameters"]][[x]]) %>%
      drop_na %>%
      select(-known) %>%
      as.matrix

    y_val <- secretome_aa %>%
      filter(gene %in% c(known_genes[[tag]][!known_samp], ctrl_genes_sub[!ctrl_samp])) %>%
      select(nn[["parameters"]][[x]]) %>%
      drop_na %>%
      mutate(known = as.integer(known)) %>%
      pull(known)

    dtrain <- xgb.DMatrix(data = x_train, label = y_train)
    dvalid <- xgb.DMatrix(data = x_val,  label = y_val)

    watchlist <- list(train = dtrain,
                      eval = dvalid)
    watchlist <- list(train = dtrain)

    bst <- xgb.train(
      params = list(
        objective = "binary:logistic",
        eval_metric = "auc",
        eta = 0.03,
        subsample = 0.7,
        colsample_bytree = 0.7,
        min_child_weight = 5,
        gamma = 1,
        lambda = 1,
        alpha = 0.5,
        tree_method = "hist",
        max_depth = 6,
        eval_metric = "auc"
      ),
      data = dtrain,
      nrounds = 500,
      watchlist = watchlist,
      early_stopping_rounds = 100,
      verbose = 1
    )

    nn_input$C$val$xgboost_C <- bst %>%
      predict(xgb_input %>%
                as.matrix)

  }
}


View(nn_input$C$val)



shap_res <- xgb.importance(model = bst)

shap_res <- shap_res %>% as_tibble %>%
  mutate(index = stringr::str_extract(Feature, "_\\d+$") %>% stringr::str_remove(., "^_") %>% as.integer) %>%
  mutate(Feature = stringr::str_remove(Feature, "_\\d+$"))



shap_res %>%data.table::fwrite(., "~/AF2_analysis/xgboost_pep_logits.csv")




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







