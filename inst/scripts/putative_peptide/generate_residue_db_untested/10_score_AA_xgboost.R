



#install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))

library(xgboost)



for(x in nn[["neural_net"]]) {
  for(tag in names(known_genes)) {


  message("computing ", x)

  start <- Sys.time()

  x_train <- secretome_aa %>%
    filter(gene %in% c(known_genes[[tag]][known_samp], ctrl_genes_sub[ctrl_samp])) %>%
    select(nn[["parameters"]][[x]]) %>%
    drop_na %>%
    select(-known) %>%
    as.matrix

  y_train <- secretome_aa %>%
    filter(gene %in% c(known_genes[[tag]][known_samp], ctrl_genes_sub[ctrl_samp])) %>%
    select(nn[["parameters"]][[x]]) %>%
    drop_na %>%
    mutate(known = as.integer(known)) %>%
    pull(known)

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

bst <- xgb.train(
  params = list(
    objective = "binary:logistic",
    eta = 0.08,
    subsample = 0.5,
    min_child_weight = 0,
    gamma = 0,
    lambda = 0,
    alpha = 1,
    tree_method = "exact",
    max_depth = 6,
    eval_metric = "auc"
  ),
  data = dtrain,
  nrounds = 500,
  watchlist = watchlist,
  early_stopping_rounds = 100,
  verbose = 1
)

secretome_aa[[paste0(tag, "_", sub("^nn", "xgb", x))]] <- bst %>%
  predict(secretome_aa %>%
            select(nn[["parameters"]][[x]]) %>%
            select(-known) %>%
            as.matrix)

  }
}






xgb.importance(model = bst) %>% data.table::fwrite(., "~/AF2_analysis/xgboost_pep2.csv")









