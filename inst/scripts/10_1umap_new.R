



test2 <- nn_input$C$val %>%
  mutate(data = map(data, \(x) {x %>%
                      mutate(index = row_number()) %>%
                      pivot_wider(names_from = index, values_from = -index)}))

  test2 %>%  pull(data) %>%
  bind_rows -> test



  test <- models[["C"]] %>%
    predict(x_val) %>%
    `[[`("per_index_logits")

  nn_input$C$val <- nn_input$C$val %>%
    mutate(logits = map(1:nrow(test), \(x) {
      tmp <- test[x,,-8] %>% as_tibble
      colnames(tmp) <- names(classes)[-8]
      tmp %>%
        mutate(index = row_number()) %>%
        pivot_wider(names_from = index, values_from = -index)
    }))

xgb_input <- nn_input$C$val$logits %>% bind_rows

umap_config <- umap::umap.defaults

umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 200

umap_res <- umap::umap(d = xgb_input %>% as.matrix, config = umap_config)



nn_input$C$val[["UMAP1"]] <- umap_res[["layout"]][,1]
nn_input$C$val[["UMAP2"]] <- umap_res[["layout"]][,2]

nn_input$C$val <- left_join(nn_input$C$val, nn_input_comb %>% filter(target %in% c("C", "loop_C")) %>% distinct(pep_name, .keep_all = T) %>% select(pep_name, pred, rank, gene, location), by = "pep_name")

data.table::fwrite(nn_input$C$val %>% select(!where(is.list)), "~/AF2_analysis/window_umap_initial2_logits2.csv")


