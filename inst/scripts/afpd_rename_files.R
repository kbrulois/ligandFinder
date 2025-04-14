



parse_afpd_files <- function(input,
                             dir_name = "afpd_dir_name",
                             run_dir = "~/peptide_alg/rename_test") {

  dat <- input %>%
    mutate(files = map(!!sym(dir_name), ~list.files(paste0(run_dir, "/", .)))) %>%
    mutate(ranks = map(!!sym(dir_name), \(x) {tryCatch({
      jsonlite::read_json(paste(run_dir, x, "ranking_debug.json", sep = "/")) %>%
        as_tibble %>%
        select(order) %>%
        tidyr::unnest(order) %>%
        dplyr::rename(model = order) %>%
        mutate(rank = 1:n() - 1) %>%
        arrange(model)
    })}))

  dat <- dat %>%
    mutate(files = map2(files, ranks, \(x, y) {
      rank_mapping <- setNames(paste0("ranked_", y[["model"]]),
                               paste0("ranked_", y[["rank"]]))
      ranked_files <- grep("ranked_\\d+", x, value = TRUE)
      for(i in ranked_files) {
        rank <- stringr::str_extract(x[x == i], "ranked_\\d+")
        x[x == i] <- stringr::str_replace(x[x == i], pattern = rank, replacement = rank_mapping[rank])
      }
      x
    }))

  dat %>%
    mutate(files = map2(files, ranks, \(x, y) {
      tibble(file_name = x,
             model = stringr::str_extract(file_name, "model_\\d+_.*_pred_\\d+"),
             file_type = stringr::str_remove(file_name, "model_\\d+_.*_pred_\\d+\\.[^.]+$") %>% stringr::str_remove(., "_$"),
             file_extension = stringr::str_extract(file_name, "\\.[^.]+$"),
             model_num = stringr::str_extract(model, "model_\\d+") %>% stringr::str_remove(., "model_"),
             pred_num = stringr::str_extract(model, "pred_\\d+") %>% stringr::str_remove(., "pred_"),
             rank = if_else(!is.na(model), setNames(y[["rank"]], y[["model"]])[model], NA)) %>%
      mutate(file_type = if_else(file_type == "ranked", paste0(file_type, file_extension), file_type))
    }))

}



make_new_file_names <- function(input,
                                random_seed = 42,
                                algorithm = "AF2mV3",#options "AF3", "boltz?"
                                site = "SU",#options "AP", "SD"
                                submitter = "KB") {

  test <- table(input$files[[1]]$file_type)

  rc <- get_codes(n = nrow(input))

  input %>%
    mutate(files = map(files, \(x) {
        uni_mods <- unique(x[["model"]])


    }))

}

tmp <- parse_dirname()


tmp <- make_new_dirname(tmp)

tmp2 <- parse_afpd_files(input = tmp,
                         dir_name = "afpd_dir_name",
                         run_dir = "~/peptide_alg/rename_test")








rename_dir(run_dir = "~/peptide_alg/rename_test",
           input = tmp,
           from = "afpd_dir_name",
           to = "new_dir_name")


rename_dir(run_dir = "~/peptide_alg/rename_test",
           input = tmp,
           to = "afpd_dir_name",
           from = "new_dir_name")






tmp2 <- parse_dirname(delim_proteins = "_",
                      delim_ranges = "x",
                      delim_start_end = "x")

tmp3 <- make_new_dirname(input = tmp2)






