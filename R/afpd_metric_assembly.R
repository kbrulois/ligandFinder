




afpd_check_metrics <- function(runs) {

  tmp <- runs %>%
    mutate(parse_proteins(afpd_dir_name, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
    rename(p1_name = p1_id, p2_name = p2_id) %>%
    mutate(file_names = furrr::future_map(afpd_dir, ~fs::dir_ls(.) %>% fs::path_file(.)))


  tmp2 <- furrr::future_map(tmp[["file_names"]], \(x) {
      tibble(has_json_debug = "ranking_debug.json" %in% x,
             has_metrics_csv = "metrics_v2.csv" %in% x,
             has_contact_rds = "metrics_v2c.rds" %in% x,
             num_models = sum(grepl("_xtr_.*.pdb$", x)))
    }) %>% bind_rows()


  tmp <- bind_cols(tmp, tmp2)

  tmp$metrics <- furrr::future_map(tmp[["afpd_dir"]], \(x) {
      file <- fs::path(x, "metrics_v2.csv")
      if(file.exists(file)) {
        return(data.table::fread(file) %>% as_tibble)
      } else {
        return(NULL)
      }
    })

  metrics_good <- tmp %>%
    filter(has_metrics_csv) %>%
    filter(!map_lgl(metrics, is.null)) %>%
    filter(map_int(metrics, ~nrow(.)) == num_models) %>%
    filter(map_lgl(metrics, ~"depth" %in% colnames(.))) %>%
    pull(afpd_dir_name)


  tmp %>%
    mutate(metrics_good = afpd_dir_name %in% !!metrics_good) %>%
    mutate(num_xtr = map_int(file_names, ~stringr::str_detect(., "_xtr_") %>% sum)) %>%
    mutate(num_ark = map_int(file_names, ~stringr::str_detect(., "_ark_") %>% sum))

}

afpd_check_contacts <- function(runs) {

  tmp <- runs %>%
    mutate(contacts = furrr::future_map(afpd_dir, \(x) {
      file <- fs::path(x, "metrics_v2c.rds")
      if(file.exists(file)) {
        return(tryCatch({readRDS(file)}, error = function(e) NULL))
      } else {
        return(NULL)
      }
    }))


  contacts_good <- tmp %>%
    filter(has_contact_rds & has_metrics_csv) %>%
    filter(map2_lgl(metrics, contacts, \(x, y) {
      if(is.null(y)) {test <- rep("irrelevant", nrow(x))} else {
        test <- map_chr(y, \(z) {
          if(is.null(z)) {return("irrelevant")} else
            if(nrow(z) == 0) {return("irrelevant")} else
              if(!"area" %in% colnames(z)) {return("irrelevant")} else
              {return("relevant")}
        })
      }
      return(identical(x[["location"]], test))
    })) %>%
    pull(afpd_dir_name)

  tmp %>%
    mutate(contacts_good = afpd_dir_name %in% !!contacts_good)
}




archive_metrics <- function(runs) {

  num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) %/% 1.5

  future::plan(strategy = future::multicore(workers = num_cores))

  runs %>%
    filter(has_metrics_csv) %>%
    mutate(tar_cmd = paste("tar -rf",
                           fs::path(oak_models, run_name, afpd_dir_name, ext = "tar"),
                           "-C",
                           afpd_dir,
                           "metrics_v2.csv")) %>%
    mutate(tar_cmd = if_else(has_contact_rds, paste(tar_cmd, "metrics_v2c.rds"), tar_cmd)) %>%
    pull(tar_cmd) %>%

  furrr::future_walk(., system)

}




