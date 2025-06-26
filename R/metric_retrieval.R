get_metrics <- function(run_dir,
                        file = "metrics_v2.csv",
                        reader = data.table::fread) {



  jobs <- tibble(dir_name = fs::dir_ls(run_dir, type = "directory") %>% basename) %>%
    mutate(complete = furrr::future_map_lgl(dir_name, \(x) {
      file.exists(paste(run_dir, x, file, sep = "/"))
    })) %>%
    filter(complete) %>%
    mutate(group = paste0("job", ntile(n = num_of_grps)))

  res <- furrr::future_map(unique(jobs[["group"]]), \(job) {
    dirs <- jobs %>%
      filter(group == job) %>%
      pull(dir_name)


    metrics <- map(dirs,
                   ~reader(paste0(run_dir, "/", ., "/", file)))

    metrics <- metrics[sapply(metrics, \(x) nrow(x) != 0)]

    bind_rows(metrics)

  }
  )

  bind_rows(res) %>% as_tibble

}

get_contacts <- function(run_dir,
                         metric_file = "metrics_v2.csv",
                         contact_file = "metrics_v2c.rds",
                         reader = data.table::fread) {



  jobs <- tibble(dir_name = fs::dir_ls(run_dir, type = "directory") %>% basename) %>%
    mutate(complete = furrr::future_map_lgl(dir_name, \(x) {
      file.exists(paste(run_dir, x, file, sep = "/"))
    })) %>%
    filter(complete) %>%
    mutate(group = paste0("job", ntile(n = num_of_grps)))

  res <- furrr::future_map(unique(jobs[["group"]]), \(job) {
    dirs <- jobs %>%
      filter(group == job) %>%
      pull(dir_name)


    metrics <- map(dirs,
                   ~reader(paste0(run_dir, "/", ., "/", file)))

    metrics <- metrics[sapply(metrics, \(x) nrow(x) != 0)]

    bind_rows(metrics)

  }
  )

  bind_rows(res) %>% as_tibble

}

get_known_pairs <- function() {


  known_pairs <- list(c("GPR25", "CXL17"),
                      c("CCR9", "CCL25"),
                      c("GPR15", "GP15L"),
                      c("CML1", "RARR2"),
                      c("CML2", "RARR2"),
                      c("CCRL2", "RARR2"))

  ligand_list <- readRDS(system.file("extdata/ligand_list.rds", package = "ligandFinder"))

  known_pairs2 <- ligand_list %>%
    rowwise %>%
    mutate(known = map2(.x = uniprot_name,
                        .y = receptor,
                        .f = \(x, y) {c(x, y)})) %>%
    ungroup %>%
    mutate(known_lgl = map_lgl(known, ~any(is.na(.)))) %>%
    filter(!known_lgl) %>%
    pull(known) %>%
    unname


  id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

  chemokine_pairs <- readRDS(system.file("data/Chemokine_Pairs.rds", package = "ligandFinder")) %>%
    slice(-1) %>%
    as_tibble

  chem_pairs <- chemokine_pairs %>%
    rowwise %>%
    mutate(across(all_of(c("V1", "V2")), \(x) {

      if(x %in% id_map$`Gene Names (primary)`) {
        return(setNames(id_map$`Entry Name`, id_map$`Gene Names (primary)`)[x])
      } else {
        id_map %>%
          mutate(thing = map_lgl(`Gene Names`, \(y) {
            stringr::str_detect(y, x)
          })) %>%
          pull(thing) -> gene_alts
        if(sum(gene_alts) == 1) {
          return(id_map$`Entry Name`[gene_alts])
        } else {
          return(NA)
        }

      }})) %>%
    rowwise %>%
    mutate(known = map2(.x = V1,
                        .y = V2,
                        .f = \(x, y) {c(x, y)})) %>%
    ungroup %>%
    mutate(known_lgl = map_lgl(known, ~any(is.na(.)))) %>%
    filter(!known_lgl) %>%
    pull(known) %>%
    unname

  c(known_pairs, known_pairs2, chem_pairs)

}
