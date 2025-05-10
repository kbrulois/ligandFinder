

generate_random_codes <- function(file_path) {

  if(file.exists(file_path)) {
    message("random codes file already exists")
  } else {

  alphanumeric <- c(letters, 0:9)

  ids <- expand.grid(letters,
                     alphanumeric,
                     alphanumeric,
                     alphanumeric,
                     alphanumeric,
                     stringsAsFactors = FALSE)

  ids <- tibble::as_tibble(ids)

  ids <- ids %>%
    dplyr::mutate(id = paste0(Var5, Var4, Var3, Var2, Var1, sep = "")) %>%
    dplyr::mutate(has_numbers = dplyr::if_else(stringr::str_detect(id, "\\d"), "alphanumericeric", "letters_only")) %>%
    dplyr::arrange(desc(has_numbers)) %>%
    dplyr::mutate(usage = factor("unused", levels = c("unused", "used"))) %>%
    dplyr::select(-Var1, -Var2, -Var3, -Var4, -Var5)

  data.table::fwrite(ids, file_path)

  message("random codes saved to:\n", file_path)
}
}

get_codes <- function(n, codes_file = paste0(get_db_path(), "/random_codes.csv")) {

  if (!file.exists(codes_file)) {
    message("Generating random codes file")
    generate_random_codes(file_path = codes_file)
  } else {
    message("using random codes file:\n", codes_file)
  }

  codes_table <- data.table::fread(file = codes_file)

  unused_idx <- which(codes_table$usage == "unused")

  if (length(unused_idx) < n) {
    stop("Not enough unused codes available.")
  }

  selected_idx <- unused_idx[seq_len(n)]
  codes <- codes_table[selected_idx, id]

  codes_table[selected_idx, usage := "used"]

  data.table::fwrite(codes_table, codes_file)

  return(codes)
}

check_random_codes <- function(codes = paste0(get_db_path(), "/random_codes.csv")) {

  if(is.character(codes)) {
    codes_table <- data.table::fread(codes_file)
  } else {
    codes_table <- codes
  }

  code_usage <- table(codes_table[["usage"]])

  message(sapply(names(code_usage)[2:1], \(x) {paste0(x, ": ", round(code_usage[[x]]/1000000, 2), "M\n")}))

}




parse_ranges <- function(y,
                         delim_ranges = "_",
                         delim_start_end = "-") {

  if(delim_ranges == delim_start_end) {
    delim_ranges_old <- delim_ranges
    delim_ranges <- "ZSE76FJ"
    y <- sub(delim_ranges_old, delim_ranges, y)
  }

  y2 <- stringr::str_split(y, delim_ranges, simplify = TRUE)
  colnames(y2) <- paste0("col", 1:ncol(y2))
  y2 <- y2 %>% as_tibble
  logi <- apply(y2, 2, \(z) any(grepl(paste0("\\d+", delim_start_end, "\\d+"), z)))
  to_collapse <- names(logi)[!logi]
  y2 <- y2 %>%
    rowwise %>%
    mutate(id = paste(!!!rlang::syms(to_collapse), sep = "_")) %>%
    select(-all_of(to_collapse)) %>%
    select(id, everything()) %>%
    ungroup
  if(ncol(y2) == 1) {
    y2[["range"]] <- ""
  }
  colnames(y2) <- c("id", "range")
  y2

}


parse_p_id <- function(y) {

  tmp <- tibble(prefix = stringr::str_extract(y, "^[a-z]"),
                cwkov = stringr::str_remove(y, "^[a-z]"),
                suffix = stringr::str_extract(cwkov, "[a-z].+$")
  )

  tmp[[2]] <- stringr::str_remove(tmp[[2]], paste0(tmp[["suffix"]], "$"))

  tmp
}


parse_proteins <- function(file_name,
                           delim_proteins = "_and_",
                           delim_ranges = "_",
                           delim_start_end = "-",
                           protein_names = paste0("p", 1:(1 + max(stringr::str_count(file_name, delim_proteins))))
                           ) {

  tmp <- stringr::str_split(file_name, delim_proteins, simplify = TRUE)
  colnames(tmp) <- protein_names
  tmp %>%
    as_tibble %>%
    mutate(across(everything(), ~parse_ranges(., delim_ranges, delim_start_end), .unpack = TRUE)) %>%
    mutate(across(matches("p\\d_id"), ~parse_p_id(.))) %>%
    select(-all_of(protein_names)) %>%
    tidyr::unnest(matches("p\\d_id"), names_sep = "") %>%
    dplyr::rename_with(.fn = ~stringr::str_remove(., "cwkov$"), .cols = matches("p\\d_id"))
}


split_name <- function(x) {
  tmp <- stringr::str_split(x, "_", simplify = TRUE)
  colnames(tmp) <- c("protein", "annotation")
  tmp %>% as_tibble
}


make_model_names <- function(x, gpcr_mode = TRUE) {

  if(gpcr_mode) {

    model_type <- x %>%
                    filter(protein == "p1" & annotation == "idsuffix") %>%
                    pull(value)

    p1_name <- x %>%
                    filter(protein == "p1" & annotation == "name") %>%
                    pull(value)

    if(length(model_type) == 0) {
      x <- gpcr_list %>%
              filter(uniprot_name == p1_name) %>%
              select(first_AA, last_AA) %>%
              mutate(protein = "p1",
                     annotation = "range",
                     value = paste(first_AA, last_AA, sep = "-")) %>%
              select(protein, annotation, value) %>%
              bind_rows(x, .)
    }



  }

  x <- x %>%
    filter(!annotation %in% c("idprefix", "idsuffix"))

  bind_cols(
    lapply(c("name", "id"), \(z) {
      new_col_name <- paste0("model_", z)
      x %>%
        filter(annotation %in% c(z, "range")) %>%
        group_by(protein) %>%
        summarise(value = paste(value, collapse = ",")) %>%
        summarise(!!new_col_name := paste(value, collapse = ";"))
    })
  )

}


parse_dirname <- function(run_dir = "~/peptide_alg/rename_test",
                          pairing_dir = NULL,
                          ...) {

  id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

  if(!is.null(pairing_dir)) {
    tmp <- tibble(afpd_dir_name = pairing_dir)
  } else {
    tmp <- tibble(afpd_dir_name = fs::dir_ls(run_dir, type = "directory") %>% basename)
  }

  tmp <- tmp %>%

    mutate(parse_proteins(file_name = afpd_dir_name, ...))

  all_ids <- unique(do.call(c, tmp %>% select(ends_with("_id"))))

  all_ids_len <- length(all_ids)

  id_types <- c("Entry Name", "Entry", "Gene Names (primary)")

  id_types <- Map(\(x) {percent <- 100 * sum(all_ids %in% id_map[[x]])/all_ids_len
  message(x, ": ", percent, "%")
  return(percent)}, id_types)

  primary_id <- names(which.max(id_types[names(id_types) != "Gene Names (primary)"]))

  non_primary_id <- names(id_types)[!names(id_types) %in% c("Gene Names (primary)", primary_id)]

  if(id_types[[primary_id]] < id_types[["Gene Names (primary)"]]) {
    stop("Expecting Uniprot ID or Name but detected Gene symbols. Converting to Uniprot Entry Name")
    #tmp <- tmp %>%
    #        mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Gene Names (primary)"]])[.]))
  }

  if(id_types[[non_primary_id]] > 10) {
    warning("The primary identifier detected was '", primary_id, "'\nBut many protein identifiers also mapped to '", non_primary_id, "'")
  }

  converter <- setNames(id_map[[non_primary_id]], id_map[[primary_id]])

  model_type <- c(`Entry` = "id", `Entry Name` = "name")[c(primary_id, non_primary_id)]
  names(model_type) <- c("id", "idnp")

  for(protein in grep("_id$", colnames(tmp), value = TRUE)) {
    tmp <- tmp %>%
      mutate(!!paste0(protein, "np") := converter[!!sym(protein)], .after = protein)
  }


  tmp <- tmp %>%

    tidyr::pivot_longer(matches("p\\d_"), names_to = "name_type") %>%

    tidyr::nest(.key = "parsed_pair", .by = afpd_dir_name) %>%

    mutate(parsed_pair = map(parsed_pair, \(x) {
      x %>%
        mutate(protein = stringr::str_extract(name_type, "^p\\d")) %>%
        mutate(annotation = stringr::str_remove(name_type, "^p\\d_")) %>%
        arrange(protein, annotation) %>%
        filter(value != "") %>%
        mutate(annotation = if_else(annotation %in% c("id", "idnp"), model_type[annotation], annotation)) %>%
        mutate(annotation = factor(annotation, levels = c("idprefix", "name", "id", "idsuffix", "range"))) %>%
        mutate(value = if_else(annotation == "range", sub("\\D+", "-", value), value)) %>%
        select(any_of(c("protein", "annotation", "value")))
    })) %>%
    mutate(map_df(parsed_pair, make_model_names))

  return(tmp)
}

make_new_dirname <- function(input,
                             delim_proteins = "_",
                             delim_ranges = "x",
                             delim_start_end = "x",
                             p1_prefix = "h",
                             p1_suffix = NA,
                             p2_prefix = "h",
                             exclude_p1_range = TRUE) {

  input %>%

    mutate(new_dir_name = map_chr(parsed_pair, \(x) {

      if(exclude_p1_range) {
        x <- x %>%
          filter(!(protein == "p1" & annotation == "range"))
      }

      ps_override <- tibble(p1_idprefix = p1_prefix,
                            p1_idsuffix = p1_suffix,
                            p2_idprefix = p2_prefix) %>%
        tidyr::pivot_longer(everything(), names_to = "name_type", values_to = "override") %>%
        mutate(protein = stringr::str_extract(name_type, "^p\\d")) %>%
        mutate(annotation = stringr::str_remove(name_type, "^p\\d_")) %>%
        select(-name_type)

      full_join(x, ps_override, by = c("protein", "annotation")) %>%
        mutate(annotation = factor(annotation, levels = c("idprefix", "name", "id", "idsuffix", "range"))) %>%
        arrange(protein, annotation) %>%
        mutate(annotation_type = if_else(annotation %in% c("idprefix", "name", "idsuffix"), "name", "range")) %>%
        mutate(value = if_else(is.na(value) & !is.na(override), override, value)) %>%
        filter(annotation != "id" & !is.na(value)) %>%
        mutate(value = if_else(annotation == "range", sub("\\D+", delim_start_end, value), value)) %>%
        group_by(protein, annotation_type) %>%
        summarise(value = paste(value, collapse = ""), .groups = "drop") %>%
        group_by(protein) %>%
        summarise(value = paste(value, collapse = delim_ranges)) %>%
        summarise(value = paste(value, collapse = delim_proteins)) %>%
        pull(value)
    })) %>%

    select(afpd_dir_name, new_dir_name, starts_with("model"), parsed_pair)

}

rename_dir <- function(run_dir = "~/peptide_alg/rename_test",
                       input = tmp,
                       from = "afpd_dir_name",
                       to = "new_dir_name") {

  tmp <- input %>%
    mutate(rename_status = file.rename(paste0(run_dir, "/", !!sym(from)),
                                       paste0(run_dir, "/", !!sym(to))))

  message(sum(tmp[["rename_status"]]), " of ", nrow(tmp), " directories renamed")

  invisible(tmp)

}

safe_fromJSON <- function(txt, encoding = "UTF-8") {

  obj_lines <- readLines(con = txt, encoding = encoding)

  check_json <- jsonlite::validate(txt = obj_lines)
  if(check_json == TRUE ) {
    obj_json <- jsonlite::fromJSON(txt = obj_lines)
  } else {
    obj_json <- jsonlite::fromJSON(txt = "[]") # "[]" represents an empty json.
  }

  obj_json
}

parse_afpd_files <- function(input,
                             dir_name = "afpd_dir_name",
                             run_dir = "~/peptide_alg/rename_test") {


  dat <- input %>%
    mutate(files = map(!!sym(dir_name), ~list.files(paste0(run_dir, "/", .)))) %>%
    mutate(raw_json = map(!!sym(dir_name), \(x) {
      tryCatch({jsonlite::fromJSON(paste(run_dir, x, "ranking_debug.json", sep = "/"))},
               error = function(e) return("problem JSON"))})) %>%
    filter(raw_json != "problem JSON") %>%
    mutate(ranks = map(raw_json, \(x) {
        x %>%
        as_tibble %>%
        select(order) %>%
        tidyr::unnest(order) %>%
        dplyr::rename(model = order) %>%
        mutate(rank = 1:n() - 1) %>%
        arrange(model)
    }))

  dat <- dat %>%
    mutate(files2 = map2(files, ranks, \(x, y) {
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
    mutate(files = pmap(list(files, files2, ranks), \(x, y, z) {
      tmp <- tibble(og_file_name = x,
                    file_mod = y,
                    model = stringr::str_extract(file_mod, "model_\\d+_.*_pred_\\d+"),
                    file_type = stringr::str_remove(file_mod, "model_\\d+_.*_pred_\\d+\\.[^.]+$") %>% stringr::str_remove(., "_$"),
                    file_extension = stringr::str_extract(file_mod, "\\.[^.]+$"),
                    model_num = stringr::str_extract(model, "model_\\d+") %>% stringr::str_remove(., "model_"),
                    pred_num = stringr::str_extract(model, "pred_\\d+") %>% stringr::str_remove(., "pred_"),
                    rank = if_else(!is.na(model), setNames(z[["rank"]], z[["model"]])[model], NA)) %>%
             mutate(rlx = if_else(file_type %in% c("unrelaxed", "relaxed"), file_type, ""))

      rlx_mods <- tmp %>% filter(rlx == "relaxed") %>% pull(model)

      tmp <- tmp %>%
        mutate(rlx = case_when(grepl("ranked", file_type) & model %in% rlx_mods ~ "relaxed",
                               grepl("ranked", file_type) & !model %in% rlx_mods ~ "unrelaxed",
                               !grepl("ranked", file_type) & !file_type %in% c("unrelaxed", "relaxed") ~ NA,
                               TRUE ~ rlx)) %>%
        mutate(model_rlx = paste(model, rlx, sep = "_"))

      ranked_mods <- tmp %>% filter(file_type == "ranked") %>% pull(model_rlx) %>% unique(.)

      tmp %>%
        mutate(rank2 = if_else(model_rlx %in% ranked_mods, paste0("r", as.character(rank)), NA))

    }))

}






make_new_file_names <- function(input,
                                dir_name = "new_dir_name",
                                run_name = "run12",
                                site = "SU",
                                submitter = "KB",
                                algorithm = "AF2v3",
                                random_seed = 42) {

  rc <- get_codes(n = nrow(input))

  input <- bind_cols(input, random_code = rc)

  input %>%
    mutate(files = pmap(list(files, random_code, !!sym(dir_name)), \(x, y, z) {

      uni_mods <- unique_non_na(x %>% arrange(rank) %>% pull(model))
      num_models <- length(uni_mods)
      mc <- setNames(mod_codes[1:num_models], uni_mods)

      ms <- random_seed * num_models
      ms <- seq.int(from = ms, length.out = num_models)
      names(ms) <- unique_non_na(x %>% pull(model))

      x %>%
        mutate(dir_name = z,
               random_code = y,
               model_code = mc[model],
               model_seed = ms[model]) %>%
        mutate(file_type_og = file_type) %>%
        mutate(file_type = file_type_conv[file_type]) %>%
        mutate(rlx = rlx_conv[rlx]) %>%
        mutate(final_code = paste0(site,
                                   submitter,
                                   random_code,
                                   model_code)) %>%
        rowwise %>%
        mutate(new_file_name = if_else(is.na(model),
                                       og_file_name,
                                       stringr::str_flatten(c(dir_name,
                                                              run_name,
                                                              algorithm,
                                                              file_type,
                                                              rlx,
                                                              final_code,
                                                              paste0("s", model_seed,
                                                                     "m", model_num,
                                                                     "p", pred_num),
                                                              rank2),
                                                            collapse = "_",
                                                            na.rm = TRUE)), .after = "og_file_name") %>%
        mutate(new_file_name = if_else(is.na(model),
                                       new_file_name,
                                       paste0(new_file_name, file_extension)))


    }))

}




rename_files <- function(run_dir = "~/peptide_alg/rename_test",
                         input,
                         dir_name = "new_dir_name",
                         from = "og_file_name",
                         to = "new_file_name") {

  tmp <- input %>%
    mutate(rename_status = map2(!!sym(dir_name), files, \(x, y) {
      file.rename(paste(run_dir, x, y[[from]], sep = "/"),
                  paste(run_dir, x, y[[to]], sep = "/"))
      data.table::fwrite(y, file = paste(run_dir, x, "file_name_log.csv", sep = "/"))}))

  invisible(tmp)

}



download_rename_demo <- function(dest_dir = get_db_path()) {

  rd_url <- "https://stacks.stanford.edu/file/sc075gg6264/rename_test.tar.gz"

  db_path <- fs::path_expand(paste0(dest_dir, "/rename_test.tar.gz"))

  downloaded_path <- paste0(dest_dir, "/rename_test")

  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }

  if (!file.exists(downloaded_path)) {
    message("Downloading rename demo data...")
    download.file(rd_url, db_path, mode = "wb")
    untar(db_path, exdir = fs::path_expand(dest_dir))
    fs::file_delete(db_path)
    message("Rename demo data downloaded to: ", downloaded_path)
  } else {
    message("Rename demo data already exists at: ", downloaded_path)
  }

  return(invisible(downloaded_path))
}

condense_model_names <- function(model_names) {

  tibble(model_names = model_names,
         model_num = stringr::str_extract(model_names, "model_\\d+") %>% stringr::str_remove(., "model_"),
         pred_num = stringr::str_extract(model_names, "pred_\\d+") %>% stringr::str_remove(., "pred_"),
         model_names_c = paste0("m", model_num, "p", pred_num)) %>%
    pull(model_names_c)

}




modify_file_names <- function(input_path_models,
                              dir_name,
                              run_name = "deepX14",
                              algorithm = "AF2v3",
                              metrics) {

  rename_data <- data.table::fread(paste(input_path_models, dir_name, "file_name_log.csv", sep = "/")) %>% as_tibble

  if(!"mod_file_name" %in% colnames(rename_data)) {
  rename_data <- left_join(rename_data, metrics %>% rename(model = model_e) %>% select(model, lig1_location, lig1_end), by = "model") %>%
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    rowwise %>%
    mutate(mod_file_name = if_else(is.na(model),
                                   og_file_name,
                                   stringr::str_flatten(c(dir_name,
                                                          run_name,
                                                          algorithm,
                                                          file_type,
                                                          lig1_location,
                                                          rlx,
                                                          final_code,
                                                          paste0("s", model_seed,
                                                                 "m", model_num,
                                                                 "p", pred_num),
                                                          rank2,
                                                          lig1_end),
                                                        collapse = "_",
                                                        na.rm = TRUE)), .after = "new_file_name") %>%
    mutate(mod_file_name = if_else(is.na(model),
                                   mod_file_name,
                                   paste0(mod_file_name, file_extension)))

  data.table::fwrite(rename_data, paste(input_path_models, dir_name, "file_name_log.csv", sep = "/"))

  rename_data <- rename_files(run_dir = input_path_models,
                              input = tibble(new_dir_name = dir_name, files = list(as_tibble(rename_data))),
                              from = "new_file_name",
                              to = "mod_file_name")

  }

}

add_new_file_type <- function(file_type = "spc") {

  file_log <- data.table::fread(paste(input_path_models, dir_name, "file_name_log.csv", sep = "/")) %>%
                  as_tibble

  new_files <- tibble()

  new_file <- file_log %>%
                filter(file_type == "res") %>%
                mutate(file_type = !!sym(file_type))


}


