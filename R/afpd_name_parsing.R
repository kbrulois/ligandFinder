

generate_random_codes <- function() {

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

  file_path <- paste0(system.file("extdata", package = "ligandFinder"), "/random_codes.csv")

  data.table::fwrite(ids, file_path)

  message("random codes file saved:\n", file_path)

}

get_codes <- function(n, codes_file = system.file("extdata/random_codes.csv", package = "ligandFinder")) {

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

check_random_codes <- function(codes = system.file("extdata/random_codes.csv", package = "ligandFinder")) {

  if(is.character(codes)) {
    codes_table <- data.table::fread(codes_file)
  } else {
    codes_table <- codes
  }

  code_usage <- table(codes_table[["usage"]])

  message(sapply(names(code_usage)[2:1], \(x) {paste0(x, ": ", round(code_usage[[x]]/1000000, 2), "M\n")}))

}



.onAttach <- function(libname = .libPaths(), pkgname = "ligandFinder") {
  file_path <- system.file("extdata/random_codes.csv", package = "ligandFinder")
  if (file_path == "") {
    message("Generating random codes file")
    generate_random_codes()
  } else {
    message("random codes file:\n", file_path)
  }
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


make_model_names <- function(x) {
  x <- x %>%
    filter(!annotation %in% c("idprefix", "idsuffix"))

  bind_cols(
    lapply(unname(model_type), \(z) {
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
                          ...) {

  id_map <- readRDS(system.file("data/id_mapping.rds", package = "ligandFinder"))

  tmp <- tibble(afpd_dir_name = list.files(run_dir)) %>%

    mutate(parse_proteins(afpd_dir_name,
                          delim_proteins = "_",
                          delim_ranges = "x",
                          delim_start_end = "x"))

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

  tmp <- tmp <- input %>%
    mutate(rename_status = file.rename(paste0(run_dir, "/", !!sym(from)),
                                       paste0(run_dir, "/", !!sym(to))))

  message(sum(tmp[["rename_status"]]), " of ", nrow(tmp), " directories renamed")

  invisible(tmp)

}





