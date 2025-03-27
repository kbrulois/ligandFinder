

generate_random_codes <- function() {

  alphanumeric <- c(letters, 0:9)

  ids <- expand.grid(letters,
                     alphanumeric,
                     alphanumeric,
                     alphanumeric,
                     alphanumeric,
                     stringsAsFactors = FALSE)

  ids <- dtplyr::lazy_dt(ids)

  ids <- ids %>%
    mutate(id = paste0(Var5, Var4, Var3, Var2, Var1, sep = "")) %>%
    mutate(has_numbers = if_else(str_detect(id, "\\d"), "alphanumericeric", "letters_only")) %>%
    arrange(desc(has_numbers)) %>%
    as_tibble

  data.table::fwrite(ids[["id"]], "~/random_codes.txt")


}


get_codes <- function(n, codes_file = "~/random_codes.txt") {
  codes <- read_lines(codes_file)
  write_lines(codes[-c(1:n)], codes_file)
  return(codes[1:n])
}

# .onLoad <- function() {
#   file_path <- system.file("data/random_codes.txt", package = "ligandFinder")
#   if (file_path == "") {
#     message("random_codes.txt not initialized. Generating now...")
#     generate_random_codes()
#   } else {
#     message()
#   }
# }




parse_ranges <- function(y,
                         delim_ranges = "_",
                         delim_start_end = "-") {

  if(delim_ranges == delim_start_end) {
    delim_ranges_old <- delim_ranges
    delim_ranges <- "ZSE76FJ"
    y <- sub(delim_ranges_old, delim_ranges, y)
  }

  y2 <- stringr::str_split(y, delim_ranges, simplify = TRUE) %>% as_tibble
  logi <- apply(y2, 2, \(z) any(grepl(paste0("\\d+", delim_start_end, "\\d+"), z)))
  to_collapse <- names(logi)[!logi]
  y2 <- y2 %>%
    rowwise %>%
    mutate(id = paste(!!!rlang::syms(to_collapse), sep = "_")) %>%
    select(-all_of(to_collapse)) %>%
    select(id, everything()) %>%
    ungroup
  colnames(y2) <- c("id", "range")[1:ncol(y2)]
  y2

}

parse_proteins <- function(file_name,
                           delim_proteins = "_and_",
                           delim_ranges = "_",
                           delim_start_end = "-",
                           protein_names = paste0("p", 1:(1 + max(stringr::str_count(file_name, delim_proteins)))),
                           p1_range_type = "") {


  tmp <- stringr::str_split(file_name, delim_proteins, simplify = TRUE)
  colnames(tmp) <- protein_names
  tmp %>%
    as_tibble %>%
    mutate(across(everything(), ~parse_ranges(., delim_ranges, delim_start_end), .unpack = TRUE)) %>%
    mutate(p1_range_type = p1_range_type, .before = "p2_id") %>%
    select(-all_of(protein_names))

}



split_name <- function(x) {
  tmp <- stringr::str_split(x, "_", simplify = TRUE)
  colnames(tmp) <- c("protein", "annotation")
  tmp %>% as_tibble
}



make_names <- function(prefix = "h",
                       name,
                       value,
                       custom_suffix = list(p1 = ""), #to be pre-pended to residue ranges
                       include_ranges = paste0("p", 2:length(unique(name))), #which proteins to include ranges
                       type = "file_name", #file_name or model_name (afpd model)
                       delim = list(file_name = list(protein = "_",
                                                     annotation = ""),
                                    model_name = list(protein = ";",
                                                      annotation = ","))
) {

  delim <- delim[[type]]

  if(type == "model_name") {
    prefix <- ""
  }

  names(custom_suffix) <- paste0(names(custom_suffix), "_range")

  tibble(name = name,
         value = value) %>%
    mutate(split_name(name)) %>%
    {left_join(., tibble(name = names(custom_suffix), custom_suffix = unlist(custom_suffix)), by = "name")} %>%
    mutate(custom_suffix = tidyr::replace_na(custom_suffix, "")) %>%
    {if(type == "file_name") {
      mutate(., value = if_else(grepl("_range$", name) & grepl("\\d+-\\d+", value) & protein %in% include_ranges,
                                paste0(custom_suffix, "x", sub("-", "x", value)),
                                if_else(grepl("_range$", name) & grepl("\\d+-\\d+", value) & ! protein %in% include_ranges,
                                        custom_suffix,
                                        value))) } else if(type == "model_name") {.} } %>%
    group_by(protein) %>%
    summarise(value = stringr::str_c(value[value != ""], collapse = delim[["annotation"]])) %>%
    mutate(prefix = prefix) %>%
    rowwise %>%
    mutate(value = stringr::str_c(prefix, value, collapse = "")) %>%
    pull(value) %>%
    paste(., collapse = delim[["protein"]])

}
