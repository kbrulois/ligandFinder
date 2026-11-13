
get_nonNA_ranges <- function(x, column = "AA") {

  tmp <- x %>%
    mutate(id = 1:n()) %>%
    mutate(non_na_group = cumsum(!is.na(!!sym(column)) & lag(is.na(!!sym(column)), default = TRUE))) %>%
    filter(!is.na(!!sym(column))) %>%
    group_by(non_na_group) %>%
    summarize(start = min(id), end = max(id)) %>%
    select(-non_na_group)

  if(nrow(tmp) > 1) {
    tmp <- tibble(start = NA, end = NA)
  }
  tmp

}


