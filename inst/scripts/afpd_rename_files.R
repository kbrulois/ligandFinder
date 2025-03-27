

#mutate(p2_id = paste(p2_id, p3_id, sep = "_")) %>%
#mutate(p3_id = paste(p4_id, p5_id, sep = "_")) %>%
#select(-p4_id, -p5_id) %>%
#mutate(across(ends_with("_id"), ~gsub("[^A-Z0-9]", "", .))) %>%
#mutate(across(ends_with("_id"), ~setNames(id_map[["Entry Name"]], id_map[["Gene Names (primary)"]])[.])) %>%
#mutate(new_file_name = paste(p1_id, paste(p2_id, p2_range, sep = "x"), sep = "_"))







# ids <- ids %>%
#   mutate(has_numbers = if_else(str_detect(id, "\\d"), "alphanumericeric", "letters_only")) %>%
#   arrange(desc(has_numbers))

chars <- c(letters, 0:9)

ids <- expand.grid(letters, letters, letters, letters, letters, stringsAsFactors = FALSE)

library(dtplyr)

ids <- dtplyr::lazy_dt(ids)

ids <- ids %>%
  mutate(id = paste0(Var5, Var4, Var3, Var2, Var1, sep = "")) %>%
  as_tibble

write_lines(ids$id, "~/random_codes.txt")

get_codes <- function(n, codes_file = "~/random_codes.txt") {
  codes <- read_lines(codes_file)
  write_lines(codes[-c(1:n)], codes_file)
  return(codes[1:n])
}

get_codes(10)

get_codes(10)








