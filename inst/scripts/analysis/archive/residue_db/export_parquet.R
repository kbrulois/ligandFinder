


pq_path <- paste0(s_localDir, "/processed/secretome_parquet3")

dir.create(pq_path)

lists_to_unpack <- c("features", "cons", "dssp", "af_missense", "af_xyz")

secretome <- secretome %>%
  mutate(across(all_of(lists_to_unpack), 
                .fns = list(tagtoremove = ~map(.x, .f = ~`[[`(., "ms")),
                            score = ~map_dbl(.x, .f = ~`[[`(., "score"))),
                .unpack = TRUE)) %>%
  select(-all_of(lists_to_unpack)) %>%
  rename_with(.cols = ends_with("_tagtoremove"), .fn = ~sub("_tagtoremove", "", .))


add_missing_columns <- function(df_list) {
  # Get all unique column names across all data frames
  all_columns <- unique(unlist(lapply(df_list, colnames)))
  
  # Determine correct data types based on existing columns
  column_types <- sapply(all_columns, function(col) {
    existing_types <- unique(unlist(lapply(df_list, function(df) {
      if (col %in% colnames(df)) class(df[[col]]) else NULL
    })))
    
    # Prioritize data types: numeric > character > logical
    if ("numeric" %in% existing_types) return("numeric")
    if ("integer" %in% existing_types) return("integer")
    if ("character" %in% existing_types) return("character")
    if ("factor" %in% existing_types) return("factor")
    return("logical")  # Default if no type found
  })
  
  # Function to create an NA value of the correct type
  create_na <- function(type) {
    switch(type,
           numeric = as.numeric(NA),
           integer = as.integer(NA),
           character = as.character(NA),
           factor = factor(NA, levels = character()),  # Empty factor
           NA)  # Default to logical NA
  }
  
  # Add missing columns to each data frame
  df_list <- lapply(df_list, function(df) {
    missing_cols <- setdiff(all_columns, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- create_na(column_types[col])  # Add NA with correct type
    }
    return(df[, all_columns, drop = FALSE])  # Reorder columns
  })
  
  return(df_list)
}

secretome %>%
  mutate(across(lists_to_unpack, .fns = clean_list_cols)) %>%
  mutate(gene_grp = stringr::str_sub(gene, 1, 1)) %>%
  group_by(gene_grp) %>%
  select(-annotations) %>%
  arrow::write_dataset(path = pq_path, format = "parquet")








residue_db %>%
  mutate(across(lists_to_unpack, .fns = clean_list_cols)) %>%
  mutate(gene_grp = stringr::str_sub(uni_gene, 1, 1)) %>%
  group_by(gene_grp) %>%
  select(-annotations) %>%
  arrow::write_dataset(path = "/oak/stanford/groups/ebutcher/deorphan-AI-ze", format = "parquet")

