
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

get_unique_or_collapse <- function(x, y) {
  uni_res <- unique(c(x,y))
  uni_res <- uni_res[!is.na(uni_res)]
  if(length(uni_res) == 1) {
    return(uni_res)
  } else if(length(uni_res) == 0) {
    return(NA)
  } else {
    return(paste(c(x,y), collapse = " "))
  }
}

get_unique_or_collapse_v <- Vectorize(get_unique_or_collapse)

get_unique_or_collapse <- function(x, delim = "; ") {
  if(length(unique(x)) == 1) {
    return(x[1])
  } else {
    return(paste(unique(x), collapse = delim))
  }
}


getBallesterosFromGPCRDB <- function(gene_name,
                                     gpcrdb_url = "https://gpcrdb.org/services/residues/extended/") {

  url_link <- paste0(gpcrdb_url, gene_name)

  dat <- jsonlite::read_json(url_link)

  tibble(AA = sapply(dat, `[[`, "amino_acid"),
         index = sapply(dat, `[[`, "sequence_number"),
         protein_segment = sapply(dat, `[[`, "protein_segment"),
         display_generic_number = sapply(dat, `[[`, "display_generic_number"),
         BW = sapply(dat, \(x) {tryCatch(x[["alternative_generic_numbers"]][[1]][["label"]],
                                         error = function(e) {
                                           message("Index out of bounds! Returning NA.")
                                           return(NA)
                                         })}))
}


get_lengths <- function(x, params = c("N-term", "C-term", paste0("ECL", 1:3), paste0("ICL", 1:3), paste0("TM", 1:7))) {
  dummy_table <- tibble(!!!setNames(rep(list(NA), length(params)), paste("bw: length", params)))
  if(nrow(x) == 0) {
    return(dummy_table)

  } else {
    tmp <- table(x[["protein_segment"]])
    names(tmp) <- paste("bw: length", names(tmp))
    tmp <- bind_rows(dummy_table, tmp) %>% filter(!if_all(everything(), is.na)) %>% select(all_of(paste("bw: length", params)))
    return(tmp)
  }
}

get_bw_inds <- function(x, params = c("1.50", "7.50")) {
  dummy_table <- tibble(!!!setNames(rep(list(NA), length(params)), paste("bw: indices", params)))
  if(nrow(x) == 0) {
    return(dummy_table)
  } else {
    tmp <- Map(\(y) x[["index"]][x[["BW"]] %in% y], params)
    names(tmp) <- paste("bw: indices", names(tmp))
    tmp <- bind_rows(dummy_table, tmp) %>% filter(!if_all(everything(), is.na))
    if(nrow(tmp) == 1) {
      return(tmp)
    } else {
      return(dummy_table)
    }
  }
}


get_sp <- function(x) {

  sp <- x %>%
    filter(type == "signal peptide" & source == "uniprot")
  if(nrow(sp) == 0) {
    return(as.integer(NA))
  } else {
    return(as.integer(sp[["end"]] + 1))
  }

}

get_lastAA <- function(x) {

  sp <- x %>%
    filter(type == "c-terminus" & source == "manual")
  if(nrow(sp) == 0) {
    return(as.integer(NA))
  } else {
    return(as.integer(sp[["end"]]))
  }

}

clean_feats <- function(x) {
  if(is.null(x)) {
    return(tibble(!!!setNames(rep(list(NA), 6), c("type", "evidence", "start", "end", "source", "priority"))))
  } else {
    return(x)
  }
}

extend_bw_notation <- \(x) {

  n_term_residues <- x[["protein_segment"]] == "N-term"
  x[["BW"]][n_term_residues] <- paste0("N", sum(n_term_residues):1)

  c_term_residues <- x[["protein_segment"]] == "C-term"
  x[["BW"]][c_term_residues] <- paste0("C", 1:sum(c_term_residues))

  x

}


summarize_bw <- function(gpcr_list = system.file("extdata/gpcr_list.rds", package = "ligandFinder")) {

  gpcr_list <- readRDS(gpcr_list)

  bw_align <- map(gpcr_list$`bw: full_table`, \(x) {tmp <- unique(x[["BW"]]); return(tmp[!is.na(tmp)])})
  bw_align <- table(unlist(bw_align))
  bw_align <- tibble(BW = names(bw_align),
                     prop = as.integer(unname(bw_align)))
  bw_align <- bw_align %>%
    mutate(CP = if_else(BW %in% cp_bw, "CP", "")) %>%
    mutate(prop = round(100 * prop/805, 0)) %>%
    mutate(name = paste(BW, prop, CP, sep = "_")) %>%
    mutate(name = stringr::str_remove(name, "_$")) %>%
    slice(gtools::mixedorder(name))

  bw_align %>% filter(prop > 20)

}


