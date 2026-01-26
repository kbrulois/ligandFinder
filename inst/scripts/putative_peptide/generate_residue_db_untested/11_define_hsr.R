



###add KK/KR and other features

secretome_aa <- readRDS(paste0(s_localDir, "/processed/secretome_aa.rds"))

secretome_aa %>%
  select(!matches("_lead|_lag")) %>% colnames()




secretome <- left_join(secretome,
                       secretome_aa %>%
                         select(!matches("_lead|_lag")) %>%
                         group_by(accession) %>%
                         nest(.key = "aa_scores"),
                       by = "accession")



gc()
gc()
start <- Sys.time()
future::plan(strategy = future::sequential())

secretome[["aa_scores"]] <- furrr::future_map2(secretome[["aa_scores"]], secretome[["sequence_uni"]], \(x, y) {

  if(is.null(x)) {return(x)} else {
  tmp <- tibble(index = 1:nchar(y),
                AA = str_split(y, "", simplify = TRUE) %>% `c`)

  tmp <- left_join(tmp, x, by = c("index", "AA"))

  return(tmp)
  }

})

end <- Sys.time()
end - start


features <- secretome %>% filter(accession == "P01275") %>% pull(features) %>% `[[`(1)

aa_scores <- secretome %>% filter(accession == "P01275") %>% pull(aa_scores) %>% `[[`(1)


sequence_uni <- secretome %>% filter(accession == "P01275") %>% pull(sequence_uni)



add_hsr <- function(features, aa_scores, score = "pep_xgb3c_s6", cutoff = 0.6) {

  x <- aa_scores[[score]]
  x[is.nan(x)] <- 0
  #x <- x - min(x, na.rm = TRUE)
  #x <- x/max(x, na.rm = TRUE)

  vec <- which(x > cutoff)

  if(length(vec) > 0) {

    breaks <- c(0, which(diff(vec) != 1), length(vec))

    sequences <- lapply(seq_along(breaks[-1]), \(i) vec[(breaks[i] + 1):breaks[i + 1]])

    hsr_dat <- bind_rows(lapply(sequences, \(y) tibble(type = paste0("hsr_", min(y, na.rm = TRUE), "-", max(y, na.rm = TRUE)),
                                                       evidence = paste0(score, "_", cutoff),
                                                       start = min(y, na.rm = TRUE),
                                                       end = max(y, na.rm = TRUE),
                                                       source = "hsr")
    ))
  } else {
    hsr_dat <- tibble(type = character(),
                      evidence = character(),
                      start = integer(),
                      end = integer(),
                      source = character())
  }

  bind_rows(features, hsr_dat)
}

combine_phs <- function(x) {

  comb_dat <- x %>%
    filter(source %in% c("phs_n-term-Cys",
                         "phs_c-term-Cys",
                         "phs_n-term",
                         "phs_c-term",
                         "phs_dbr")) %>%
    arrange(start) %>%
    select(start, end) %>%
    drop_na() %>%
    {IRanges::IRanges(start = .[["start"]], end = .[["end"]])} %>%
    IRanges::reduce(., min.gapwidth = 0) %>%
    as.data.frame() %>%
    as_tibble()

  if(nrow(comb_dat) > 0) {
    x <- comb_dat %>%
      mutate(type = paste0("phs_", start, "-", end),
             evidence = "phs", .before = everything()) %>%
      select(-width) %>%
      mutate(source = "phs") %>%
      plyr::rbind.fill(x, .)
  }
  return(x)
}

x <- secretome %>%
  filter(gene == "SCGB3A2") %>%
  pull(features) %>%
  `[[`(1)

topo <- secretome %>%
  slice(132) %>%
  pull(topo) %>%
  `[[`(1)

topo <- secretome %>%
          filter(gene == "CCR10") %>%
            pull(topo) %>%
            `[[`(1)


limit_phs <- function(x,
                      topo,
                      offlimits = c("signal peptide", "E", "Dibasic"),
                      offlimits_alt = c("signal peptide", "strand", "Dibasic")) {

  if(!"alpha fold" %in% x[["source"]]) {
    offlimits <- offlimits_alt
  }

  comb_dat <- x %>%
    filter(type %in% offlimits | source == "phs")



  if(!is.na(topo)) {

      topo <- stringr::str_split(topo, "", simplify = TRUE) %>% `c`

      if(!all(topo == topo[1])) {

      tmp <- factor_to_uniprotFeature(tibble(topo = topo),
                                        source = "uni_topo", evidence = "uni_topo")

      tmp <- tmp %>%
        filter(!type %in% c("e", "-"))

      comb_dat <- bind_rows(comb_dat, tmp)
      }

  }



  cond <- sum(comb_dat[["type"]] %in% offlimits) > 0 &
    sum(comb_dat[["source"]] == "phs") > 0

  if(cond) {
    comb_dat <- comb_dat %>%
      arrange(start) %>%
      select(start, end, type, source) %>%
      drop_na() %>%
      mutate(splitter = ifelse(source == "phs", "phs", "ol")) %>%
      split(., .[["splitter"]]) %>%
      {cross_join(.[["phs"]], .[["ol"]], suffix = c("_phs", "_ol"))} %>%
      {query <- IRanges::IRanges(start = .[["start_phs"]], end = .[["end_phs"]], name = .[["type_phs"]])
      subject <- IRanges::IRanges(start = .[["start_ol"]], end = .[["end_ol"]], name = .[["type_ol"]])
      IRanges::setdiff(query, subject)} %>%
      as.data.frame() %>%
      as_tibble()
  } else {
    comb_dat <- tibble(start = integer())
  }

  if(nrow(comb_dat) > 0) {
    x <- comb_dat %>%
      mutate(type = paste0("phs_", start, "-", end),
             evidence = "phs", .before = everything()) %>%
      select(-width) %>%
      mutate(source = "phs") %>%
      plyr::rbind.fill(x %>% filter(source != "phs"), .)
  }

  return(x)
}

combine_phs_hsr <- function(x) {

  comb_dat <- x %>%
    filter(source %in% c("phs", "hsr")) %>%
    arrange(start) %>%
    select(start, end) %>%
    drop_na() %>%
    {IRanges::IRanges(start = .[["start"]], end = .[["end"]])} %>%
    IRanges::reduce(., min.gapwidth = 0) %>%
    as.data.frame() %>%
    as_tibble()

  if(nrow(comb_dat) > 0) {
    x <- comb_dat %>%
      mutate(type = paste0("phs_hsr_", start, "-", end),
             evidence = "phs_hsr", .before = everything()) %>%
      select(-width) %>%
      mutate(source = "phs_hsr") %>%
      plyr::rbind.fill(x, .)
  }
  return(x)
}

combine_a_b <- function(x, a, b, tag) {

  comb_dat <- x %>%
    filter(source %in% c(a, b)) %>%
    arrange(start) %>%
    select(start, end) %>%
    drop_na() %>%
    {IRanges::IRanges(start = .[["start"]], end = .[["end"]])} %>%
    IRanges::reduce(., min.gapwidth = 0) %>%
    as.data.frame() %>%
    as_tibble()

  if(nrow(comb_dat) > 0) {
    x <- comb_dat %>%
      mutate(type = paste0(tag, "_", start, "-", end),
             evidence = tag, .before = everything()) %>%
      select(-width) %>%
      mutate(source = tag) %>%
      plyr::rbind.fill(x, .)
  }
  return(x)
}

define_phs_hsr_gtp_ol <- function(x) {

  if("gpcrdb_gtp" %in% x[["source"]] & "phs_hsr" %in% x[["source"]]) {

    ol_dat <- x %>%
      filter(source %in% c("gpcrdb_gtp", "phs_hsr")) %>%
      distinct(start, end, .keep_all = TRUE) %>%
      arrange(start) %>%
      select(type, start, end, source) %>%
      drop_na() %>%
      split(., .[["source"]]) %>%
      {cross_join(.[["gpcrdb_gtp"]], .[["phs_hsr"]], suffix = c("_gpcrdb_gtp", "_phs_hsr"))} %>%
      mutate(overlap_start = pmax(start_gtp, start_phs_hsr),
             overlap_end = pmin(end_gtp, end_phs_hsr)) %>%
      filter(overlap_start <= overlap_end) %>%
      mutate(overlap_length = overlap_end - overlap_start + 1,
             percent_ol_gtp = 100 * overlap_length/(end_gtp - start_gtp + 1),
             percent_ol_phs_hsr = 100 * overlap_length/(end_phs_hsr - start_phs_hsr + 1)) %>%
      filter(percent_ol_gtp == max(percent_ol_gtp), .by = type_gtp) %>%
      select(type_gtp, type_phs_hsr, percent_ol_gtp, percent_ol_phs_hsr)

    if(nrow(ol_dat) > 0) {

      bind_rows(
        ol_dat %>% dplyr::rename(type = type_gtp, overlap_region = type_phs_hsr),
        ol_dat %>% dplyr::rename(type = type_phs_hsr, overlap_region = type_gtp)
      ) %>%
        filter((percent_ol_gtp + percent_ol_phs_hsr) == max(percent_ol_gtp + percent_ol_phs_hsr), .by = type) %>%
        left_join(x, ., by = "type") -> x

    }
  }

  return(x)

}

define_ol <- function(x, a, b) {

  start_a <- paste0("start_", a)
  end_a <- paste0("end_", a)
  start_b <- paste0("start_", b)
  end_b <- paste0("end_", b)

  if(b == "uni_pep") {
    x <- x %>%
      mutate(source = if_else(type == "peptide" & source == "uniprot", "uni_pep", source))
  }

  if(a %in% x[["source"]] & b %in% x[["source"]]) {

    ol_dat <- x %>%
      filter(source %in% c(a, b)) %>%
      group_by(source) %>%
      distinct(start, end, .keep_all = TRUE) %>%
      ungroup %>%
      arrange(start) %>%
      select(type, start, end, source) %>%
      drop_na() %>%
      split(., .[["source"]]) %>%
      {cross_join(.[[a]], .[[b]], suffix = c(paste0("_", a), paste0("_", b)))} %>%
      mutate(overlap_start = pmax(!!sym(start_a), !!sym(start_b)),
             overlap_end = pmin(!!sym(end_a), !!sym(end_b))) %>%
      filter(overlap_start <= overlap_end) %>%
      mutate(!!paste0("overlap_length_", a, ":", b) := overlap_end - overlap_start + 1,
             !!paste0("percent_ol_", a, ":", b) := 100 * !!sym(paste0("overlap_length_", a, ":", b))/(!!sym(end_a) - !!sym(start_a) + 1),
             !!paste0("percent_ol_", b, ":", a) := 100 * !!sym(paste0("overlap_length_", a, ":", b))/(!!sym(end_b) - !!sym(start_b) + 1)) %>%
      filter(!!sym(paste0("percent_ol_", a, ":", b)) == max(!!sym(paste0("percent_ol_", a, ":", b))), .by = !!sym(paste0("type_", a))) %>%
      select(starts_with("type_"), starts_with("percent_ol_"))

    if(nrow(ol_dat) > 0) {

      bind_rows(
        ol_dat %>% dplyr::rename(type = !!sym(paste0("type_", a)), !!sym(paste0("overlap_region_", b, ":", a)) := !!sym(paste0("type_", b))),
        ol_dat %>% dplyr::rename(type = !!sym(paste0("type_", b)), !!sym(paste0("overlap_region_", a, ":", b)) := !!sym(paste0("type_", a)))
      ) %>%
        filter((!!sym(paste0("percent_ol_", a, ":", b)) + !!sym(paste0("percent_ol_", b, ":", a))) == max(!!sym(paste0("percent_ol_", a, ":", b)) + !!sym(paste0("percent_ol_", b, ":", a))), .by = type) %>%
        left_join(x, ., by = "type") -> x

    }
  }

  return(x)

}




secretome <- secretome %>%
  mutate(features = pmap(list(features, aa_scores), .f = add_hsr)) %>%
  mutate(features = map(features, .f = combine_phs)) %>%
  mutate(features = map2(features, topo, .f = \(x, y) limit_phs(x, y))) %>%
  mutate(features = map(features, .f = combine_phs_hsr))


secretome <- secretome %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "gpcrdb_gtp"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs", b = "gpcrdb_gtp"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "phs"))) %>%
  mutate(features = map(features, .f = ~combine_a_b(x = ., a = "phs_c-term", b = "phs_c-term-Cys", tag = "phs_C"))) %>%
  mutate(features = map(features, .f = ~combine_a_b(x = ., a = "phs_n-term", b = "phs_n-term-Cys", tag = "phs_N"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "phs_N"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "phs_C")))


secretome <- secretome %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs", b = "top200NC")))


secretome <- secretome %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs", b = "uni_pep")))


score_regions_h <- function(start, end, to_score, aa_scores) {

  entire <- start:end
  len <- end - start


  mid <- length(entire) %/% 2
  mid <- start + mid


  bind_cols(
    lapply(to_score, function(y) {
      to_average <- as.numeric(aa_scores[[y]])
      tmp <- tibble(!!paste0(y,"_entire") := mean(to_average[entire], na.rm = TRUE),
             !!paste0(y, "_nterm") := mean(to_average[start:mid], na.rm = TRUE),
             !!paste0(y, "_cterm") := mean(to_average[mid:end], na.rm = TRUE))
      if(do.call(c, tmp) %>% {all(is.na(.))}) {
        tmp[[paste0(y, "_max")]] <- NA
      } else {
        tmp <- tmp %>%
             mutate(!!paste0(y, "_max") := max(c_across(everything()), na.rm = TRUE))
      }
      return(tmp)
    }
    ))


}

score_regions <- function(features,
                          aa_scores,
                          to_score = c(basic, angles, energy, "DB_Dibasic",  "SV_sequence variant",
                                       paste0(c("chem", "pep"), "_", nn[["neural_net"]]),
                                       paste0(c("chem", "pep"), "_", sub("nn", "xgb", nn[["neural_net"]]))),
                          n_v_c = "pep_xgb4c") {

  roi_to_score <- c("gpcrdb_gtp", "sven", "phs", "top200NC", "peptide")

  if("peptide" %in% roi_to_score) {
    rts <- bind_rows(features %>%
                       filter(source == "uniprot" & type == "peptide")  %>%
                       drop_na(start,end),
                     features %>%
                       filter(source %in% roi_to_score) %>%
                       drop_na(start,end)
    )

    tmp2 <- features %>% filter(!(source == "uniprot" & type == "peptide")) %>% filter(!source %in% roi_to_score)
  } else {
    rts <- features %>%
      filter(source %in% roi_to_score) %>%
      drop_na(start,end)

    tmp2 <- features %>%
      filter(!source %in% roi_to_score)
  }


  if(nrow(rts) > 0) {
    features <- bind_rows(tmp2,
                          rts %>%
                            rowwise %>%
                            mutate(score_regions_h(start, end,
                                                   to_score = to_score,
                                                   aa_scores = aa_scores)) %>%
                            ungroup %>%
                            mutate(max_score = pmax(!!!syms(paste0(n_v_c, c("_entire", "_nterm", "_cterm"))))) %>%
                            mutate(hsterm = ifelse(.data[[paste0(n_v_c, "_nterm")]] > .data[[paste0(n_v_c, "_cterm")]], "N", "C")) %>%
                            mutate(dterm = abs(.data[[paste0(n_v_c, "_nterm")]] - .data[[paste0(n_v_c, "_cterm")]])))
  }

  return(features)

}

gc()
gc()
gc()
gc()
start <- Sys.time()

future::plan(strategy = future::multisession(workers = 8))

secretome[["features"]] <- furrr::future_map2(secretome[["features"]], secretome[["aa_scores"]], .f = score_regions)

end <- Sys.time()
end - start

secretome <- secretome %>%
              mutate(features = map(features, \(x) x %>%
                                                    mutate(roi_length = end - start + 1)))





features <- secretome %>% filter(accession == "P01275") %>% pull(features) %>% `[[`(1)

sequence_uni <- secretome %>% filter(accession == "P01275") %>% pull(sequence_uni)

gene <- secretome %>% filter(accession == "P01275") %>% pull(gene)

features <- secretome[["features"]][[2]]

sequence_uni <- secretome[["sequence_uni"]][[2]]


NTC = list(`P9` = c("P", "L", "S", "G", "E", "A"),
           `P8` = c("A", "L", "S", "E", "R"),
           `P7` = c("G", "E", "L", "S", "P", "K", "A", "R"),
           `P4` = c("P", "E", "L", "S", "R"),
           `P3` = c("L", "R", "A"),
           `P2` = c("G", "E", "S", "A", "Q", "P", "V", "L"),
           `P3'` = c("R", "S", "K", "L"),
           `P4'` = c("G", "R", "S", "P", "D", "E", "N"))

CTC = list(`P7` = c("G", "L", "E"),
           `P5` = c("G"),
           `P4` = c("L", "T", "P"),
           `P2` = c("G", "L", "S", "E", "Q", "A", "V", "T", "M"),
           `P3'` = c("A", "L", "G", "S", "R", "E"),
           `P4'` = c("G", "E", "P", "S", "R", "A", "T"),
           `P5'` = c("F", "G", "S", "P", "E", "D", "R", "K", "L"),
           `P6'` = c("L", "E", "S", "K", "V"),
           `P7'` = c("R", "L", "S", "A", "K", "E", "Q", "P"),
           `P8'` = c("R", "K", "L", "E", "A", "P"),
           `P9'` = c("R", "E", "S", "P", "A"),
           `P10'` = c("K", "G", "D", "P", "R", "L", "E", "Q"))

score_dbc <- function(seq_dbcenter, start_dbc, shift) {

  if(is.na(seq_dbcenter)) {
    return(tibble(index = NA))
  } else {



  tmp <- tibble(index = c(paste0("P", 10:1), paste0("P", 1:10, "'")),
                     AA = stringr::str_split(seq_dbcenter, "", simplify = TRUE) %>% c) %>%
                      mutate(uniprot_index = (start_dbc - shift - 1) + row_number())

  tmp <- tmp %>%
      rowwise() %>%
      mutate(NTC = case_when(!index %in% names(NTC) ~ NA,
                             AA %in% NTC[[index]] ~ 1,
                             !AA %in% NTC[[index]] ~ 0)) %>%
      mutate(CTC = case_when(!index %in% names(CTC) ~ NA,
                             AA %in% CTC[[index]] ~ 1,
                             !AA %in% CTC[[index]] ~ 0)) %>%
      ungroup() %>%
      pivot_longer(cols = c("NTC", "CTC")) %>%
      filter(!is.na(value))

  pep_dat <- tmp %>%
      select(-AA, -uniprot_index) %>%
      pivot_wider(names_from = c("name", "index"), values_from =  value) %>%
      mutate(NTC_tot = sum(c_across(starts_with("NTC_")))/8) %>%
      mutate(CTC_tot = sum(c_across(starts_with("CTC_")))/12) %>%
      mutate(dbc_term = if_else(NTC_tot >= CTC_tot, "N", "C"))

  pep_dat[["dbc_data"]] <- tmp %>%
                            filter(value == 1)  %>%
                            list()

  return(pep_dat)
  }
}


get_db_sites <- function(features, sequence_uni,
                         db_pep_source = c("gpcrdb_gtp", "phs", "sven", "peptide"), wN = 10, wC = 10) {


  shift <- max(c(wN, wC))

  n_trunc <- features %>%
    filter(type == "signal peptide") %>%
    pull(end)

  if(length(n_trunc) == 0) {
    n_trunc <- 1
  }
  if(is.na(n_trunc)) {
    n_trunc <- 1
  }
  if(length(n_trunc) > 1) {
    n_trunc <- max(n_trunc, na.rm = TRUE)
  }

  c_trunc <- nchar(sequence_uni)

  AA_seq <- stringr::str_sub(string = sequence_uni, start = n_trunc + 1, end = c_trunc) %>%
    stringr::str_c(stringr::str_c(rep("-", n_trunc + shift), collapse = ""), ., collapse = "") %>%
    stringr::str_c(., stringr::str_c(rep("-", shift), collapse = ""))


  if("peptide" %in% db_pep_source) {
    tmp <- bind_rows(features %>%
                      filter(source == "uniprot" & type == "peptide"),
                     features %>%
                       filter(source %in% db_pep_source)
    )

    tmp2 <- features %>% filter(!(source == "uniprot" & type == "peptide")) %>% filter(!source %in% db_pep_source)
  } else {
  tmp <- features %>%
    filter(source %in% db_pep_source)

  tmp2 <- features %>%
    filter(!source %in% db_pep_source)
  }

  if(nrow(tmp) > 0) {

    tmp <- tmp %>%
      pivot_longer(cols = c("start", "end"), values_to = "position", names_to = "terminus")


    tmp <- tmp %>%
      mutate(case_when(terminus == "start" ~ tibble(start_db = position - (wN - shift - 1),
                                                    end_db = position + (shift + wC)),
                       terminus == "end" ~ tibble(start_db = position - (wN - shift - 1),
                                                  end_db = position + (shift + wC)))) %>%
      mutate(seq_cleavage = map2_chr(start_db, end_db, \(x, y) {stringr::str_sub(AA_seq, start = x, end = y)})) %>%
      mutate(has_db = map_lgl(seq_cleavage, ~stringr::str_detect(., "KR|RK|RR|KK"))) %>%
      mutate(has_db3 = map_lgl(seq_cleavage, ~stringr::str_sub(., 7, -7) %>% stringr::str_detect(., "KR|RK|RR|KK"))) %>%
      mutate(db_loc = stringr::str_locate_all(seq_cleavage, "KR|RK|RR|KK"))


    center_db <- function(db_loc) {

      if(nrow(db_loc) == 0) {
        return(NA)
      } else {

        return(db_loc %>%
                 as_tibble %>%
                 mutate(offset = min(abs(((shift / 2) + 0.5) - start), abs(((shift / 2) + 0.5) - end))) %>%
                 filter(offset == min(offset)) %>%
                 mutate(nudge = start - shift) %>%
                 slice(1) %>%
                 pull(nudge)
        )

      }
    }

    tmp <- tmp %>%
      mutate(nudge = map_dbl(db_loc, center_db)) %>%
      mutate(case_when(terminus == "start" ~ tibble(start_dbc = position - (wN - shift - 1) + nudge,
                                                    end_dbc = position + (shift + wC) + nudge),
                       terminus == "end" ~ tibble(start_dbc = position - (wN - shift - 1) + nudge,
                                                  end_dbc = position + (shift + wC) + nudge))) %>%
      mutate(seq_dbcenter = map2_chr(start_dbc, end_dbc, \(x, y) {stringr::str_sub(AA_seq, start = x, end = y)}))

    tmp <- tmp %>%
      rowwise() %>%
      mutate(bind_rows(score_dbc(seq_dbcenter, start_dbc, shift))) %>%
      ungroup()

    tmp <- tmp %>%
      {cols_after_pos <<- names(.)[(match("position", names(.))):ncol(.)]; .} %>%
      pivot_wider(names_from = "terminus", values_from = all_of(cols_after_pos))

    tmp <- tmp %>%
      dplyr::rename(start = position_start,
                    end = position_end)

    new_feats <- tmp %>%
      filter(has_db_start | has_db_end)

    if(nrow(new_feats) > 0) {
    new_feats <- new_feats %>%
      mutate(dbc_data = map2(dbc_data_start, dbc_data_end, \(x, y) {

        if(!is.null(x)) x <- x %>% mutate(terminus = "start")
        if(!is.null(y)) y <- y %>% mutate(terminus = "end")

        bind_rows(x, y)

      })) %>%
      select(type, evidence, source, description, start, end, dbc_data) %>%
      dplyr::rename_with(.fn = ~paste0(., "_peptide"))

    new_feats <- new_feats %>%
      unnest(dbc_data_peptide) %>%
      distinct(uniprot_index, name, terminus) %>%
      dplyr::rename(start = uniprot_index) %>%
      mutate(type = paste(name, terminus, sep = "_"),
             end = NA,
             evidence = "dbc",
             source = "dbc")
    } else {
      new_feats <- NULL
    }

    return(bind_rows(tmp2, tmp, new_feats))


  } else {
    return(features)
  }


}


get_db_sites(features, sequence_uni)

gc()
gc()
start <- Sys.time()
future::plan(strategy = future::multisession(workers = 8))

secretome[["features"]] <- furrr::future_map2(.x = secretome[["features"]],
                                              .y = secretome[["sequence_uni"]], .f = \(x, y) get_db_sites(features = x, sequence_uni = y))

end <- Sys.time()
end - start



secretome <- secretome %>%
  mutate(bind_rows(map2(features, sequence_uni, \(x, y) {

    tibble(tot_db = x %>%
              filter(source == "sites" & type == "Dibasic") %>%
                nrow) %>%
      mutate(db_per_AA = tot_db/nchar(y))

  })))






















