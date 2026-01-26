

###add KK/KR and other features


secretome <- readRDS("~/peptide_alg/secretome_5.rds")

secretome_aa <- readRDS("~/peptide_alg/secretome_aa.rds")

to_join <- secretome_aa %>%
  group_by(accession) %>%
  nest(.key = "aa_scores")

secretome_aa2 <- readRDS("~/peptide_alg/og_secretome/secretome_aa.rds")

to_join2 <- secretome_aa2 %>%
  group_by(accession) %>%
  nest(.key = "aa_scores2")

to_join <- left_join(to_join, to_join2, by = "accession") %>%
  mutate(aa_scores = map2(aa_scores, aa_scores2, .f = \(x, y) {
    if(!is.null(y)) {
    bind_cols(x, y %>%
                select(starts_with("score_nn4")) %>%
                rename_with(.cols = starts_with("score_nn4"), .fn = ~sub("score_nn4", "score_nn4og", .)) %>%
                mutate(across(starts_with("score_nn"), .fns = ~scales::rescale(., to = c(0,1)), .unpack = TRUE)))
    } else {
      x
    }
  })) %>%
  select(-aa_scores2)



secretome <- left_join(secretome, to_join, by = "accession")



add_hsr <- function(features, aa_scores, score = roi_color, cutoff = 0.75) {

  x <- aa_scores[[score]]
  x[is.nan(x)] <- 0
  x <- x - min(x, na.rm = TRUE)
  x <- x/max(x, na.rm = TRUE)

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

limit_phs <- function(x,
                      offlimits = c("signal peptide", "signal peptide (manual)", "E", "Dibasic"),
                      offlimits_alt = c("signal peptide", "signal peptide (manual)", "strand", "Dibasic")) {
  if(!"alpha fold" %in% x[["source"]]) {
    offlimits <- offlimits_alt
  }

    comb_dat <- x %>%
    filter(type %in% offlimits | source == "phs")


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

    if("gtp" %in% x[["source"]] & "phs_hsr" %in% x[["source"]]) {

    ol_dat <- x %>%
    filter(source %in% c("gtp", "phs_hsr")) %>%
    distinct(start, end, .keep_all = TRUE) %>%
    arrange(start) %>%
    select(type, start, end, source) %>%
    drop_na() %>%
    split(., .[["source"]]) %>%
    {cross_join(.[["gtp"]], .[["phs_hsr"]], suffix = c("_gtp", "_phs_hsr"))} %>%
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
  mutate(features = map(features, .f = limit_phs)) %>%
  mutate(features = map(features, .f = combine_phs_hsr))


secretome <- secretome %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "gtp"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "phs"))) %>%
  mutate(features = map(features, .f = ~combine_a_b(x = ., a = "phs_c-term", b = "phs_c-term-Cys", tag = "phs_C"))) %>%
  mutate(features = map(features, .f = ~combine_a_b(x = ., a = "phs_n-term", b = "phs_n-term-Cys", tag = "phs_N"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "phs_N"))) %>%
  mutate(features = map(features, .f = ~define_ol(x = ., a = "phs_hsr", b = "phs_C")))






score_regions_h <- function(start, end, to_score, aa_scores) {

  entire <- start:end
  len <- end - start


  mid <- length(entire) %/% 2
  mid <- start + mid


  bind_cols(
    lapply(to_score, function(y) {
      to_average <- scales::rescale(aa_scores[[y]], to = c(0,1))
      tibble(!!paste0(y,"_entire") := mean(to_average[entire], na.rm = TRUE),
             !!paste0(y, "_nterm") := mean(to_average[start:mid], na.rm = TRUE),
             !!paste0(y, "_cterm") := mean(to_average[mid:end], na.rm = TRUE))
    }
    ))


}

score_regions <- function(features,
                          aa_scores,
                          to_score = c("conservation",
                                       "relASA",
                                       "pathogenicity",
                                       roi_color),
                          n_v_c = roi_color) {


  rts <- features %>%
    filter(source %in% c("gtp", "phs_hsr")) %>%
    drop_na(start,end)

  if(nrow(rts) > 0) {
    features <- rts %>%
    rowwise %>%
    mutate(score_regions_h(start, end,
                           to_score = to_score,
                           aa_scores = aa_scores)) %>%
    ungroup %>%
    mutate(max_score = pmax(!!!syms(paste0(n_v_c, c("_entire", "_nterm", "_cterm"))))) %>%
    mutate(hsterm = ifelse(.data[[paste0(n_v_c, "_nterm")]] > .data[[paste0(n_v_c, "_cterm")]], "N", "C")) %>%
    mutate(dterm = abs(.data[[paste0(n_v_c, "_nterm")]] - .data[[paste0(n_v_c, "_cterm")]])) %>%
    plyr::rbind.fill(features %>% filter(!source %in% c("gtp", "phs_hsr")), .)
  }

  return(features)

}

secretome <- secretome %>%
  mutate(features = pmap(list(features, aa_scores), .f = score_regions))


saveRDS(secretome, "~/peptide_alg/secretome_5.2.rds")




all_feat_cols <- unique(do.call(c, lapply(secretome$features, \(x) colnames(x))))



extract_roi <- \(x, y) {
  rois <- x %>%
    filter(source %in% c("gtp", "phs_hsr")) %>%
    dplyr::select(-start, -end, -evidence) %>%
    dplyr::rename("roi_name" = "type", "roi_type" = "source")

  dummy_table <- tibble(accession = y, roi_name = NA, roi_type = NA)
  if(nrow(rois) > 0) {
    return(bind_cols(accession = rep(y, nrow(rois)), rois))
  } else {
    return(dummy_table)
  }
}

left_side <- secretome %>%
  rowwise() %>%
  reframe(extract_roi(features, accession))

secretome_roi <- left_join(left_side, secretome, by = "accession")

secretome_roi <- secretome_roi %>%
  drop_na(roi_name) %>%
  mutate(roi_length = map2_int(features, roi_name, .f = \(features, roi) {
    to_ret <- features %>%
      filter(type == roi) %>%
      mutate(dif = end - start + 1) %>%
      pull(dif)
    return(to_ret[1])
  }))

secretome_roi <- secretome_roi %>%
  mutate(known_ligand = droplevels(factor(ifelse(roi_type == "gtp", "known", "putative"))))


















add_site_stats <- function(features, cons, dssp, sequence_uni) {

  regions <- grep("^dbr_", features[["type"]], value = TRUE)

  regions2 <- features %>%
    filter(source == "gtp") %>%
    pull(type)

  regions <- c(regions, regions2)

  if(length(regions) == 0) {
    return(features)
  } else {

    aa_sequence <- strsplit(sequence_uni, "")[[1]]

    if("frequency" %in% names(cons$ms)) {
      conservation <- cons$ms$frequency
    } else {
      conservation <- rep(NA, length(aa_sequence))
    }

    if("relASA" %in% names(dssp$ms)) {
      asa <- dssp$ms$relASA
    } else {
      asa <- rep(NA, length(aa_sequence))
    }

    fois <- features %>% filter(type == "W or Y")

    for(y in regions) {
      site <- unlist(features[features[["type"]] == y,c("start", "end")], use.names = FALSE)

      res <- fois %>%
        mutate(overlap = start <= site[2] & end >= site[1]) %>%
        filter(overlap) %>%
        mutate(WY_cons = conservation[start],
               WY_asa = asa[start],
               WY_total = rep(nrow(.), nrow(.)),
               WY_site = paste0(aa_sequence[start], start))
      if(nrow(res) > 0) {
        res <- res %>%
          filter(WY_cons == min(WY_cons)) %>%
          dplyr::slice(1) %>%
          dplyr::select(WY_cons, WY_asa, WY_total, WY_site)
      }

      if(nrow(res) == 1) {
        features[features[["type"]] == y, c("WY_cons", "WY_asa", "WY_total", "WY_site")] <- res
      }

      if(grepl("^dbr_", y)) {
        db_inds <- as.integer(strsplit(strsplit(y, "_")[[1]][2], "-")[[1]])

        features[features[["type"]] == y, "Dibasic_cons"] <- mean(conservation[db_inds], na.rm = TRUE)
        features[features[["type"]] == y, "Dibasic_asa"] <-  mean(asa[db_inds], na.rm = TRUE)
      }
    }

  }

  return(features)

}

uniprot_t <- uniprot_t %>%
  mutate(features = pmap(list(features, cons, dssp, sequence_uni), .f = add_site_stats))



