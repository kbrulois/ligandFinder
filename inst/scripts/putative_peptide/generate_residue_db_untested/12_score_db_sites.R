




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












secretome <- secretome %>%
  mutate(features = pmap(list(features, cons, dssp, sequence_uni), .f = add_site_stats))


