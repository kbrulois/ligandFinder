

pq_path <- paste0(get_db_path(), "/residue_db")


secretome <- arrow::open_dataset(source = pq_path)

combine_scores <- function(scoreA, scoreB) {
  penalty <- ifelse(scoreB < 0.75, (1.6*(0.75 - scoreB))^2,
                    ifelse(scoreB > 0.9, (1.6*(scoreB - 0.9))^2, 0))

  combined_score <- scoreA * (1 - penalty)
  return(combined_score)
}


ss_features_coded <- 1:length(ss_features)
names(ss_features_coded) <- names(ss_features)
ss_features_coded <- c(setNames(0, "-"), ss_features_coded)


offlimits_features <- c("signal peptide", "E")

black_balled <- c("thrombin light chain", "epiregulin", "INSL5 (A chain)")

to_expand <- c("gtp", "signal peptide", "phs_dbr", "sites")

interval_to_factor <- function(type, start, end, dat_size) {
  to_return <- rep(NA, dat_size)
  for(i in seq_along(type)) {
    start_na <- is.na(start[i]) | start[i] > dat_size
    end_na <- is.na(end[i]) | end[i] > dat_size
    if(!start_na & end_na) {
      to_return[start[i]] <- type[i]
    } else if(!start_na & !end_na) {
      to_return[start[i]:end[i]] <- type[i]
    }
  }
  return(to_return)
}


secretome %>%
  mutate(factorized_features = map2(.x = features, .y = sequence_uni, .f = \(.x, .y) {

    dat <- tibble(.rows = nchar(.y))

    for(feat in to_expand) {
      dat[[feat]] <- .x %>%
        filter(source == feat) %>%
        {factorize_intervals_c(.[["type"]], .[["start"]], .[["end"]], nrow(dat))}
    }

    dat

  })) %>%
  select(29) -> test

secretome %>%
  select(-af_xyz) %>%
  unnest(cols = any_of(lists_to_unpack), names_sep = "_") %>%
  select(c(21,25,30,32, 40,49,50,60, 61, 62))


expand_by_AA <- \(sequence_uni, accession, gene, af_mapped, features, cons, af_missense, af_xyz) {

  aa_sequence <- strsplit(sequence_uni, "")[[1]]

  dat_size <- length(aa_sequence)

  if("frequency" %in% names(cons$ms)) {
    conservation <- cons$ms
  } else {
    dd <- rep(NA, nchar(seq2))
    catch_mapping <- to_map %>% reframe(across(everything(), \(x) dd))
    conservation <- rep(NA, dat_size)
  }

  if("relASA" %in% names(af_mapped$ms)) {
    asa <- af_mapped$ms %>%
      select(-AA, -index, -relASA_s, -relASA_ss, -relASA_sss) %>%
      dplyr::rename_with(\(x) {
        ifelse(!x %in% c("SS", "relASA"), paste0("AF_", x), x)
      })
  } else {
    asa <- Map(\(x) {rep(NA, dat_size)},
               c("SS", "relASA", "AF_Phi", "AF_Psi", "AF_NH->O_1_relidx", "AF_NH->O_1_energy",
                 "AF_O->NH_1_relidx", "AF_O->NH_1_energy", "AF_NH->O_2_relidx",
                 "AF_NH->O_2_energy", "AF_O->NH_2_relidx", "AF_O->NH_2_energy")) %>%
      as_tibble
  }

  if("mean_af_missense" %in% names(af_missense[["ms"]])) {
    afm <- af_missense$ms$mean_af_missense
  } else {
    afm <- rep(NA, dat_size)
  }

  if("AA" %in% names(af_xyz[["ms"]])) {
    afxyz <- af_xyz[["ms"]] %>% select(-AA)
  } else {
    afxyz <- rep(NA, dat_size)
  }

  dat <- tibble(AA = aa_sequence,
                conservation = as.numeric(conservation),
                pathogenicity = as.numeric(afm),
                accession = accession,
                gene = gene,
                known_peptide = 0,
                known_peptide_n = 0,
                known_peptide_c = 0,
                signal_peptide_or_Strand = 0,
                dibasic_Cys_W_T = as.character(NA))

  dat <- bind_cols(dat, asa, afxyz)

  gtp <- features %>%
    filter(source == "gtps")

  if(nrow(gtp) > 0) {
    gtp_pep <- gtp %>%
      rowwise() %>%
      mutate(locations = map2(start, end, .f = \(x, y) x:y)) %>%
      pull(locations) %>%
      do.call(c, .) %>%
      unique(.)

    dat[["known_peptide"]][gtp_pep] <- 1
    dat[["known_peptide_n"]][gtp[["start"]]] <- 1
    dat[["known_peptide_c"]][gtp[["end"]]] <- 1
  }

  # offlimits <- features %>%
  #   filter(type == "signal peptide" | type == "E") %>%
  #   rowwise() %>%
  #   mutate(locations = map2(start, end, .f = \(x, y) {if(!is.na(x) & !is.na(y)) x:y})) %>%
  #   pull(locations) %>%
  #   do.call(c, .) %>%
  #   unique(.)
  #
  # offlimits <- offlimits[offlimits %in% 1:nrow(dat)]
  #
  # dat[["signal_peptide_or_Strand"]][offlimits] <- 1

  for(feat in to_expand) {
    dat[[feat]] <- features %>%
      filter(source == feat) %>%
      {factorize_intervals(.[["type"]], .[["start"]], .[["end"]], dat_size)}
  }

  dat

}



start <- Sys.time()

secretome_aa <- secretome %>%
  rowwise %>%
  reframe(expand_by_AA(sequence_uni, accession, gene, af_mapped, features, cons, af_missense, af_xyz))

end <- Sys.time()
end - start



secretome <- secretome %>%
  select(-af_dat) %>%
  rename(dssp = af_mapped)




curl::curl_download(
  "https://drive.usercontent.google.com/download?id=1jPkj0fK0QIQM5JKcosQ-mZ_VZjKm57_w&authuser=0&confirm=xxx",
  destfile = "~/peptide_alg/residue_data.parquet",
  quiet = FALSE)


h <- curl::new_handle()
curl::handle_setopt(
  handle = h,
  httpauth = 1,
  userpwd = "kbrulois@globusid.org:globus_123"
)

curl::curl_download(
"https://g-c32f18.06752.75bc.data.globus.org/sec_parquet.tar.gz",
destfile = "~/peptide_alg/sec_parquet1.tar.gz", quiet = F, handle = h)

test2 <- arrow::open_dataset(source = pq_path)


chain_names <- c("BDKRB2", "CXCL14") #get from file name once formalized

residue_anno <- test2 %>%
  filter(gene %in% chain_names) %>%
  group_by(gene_grp) %>%
  collect()

residue_anno %>%
  select(-af_xyz) %>%
  unnest(cols = any_of(lists_to_unpack), names_sep = "_") %>%
  select(c(21,25,30,32, 40,49,50,60, 61, 62))




