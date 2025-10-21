


map_table <- function(seq1 = paste(to_map[["AA"]], collapse = ""),
                      seq2,
                      to_map) {

  if(is.null(to_map) | seq1 == "") {
    return(list(ms = rep(NA, nchar(seq2)),
                score = NA))
  } else if(identical(seq1, seq2)) {
    return(list(ms = to_map,
                score = 5))
  } else {

    dd <- rep(NA, nchar(seq2))
    catch_mapping <- to_map %>% reframe(across(everything(), \(x) dd))

    seq1 <- Biostrings::AAString(seq1)
    seq2 <- Biostrings::AAString(seq2)

    alignment <- Biostrings::pairwiseAlignment(seq1, seq2, type = "global-local", gapOpening = -10, gapExtension = -1)

    print(alignment)

    norm_score <- Biostrings::score(alignment) / Biostrings::width(Biostrings::alignedPattern(alignment))

    aligned_seqs <- tibble(aligned_seq1 = strsplit(as.character(Biostrings::pattern(alignment)), "")[[1]],
                           aligned_seq2 = strsplit(as.character(Biostrings::subject(alignment)), "")[[1]]) %>%
      filter(aligned_seq2 != "-")

    alignment2 <- Biostrings::pairwiseAlignment(Biostrings::AAString(paste0(aligned_seqs$aligned_seq1[aligned_seqs$aligned_seq1 != "-"], collapse = "")), seq1, type = "global", gapOpening = -10, gapExtension = -1)
    aligned_seqs2 <- tibble(aligned_seq1 = strsplit(as.character(Biostrings::pattern(alignment2)), "")[[1]],
                            aligned_seq2 = strsplit(as.character(Biostrings::subject(alignment2)), "")[[1]]) %>%
      filter(aligned_seq2 != "-")

    begin1 <- Biostrings::start(Biostrings::subject(alignment2))
    final1 <- Biostrings::end(Biostrings::subject(alignment2))

    begin2 <- Biostrings::start(Biostrings::subject(alignment))
    final2 <- Biostrings::end(Biostrings::subject(alignment))

    catch_mapping[begin2:final2,][aligned_seqs$aligned_seq1 != "-", ] <- to_map[begin1:final1, ][aligned_seqs2$aligned_seq1 != "-", ]

    return(list(ms = catch_mapping,
                score = norm_score))
  }

}


map_AA_sequence_to_uniprot <- function(x,
                                       species = "human",
                                       jar_path = "~/programs/PeptideMatchCMD_1.1.jar",
                                       sprot_path = "~/peptide_alg/sprot_index",
                                       out_path = "~/peptide_alg/peptide_match_out.txt") {

  system2(command = "java",
          args = c(paste("-jar", jar_path),
                   "-a query",
                   paste("-i", sprot_path),
                   paste("-q", paste(x, collapse = ",")),
                   paste("-o", out_path))
  )

  res <- data.table::fread(out_path,
                    sep = "\t",
                    skip = 1,
                    fill = TRUE) %>%
    as_tibble %>%
    rename(sequence = `##Query`) %>%
    mutate(species = stringr::str_extract(`Subject`, "[^_]*$")) %>%
    filter(species == toupper(!!species))

  if(nrow(res) < 1) {
    res[1, ] <- NA
  }
  return(res)

}
