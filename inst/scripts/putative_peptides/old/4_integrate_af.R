
###add alphafold structural data ~20min
##https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
##use dssp.py to extract ss data

library(Biostrings)

secretome <- readRDS("~/peptide_alg/secretome.rds")

af_dir <- "~/peptide_alg/alphafold_dssp"

af_dat <- tibble(af_files = list.files(af_dir),
                 accession = str_extract(af_files, "(?<=-)[^-]+(?=-)")) %>%
  mutate(af_dat = map(af_files, \(x) {
    tmp <- data.table::fread(paste0(af_dir, "/", x))
    tmp <- as_tibble(tmp)
    colnames(tmp) <- c("index", "AA", "SS", "relASA", "Phi", "Psi", "NH->O_1_relidx", "NH->O_1_energy", "O->NH_1_relidx", "O->NH_1_energy", "NH->O_2_relidx", "NH->O_2_energy", "O->NH_2_relidx", "O->NH_2_energy")
    return(tmp)
  })) %>%
  distinct(accession, .keep_all = TRUE)


secretome <- left_join(secretome, af_dat, by = "accession")

secretome <- secretome %>%
  mutate(sequence_af = map_chr(af_dat, ~paste(.[["AA"]], collapse = ""))) %>%
  mutate(af_dat = map(.x = af_dat, .f = \(x) {
    if(!is.null(x)) {
      return(
        x %>%
          mutate(relASA_s = slider::slide_dbl(relASA, mean, .before = 5, .after = 5)) %>%
          mutate(relASA_ss = slider::slide_dbl(relASA_s, mean, .before = 3, .after = 3)) %>%
          mutate(relASA_sss = slider::slide_dbl(relASA_ss, mean, .before = 3, .after = 3))
      )
    } else {
      return(x)
    }
  })) %>%
  mutate(uni_af_exact = sequence_af == sequence_uni)


# seq1 <- pull(secretome[198,"sequence_af"])
# seq2 <- pull(secretome[198, "sequence_uni"])
# to_map <- secretome[198,"af_dat"][[1]][[1]]


start <- Sys.time()

secretome <- secretome %>%
  mutate(af_mapped = pmap(.l = list(seq1 = sequence_af, seq2 = sequence_uni, to_map = af_dat),
                          .f = map_table))

end <- Sys.time()
end - start


# features <- secretome[3, "features"][[1]][[1]]
# af_mapped <- secretome[3, "af_mapped"][[1]][[1]]

add_af_features <- function(features, af_mapped) {
  
  tf <- af_mapped[["ms"]]
  if(!"SS" %in% names(tf)) {
    return(features %>%
             mutate(source = "uniprot") %>%
             mutate(start = as.integer(start),
                    end = as.integer(end)))
  } else {
    tf <- tf[["SS"]]
    uni_tf <- unique(tf)
    af_feats <- bind_rows(
      lapply(uni_tf[!uni_tf == "-"], \(x) {
        vec <- which(tf == x)
        
        breaks <- c(0, which(diff(vec) != 1), length(vec))
        
        sequences <- lapply(seq_along(breaks[-1]), \(i) vec[(breaks[i] + 1):breaks[i + 1]])
        
        bind_rows(lapply(sequences, \(y) tibble(type = x,
                                                evidence = "af",
                                                start = min(y),
                                                end = max(y),
                                                source = "alpha fold")
        ))
        
      })
    )
    
    features <- features %>%
      mutate(source = "uniprot") %>%
      mutate(start = as.integer(start),
             end = as.integer(end))
    
    return(bind_rows(features, af_feats))
  }
}

secretome <- secretome %>%
  mutate(features = pmap(list(features, af_mapped), add_af_features))

saveRDS(secretome, "~/peptide_alg/secretome_4.rds")


