
###add alphafold structural data ~20min
##https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
##run dssp.py to extract ss data ~14h

library(Biostrings)

secretome <- readRDS(paste0(s_localDir, "/processed/secretome_1.rds"))

dssp_dir <- paste0(s_localDir, "/processed/alphafold_dssp")

dssp <- tibble(files = list.files(dssp_dir),
                 accession = str_extract(files, "(?<=-)[^-]+(?=-)")) %>%
  mutate(dssp = map(files, \(x) {
    tmp <- data.table::fread(paste0(dssp_dir, "/", x))
    tmp <- as_tibble(tmp)
    colnames(tmp) <- c("index", "AA", "SS", "relASA", "Phi", "Psi", "NH->O_1_relidx", "NH->O_1_energy", "O->NH_1_relidx", "O->NH_1_energy", "NH->O_2_relidx", "NH->O_2_energy", "O->NH_2_relidx", "O->NH_2_energy")
    return(tmp)
  })) %>%
  distinct(accession, .keep_all = TRUE)


secretome <- left_join(secretome, dssp, by = "accession")

secretome <- secretome %>%
  mutate(sequence_dssp = map_chr(dssp, ~paste(.[["AA"]], collapse = ""))) %>%
  mutate(dssp = map(.x = dssp, .f = \(x) {
    if(!is.null(x)) {
      return(
        x %>%
          mutate(relASA_s = slider::slide_dbl(relASA, mean, .before = 3, .after = 3)) %>%
          mutate(relASA_ss = slider::slide_dbl(relASA_s, mean, .before = 3, .after = 3)) %>%
          mutate(relASA_sss = slider::slide_dbl(relASA_ss, mean, .before = 3, .after = 3))
      )
    } else {
      return(x)
    }
  })) %>%
  mutate(uni_af_exact = sequence_af == sequence_uni)


start <- Sys.time()

secretome <- secretome %>%
  mutate(dssp = pmap(.l = list(seq1 = sequence_af, seq2 = sequence_uni, to_map = dssp),
                          .f = map_table))

end <- Sys.time()
end - start


add_dssp_features <- function(features, dssp) {
  
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
  mutate(features = pmap(list(features, dssp), add_dssp_features))

saveRDS(secretome, paste0(s_localDir, "/processed/secretome_2.rds"))


