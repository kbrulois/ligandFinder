###combine uniprot and aminode ~20min

non_secreted_genes <- adhesionGPCRs <- list(`Group I` = paste0("ADGRG", 1:6),
                                            `Group II` = paste0("ADGRL", 1:3),
                                            `Group III` = paste0("ADGRB", 1:3),
                                            `Group IV` = paste0("ADGRE", 1:5),
                                            `Group V` = paste0("ADGRF", 1:5),
                                            `Group VI` = paste0("ADGRA", 1:3),
                                            `Group VII` = paste0("ADGRV1"),
                                            PARS = paste0("PAR", 1:3))

non_secreted_genes <- unname(do.call(c, non_secreted_genes))



uniprot_t <- readRDS("~/peptide_alg/uniprot.rds")

secretome_og <- data.table::fread("~/peptide_alg/sa_location_Secreted HPA 2793 genes-1.tsv")

secretome <- uniprot_t %>%
  mutate(secreted_any = map_lgl(.x = annotations, .f = \(x) {
    x %>%
      filter(annotation_name == "comment" &
             annotation_type == "subcellular location" &
             name_2 == "location") %>%
      {any(.[["annotation"]] == "Secreted")}
  })) %>%
  mutate(secreted_all = map_lgl(.x = annotations, .f = \(x) {
    x %>%
      filter(annotation_name == "comment" &
               annotation_type == "subcellular location" &
               name_2 == "location") %>%
      {all(.[["annotation"]] == "Secreted")}
  })) %>%
  mutate(secreted_HPA = gene %in% secretome_og[["Gene"]]) %>%
  mutate(non_secreted_goi = gene %in% non_secreted_genes) %>%
  mutate(secreted_final = secreted_any | secreted_HPA | non_secreted_goi) %>%
  filter(secreted_final)

rm(uniprot_t, secretome_og)


aminode <- readRDS("~/peptide_alg/aminode.rds")

aminode[aminode$gene == "C10ORF99", "gene"] <- "GPR15LG"

aminode <- aminode %>%
  mutate(cons = map(.x = cons, .f = \(x) {
    x %>%
      filter(AA != "N/A" | index != "N/A") %>%
      mutate(index = as.integer(index),
             cons = as.numeric(cons)) %>%
      mutate(frequency = cons/mean(cons)) %>%
      mutate(window = slider::slide_dbl(frequency, mean, .before = 5, .after = 5)) %>%
      mutate(smooth = slider::slide_dbl(window, mean, .before = 3, .after = 3)) %>%
      mutate(doubleSmooth = slider::slide_dbl(smooth, mean, .before = 3, .after = 3))
    
  })) %>%
  mutate(sequence = map_chr(.x = cons, .f = ~ paste(.[["AA"]], collapse = "")))

secretome <- left_join(secretome, aminode, by = "gene", suffix = c("_uni", "_ami"))


library(Biostrings)

map_table <- function(seq1, seq2, to_map) {
  
  if(is.null(to_map) | seq1 == "") {
    return(list(ms = rep(NA, nchar(seq2)),
                score = NA))
  } else {
    dd <- rep(NA, nchar(seq2))
    catch_mapping <- to_map %>% reframe(across(everything(), \(x) dd))
    
    seq1 <- Biostrings::AAString(seq1)
    seq2 <- Biostrings::AAString(seq2)
    
    alignment <- Biostrings::pairwiseAlignment(seq1, seq2, type = "global")
    
    #print(alignment)
    
    norm_score <- Biostrings::score(alignment) / width(Biostrings::alignedPattern(alignment))
    
    aligned_seq1 <- strsplit(as.character(pattern(alignment)), "")[[1]]
    aligned_seq2 <- strsplit(as.character(subject(alignment)), "")[[1]]
    
    dd <- rep(NA, length(aligned_seq2))
    mapped_data <- to_map %>% reframe(across(everything(), \(x) dd))
    
    data_index <- start(pattern(alignment))
    
    for (i in seq_along(aligned_seq1)) {
      nonDash <- aligned_seq1[i] != "-" && aligned_seq2[i] != "-"
      seq2Dash <- aligned_seq1[i] != "-" && aligned_seq2[i] == "-"
      seq1Dash <- aligned_seq1[i] == "-" && aligned_seq2[i] != "-"
      if (nonDash) {
        mapped_data[i,] <- to_map[data_index,]
        data_index <- data_index + 1
      } else if(seq2Dash) {
        data_index <- data_index + 1
      } 
    }
    
    mapped_data <- mapped_data[aligned_seq2 != "-", ]
    
    alignment2 <- Biostrings::pairwiseAlignment(paste(mapped_data[["AA"]][!is.na(mapped_data[["AA"]])], collapse = ""), seq2, type = "global")
    
    
    begin <- start(subject(alignment2))
    end <- end(subject(alignment2))
    
    catch_mapping[begin:end, ] <- mapped_data
    
    return(list(ms = catch_mapping,
                score = norm_score))
  }
  
}

future::plan(strategy = future::sequential())

secretome <- secretome %>%
  mutate(cons_mapped = furrr::future_pmap(.l = list(seq1 = sequence_ami, seq2 = sequence_uni, to_map = cons),
                            .f = map_table))

gc()

secretome <- secretome %>% 
  mutate(cons_map_score = map_dbl(cons_mapped, .f = \(x) {
    if(is.null(x)) {
      return(NA)
      } else {
        return(x[["score"]])
    }})) 

saveRDS(secretome, "~/peptide_alg/secretome.rds")

