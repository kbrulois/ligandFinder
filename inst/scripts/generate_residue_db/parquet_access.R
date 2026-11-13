




pq_path <- "~/ligandFinder_data/residue_db"

res_db <- arrow::open_dataset(source = pq_path)

proteins <- c("CXL14", "TMM70", "RM27")

residue_data <- res_db %>%
  filter(uni_gene %in% proteins) %>%
  collect()

residue_data$features[[2]] %>%
  filter(source == "uniprot") %>%
  View(.)

residue_data$sequence_uni %>% nchar


test <- rep(0, 644)

test[1:18] <- 1

