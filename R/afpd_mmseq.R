


proteins <- c("CXCL14", "AGTR1", "AGTR2", "BDKRB1", "BDKRB2")

pq_path <- paste0(ligandFinder:::get_db_path(), "/residue_db")

res_db <- arrow::open_dataset(source = pq_path)


  residue_anno <- res_db %>%
    #filter(gene %in% proteins) %>%
    group_by(gene_grp) %>%
    collect() %>%
    ungroup()

  residue_anno <- residue_anno %>%
    mutate(entry_name = sub("_HUMAN$", "", entry_name))

fasta <- do.call(c, lapply(1:nrow(residue_anno), \(x) c(paste0(">", residue_anno[["entry_name"]][x]),
                                    residue_anno[["sequence_uni"]][x])))

write_lines(fasta, file = "~/Desktop/human_proteome.fa")
