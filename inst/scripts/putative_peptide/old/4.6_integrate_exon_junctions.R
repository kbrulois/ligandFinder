


secretome <- readRDS("~/peptide_alg/secretome_4.5.rds")

# split_ens <- function(x) {
#   tmp <- strsplit(x, "\\|")[[1]]
#   tibble(ens_gene = tmp[1],
#          ens_transcript = tmp[2],
#          ens_protein = tmp[3])
# }
# 
# secretome <- secretome %>%
#   rowwise() %>%
#   mutate(split_ens(ensemble))


# library(biomaRt)
# 
# ensembl <- useMart("ensembl")
# human_ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)
# 
# cdsAnnot <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","ensembl_peptide_id","strand","gene_biotype","cds_start","cds_end","ensembl_exon_id"),
#                   filters = 'uniprot_gn_id',
#                   values = secretome$accession,
#                   mart = human_ensembl)
# 
# key <- getBM(attributes = c("ensembl_gene_id","uniprot_gn_id"), 
#              filters = "uniprot_gn_id", 
#              values = secretome$accession, 
#              mart = human_ensembl)
# 
# saveRDS(cdsAnnot, "~/peptide_alg/exon_juncs.rds")
# saveRDS(key, "~/peptide_alg/exon_junc_key.rds")

cdsAnnot <- readRDS("~/peptide_alg/exon_juncs.rds")

cdsAnnot <- cdsAnnot %>%
  as_tibble %>%
  mutate(ensemble = paste(ensembl_gene_id, ensembl_transcript_id, ensembl_peptide_id, sep = "|")) %>%
  filter(ensemble %in% secretome[["ensemble"]]) %>%
  group_by(ensemble) %>%
  summarize(ensemble = unique(ensemble), exon_coords = list(tibble(ensembl_exon_id = ensembl_exon_id,
                                                      cds_start = cds_start,
                                                      cds_end = cds_end))) %>%
  mutate(exon_coords = map(exon_coords, \(x) {
    x %>%
      distinct() %>%
      mutate(cds_end_p = round(cds_end/3, 2))
  }))
  


secretome <- left_join(secretome, cdsAnnot, by = "ensemble")
  

fetch_exon_length <- function(x) {
  if(!is.null(x)) {
    to_return <- max(x[["cds_end_p"]], na.rm = TRUE)
  } else {
    to_return <- NA
  }
  return(to_return)
}

add_exon_data <- function(ens_uni_ok, features, ens_coords) {
  logi <- ens_uni_ok
  logi[is.na(logi)] <- FALSE
  if(!logi) {
    return(features)
  } else {
    return(bind_rows(features,
    tibble(type = ens_coords[["ensembl_exon_id"]],
           evidence = "ensemble",
           start = ens_coords[["cds_end_p"]],
           end = ens_coords[["cds_end_p"]],
           source = "ensemble")
    ))
  }
}


secretome <- secretome %>%
  ungroup() %>%
  mutate(exon_length = map_dbl(.x = .[["exon_coords"]], .f = fetch_exon_length)) %>%
  mutate(ens_uni_diff = abs((exon_length - 1) - nchar(sequence_uni))) %>%
  mutate(ens_uni_ok = ens_uni_diff < 3) %>%
  mutate(features = pmap(.l = list(ens_uni_ok, features, exon_coords), .f = add_exon_data))
  
  
  
saveRDS(secretome, "~/peptide_alg/secretome_4.6.rds")

  




current_gtp <- bind_rows(lapply(seq_along(secretome[["features"]]), \(x) {
  secretome[["features"]][[x]] %>%
    filter(source == "gtp") %>%
    distinct(start, end, .keep_all = TRUE) %>%
    mutate(gene = secretome[["gene"]][[x]]) %>%
    drop_na(start, end)
}))

og_gtp <- current_gtp


