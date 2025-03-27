

library(httr)
library(jsonlite)

get_uniref <- function(uniref_id) {
  url <- paste0("https://rest.uniprot.org/uniref/", uniref_id, ".json")
  response <- GET(url)
  
  if (status_code(response) == 200) {
    data <- fromJSON(content(response, as = "text"))
    data[["representativeMember"]] <- data[["representativeMember"]][!names(data[["representativeMember"]]) == "sequence"]
    data <- data[!names(data) == "members"]
    data <- purrr::flatten(data)
    return(bind_cols(tibble(query = uniref_id), as_tibble(data, .name_repair = "minimal")))
    
  } else {
    return(structure(list(query = uniref_id, id = NA, name = NA, memberCount = NA, updated = NA, 
                          entryType = NA, scientificName = NA, taxonId = NA, seedId = NA, 
                          memberIdType = NA, memberId = NA, organismName = NA, organismTaxId = NA, 
                          sequenceLength = NA, proteinName = NA, accessions = NA, uniref50Id = NA, 
                          uniref90Id = NA, uniparcId = NA, seed = NA), class = c("tbl_df", 
                                                                                 "tbl", "data.frame"), row.names = c(NA, -1L)))
  }
}

get_uniparc <- function(uniref_id) {
  url <- paste0("https://rest.uniprot.org/uniparc/", sub("UniRef100_", "", uniref_id), ".json")
  message(url)
  response <- GET(url)
  
  if (status_code(response) == 200) {
    data <- fromJSON(content(response, as = "text"))
    data <- data[!names(data) %in% c("sequenceFeatures", "sequence", "locations")]
    data2 <- data[["uniParcCrossReferences"]]
    data2_org <- data2[["organism"]]
    data2 <- bind_cols(data2[!names(data2) == "organism"], data2_org)
    data <- purrr::flatten(data[!names(data) == "uniParcCrossReferences"])
    data <- data[!names(data) %in% c("value")]
    data <- c(data, data2)
    return(bind_cols(tibble(query = uniref_id), as_tibble(data, .name_repair = "minimal")))
    
  } else {
    message("didn't work")
    return(structure(list(query = uniref_id, uniParcId = NA, 
                          database = NA, id = NA, versionI = NA, 
                          version = NA, active = NA, created = NA, lastUpdated = NA, 
                          geneName = NA, proteinName = NA, 
                          length = NA, molWeight = NA, crc64 = NA, 
                          md5 = NA, oldestCrossRefCreated = NA, 
                          mostRecentCrossRefUpdated = NA, scientificName = NA, 
                          commonName = NA, taxonId = NA), class = c("tbl_df", 
                                                                    "tbl", "data.frame"), row.names = c(NA, -1L)))
  }
}

get_uniprotkb <- function(uniref_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", sub("UniRef100_", "", uniref_id), ".json")
  message(url)
  response <- GET(url)

  
  if (status_code(response) == 200) {
    data <- fromJSON(content(response, as = "text"))
    uniref_id2 <- data[["extraAttributes"]][["uniParcId"]]
    to_return <- get_uniparc(uniref_id = uniref_id2)
    to_return[["query"]] <- uniref_id
    return(to_return)
    
  } else {
    message("didn't work")
    return(structure(list(query = uniref_id, uniParcId = NA, 
                          database = NA, id = NA, versionI = NA, 
                          version = NA, active = NA, created = NA, lastUpdated = NA, 
                          geneName = NA, proteinName = NA, 
                          length = NA, molWeight = NA, crc64 = NA, 
                          md5 = NA, oldestCrossRefCreated = NA, 
                          mostRecentCrossRefUpdated = NA, scientificName = NA, 
                          commonName = NA, taxonId = NA), class = c("tbl_df", 
                                                                    "tbl", "data.frame"), row.names = c(NA, -1L)))
  }
}

map_all_ids <- function(uniref_id) {
  
    to_return <- get_uniref(uniref_id)
    
    if(!is.na(to_return[[2]][1])) {
    to_return[["status"]] <- "uniref"
    } else {
      
      message("uniref100 mapping failed. Trying uniparc")

      to_return <- get_uniparc(uniref_id)
      
      if(!is.na(to_return[[2]][1])) {
      to_return[["status"]] <- "uniparc"
      } else {
        message("uniparc mapping failed. trying uniprotkb")
        to_return <- get_uniprotkb(uniref_id)
        to_return[["status"]] <- "uniparc2"
      }
    }
      return(to_return)
}

map_tax_id <- function(taxid) {
  url <- paste0("https://rest.uniprot.org/taxonomy/", taxid, ".json")
  response <- GET(url)
  if (status_code(response) == 200) {
    data <- fromJSON(content(response, as = "text"))
    data2 <- data[["lineage"]]
    bind_rows(data2, tibble(scientificName = data[["scientificName"]],
                            taxonId = taxid,
                            rank = "species",
                            hidden = FALSE,
                            commonName = data[["commonName"]]))
  } else {
    message("tax not found")
    tibble(scientificName = NA,
           taxonId = taxid,
           rank = c("no rank", "superkingdom", "clade", "kingdom", "clade", "clade", 
                    "clade", "phylum", "subphylum", "clade", "clade", "clade", "clade", 
                    "superclass", "class", "subclass", "infraclass", "clade", "clade", 
                    "clade", "order", "family", "genus"),
           hidden = NA,
           commonName = NA)
  }
}


do_species_mapping <- function(ids) {
  
taxmy <- bind_rows(lapply(ids, map_all_ids))

taxmy <- taxmy %>% 
          group_by(query) %>%
          nest %>%
          ungroup

taxmy <- taxmy %>%
          mutate(map_df(data, \(x) {
            uni_taxids <- unique(x[["taxonId"]])
            tibble(taxonId = uni_taxids[1],
                   uni_taxids = length(uni_taxids))}))



taxmy <- taxmy %>%
          mutate(tax = map(taxonId, map_tax_id))

taxmy <- taxmy %>%
  mutate(tax = map(tax, \(x) {
    x %>%
      mutate(scientificName = if_else(rank == "species" & str_detect(scientificName, " "), sub("^([A-Za-z])\\w*", "\\1.", scientificName), scientificName)) %>%
      mutate(commonName = if_else(is.na(commonName), scientificName, commonName))
  }))

taxmy %>%
        mutate(map_df(tax, \(x) {
      
          x %>%
            filter(rank %in% c("phylum", "class", "species")) %>%
            select(rank, commonName) %>%
            group_by(rank) %>%
            summarise(across(everything(), first)) %>%
            pivot_wider(names_from = rank, values_from = commonName)
          
          }))

}

get_qc_data <- function(tax_dat) {
  
  tax_dat %>%
    mutate(map2_df(taxonId, data, \(taxonId, data) {
      data %>%
        select(any_of(c("proteinName", "status", "database", "geneName"))) %>%
        filter(taxonId == taxonId) %>%
        summarise(across(everything(), first))
    })) %>%
    mutate(qc = if_else(grepl("LOW QUALITY", proteinName), "low_quality", 
                        if_else(grepl("Fragment", proteinName), "fragment", NA)))
}



map_gene_to_uniprot <- function(gene_symbols, organism = "9606") {
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  
  # Convert gene symbols into a query string
  gene_query <- paste(paste0("gene:", gene_symbols), collapse = " OR ")
  
  # Construct full query with organism filter (default: human, 9606)
  query <- paste0("(", gene_query, ") AND organism_id:", organism)
  
  # API request
  response <- GET(url = base_url, query = list(query = query, format = "json", size = length(gene_symbols)))
  
  if (status_code(response) != 200) {
    stop("Error fetching data from UniProt")
  }
  
  # Parse response JSON
  data <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  # Extract relevant information
  results <- data$results
  if (is.null(results)) return(data.frame(Gene = gene_symbols, UniProt_ID = NA))
  
  # Create a mapping table
  uniprot_mapping <- data.frame(
    Gene = sapply(results$genes, function(x) x[[1]]$value),  # Extract gene symbol
    UniProt_ID = results$primaryAccession,
    uniProtkbId = results[["uniProtkbId"]], # UniProt ID
    
    stringsAsFactors = FALSE
  )
  
  return(uniprot_mapping)
}
