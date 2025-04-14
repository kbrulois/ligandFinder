
###parse uniprot (~2h)

uniprot_file <- "~/peptide_alg/UP000000589_10090.xml"

uniprot <- xml2::read_xml(x = uniprot_file)

entries <- xml2::xml_find_all(uniprot, "//d1:entry")


dir.create("~/peptide_alg/chunked_xml_ms")

chunks <- 1:length(entries)

chunks2 <- split.default(chunks, cut(chunks, breaks = c(seq(0, length(chunks), by = 1000), max(chunks))))

for(c in seq_along(chunks2)) {
  
  test <- xml2::xml_new_root(.value = xml2::read_xml(paste0("<uniprot></uniprot>")))
  
  node <- xml2::xml_find_first(test, "//uniprot")
  
  for(i in chunks2[[c]]) xml2::xml_add_child(node, entries[[i]])
  
  xml2::write_xml(test, file = paste0("~/peptide_alg/chunked_xml_ms/", c, ".xml"))
  
}

uni_dir <- "~/peptide_alg/chunked_tibble_ms"

dir.create(uni_dir)

start <- Sys.time()

for(c in seq_along(chunks2)) {
  
  future::plan(strategy = future::multisession(workers = 10))
  
  test <- dplyr::bind_rows(
    furrr::future_map(seq_along(chunks2[[c]]), \(x) {
      message("processing ", x)
      uniprot <- xml2::read_xml(x = paste0("~/peptide_alg/chunked_xml_ms/", c, ".xml"))
      entries <- xml2::xml_find_all(uniprot, "//d1:entry")
      x <- entries[x]
      gene_name <- xml2::xml_children(xml2::xml_find_first(x = x, "d1:gene"))
      tibble(accession = xml2::xml_text(xml2::xml_find_first(x = x, "d1:accession")),
             gene = xml2::xml_text(gene_name[1]),
             full_name = xml2::xml_text(xml2::xml_find_first(x = x, "d1:protein")),
             sequence = xml2::xml_text(xml2::xml_find_first(x = x, "d1:sequence")),
             features = list(dplyr::bind_rows(lapply(xml2::xml_find_all(x = x, "d1:feature"), \(y) {
               
               position_nodes <- xml2::xml_children(y)
               if("position" %in% xml2::xml_name(xml2::xml_children(position_nodes))) {
                 start <- xml2::xml_attr(xml2::xml_find_first(x = position_nodes, "d1:position"), "position")
               } else {
                 start <- xml2::xml_attr(xml2::xml_find_first(x = position_nodes, "d1:begin"), "position")
               }
               
               
               tibble(type = xml2::xml_attr(y, attr = "type"),
                      evidence = paste0(xml2::xml_attr(y, attr = "evidence"), collapse = "_"),
                      start = start,
                      end = xml2::xml_attr(xml2::xml_find_first(x = position_nodes, "d1:end"), "position"))
             }))))
      
    })
  )
  
  saveRDS(test, paste0("~/peptide_alg/chunked_tibble_ms/", c, ".rds"))
  
  gc()
  
}

end <- Sys.time()

end - start

uniprot_t <- dplyr::bind_rows(
  lapply(list.files(uni_dir), \(x) {
    readRDS(paste0(uni_dir, "/", x))
  }
  ))


saveRDS(uniprot_t, "~/peptide_alg/uniprot_ms.rds")

