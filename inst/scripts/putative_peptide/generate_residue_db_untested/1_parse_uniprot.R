
###parse uniprot (~4.5h)

s_localDir <- "~/peptide_alg/build_residue_db"

dir.create(paste0(s_localDir, "/processed"))

uniprot_file <- paste0(s_localDir, "/raw/UP000005640_9606.xml")

uniprot <- xml2::read_xml(x = uniprot_file)

entries <- xml2::xml_find_all(uniprot, "//d1:entry")

uniprot_processed <- paste0(s_localDir,"/processed/chunked_xml")

dir.create(uniprot_processed)

chunks <- 1:length(entries)

chunks2 <- split.default(chunks, cut(chunks, breaks = c(seq(0, length(chunks), by = 1000), max(chunks))))

for(c in seq_along(chunks2)) {

  test <- xml2::xml_new_root(.value = xml2::read_xml(paste0("<uniprot></uniprot>")))

  node <- xml2::xml_find_first(test, "//uniprot")

  for(i in chunks2[[c]]) xml2::xml_add_child(node, entries[[i]])

  xml2::write_xml(test, file = paste0(uniprot_processed, "/", c, ".xml"))

}

uni_dir <- paste0(s_localDir, "/chunked_tibble")

dir.create(uni_dir)

start <- Sys.time()

for(c in seq_along(chunks2)) {

  future::plan(strategy = future::multisession(workers = 10))

  test <- dplyr::bind_rows(
    furrr::future_map(seq_along(chunks2[[c]]), \(x) {

      message("processing ", x)

      uniprot <- xml2::read_xml(x = paste0(uniprot_processed, "/", c, ".xml"))
      entries <- xml2::xml_find_all(uniprot, "//d1:entry")
      x <- entries[x]

      gene_name <- xml2::xml_children(xml2::xml_find_first(x = x, "d1:gene"))

      meta_dat <- list()

      comment_nodes <- xml2::xml_find_all(x = x, "d1:comment")

      meta_dat[["comment"]] <- bind_rows(lapply(comment_nodes, \(y) {

        c1 <- xml2::xml_children(y)
        c2 <- xml2::xml_children(c1)

        if(length(c2) > 0) {
        bind_rows(lapply(c1, \(a) {
          c2 <- xml2::xml_children(a)
          bind_rows(lapply(c2, \(b) {
            tibble(
                   annotation_name = xml2::xml_name(y),
                   annotation_type = xml2::xml_attr(y, "type"),
                   annotation_id = xml2::xml_attr(a, "id"),
                   name_1 = xml2::xml_name(a),
                   name_2 = xml2::xml_name(b),
                   annotation = xml2::xml_text(b))

         }))
         }))
        } else {
          bind_rows(lapply(c1, \(a) {
              tibble(annotation_name = xml2::xml_name(y),
                     annotation_type = xml2::xml_attr(y, "type"),
                     annotation_id = xml2::xml_attr(a, "id"),
                     name_1 = xml2::xml_name(a),
                     name_2 = NA,
                     annotation = xml2::xml_text(a))

          }))
        }


        }))





      dbRef_nodes <- xml2::xml_find_all(x = x, "d1:dbReference")

      meta_dat[["dbRef"]] <- bind_rows(lapply(dbRef_nodes, \(y) {

        dat <- xml2::xml_find_all(y, "d1:property")

        tibble(annotation_name = xml2::xml_name(y),
               annotation_type = xml2::xml_attr(y, "type"),
               annotation_id = xml2::xml_attr(y, "id"),
               name_1 = xml2::xml_attr(dat, "type"),
               name_2 = NA,
               annotation = xml2::xml_attr(dat, "value"))

      }))

      annotations <- bind_rows(meta_dat)




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
                      end = xml2::xml_attr(xml2::xml_find_first(x = position_nodes, "d1:end"), "position"),
                      description = xml2::xml_attr(y, attr = "description")
)
             }))),
             annotations = list(annotations))

    })
  )

  saveRDS(test, paste0(uni_dir, "/", c, ".rds"))

  gc()

}

end <- Sys.time()

end - start

uniprot_t <- dplyr::bind_rows(
  map(list.files(uni_dir), ~readRDS(paste0(uni_dir, "/", .)))
  )

saveRDS(uniprot_t, paste0(s_localDir, "/processed/uniprot.rds"))

to_save <- uniprot_t %>%
  select(-annotations, -full_name, -sequence) %>%
  unnest(features)


data.table::fwrite(to_save, "~/Desktop/uniprot_features.csv")

to_save <- uniprot_t %>%
  select(-features, -full_name, -sequence) %>%
  unnest(annotations)


data.table::fwrite(to_save, "~/Desktop/uniprot_annotations.csv")


uniprot_t %>%
  filter(gene %in% c("GPR15", "MARCH4F")) %>%
  select(-annotations, -full_name, -sequence) %>%
  unnest(features) %>%
  View


