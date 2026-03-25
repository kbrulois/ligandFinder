





metrics <- c("relASA", "mean_afm", "cons_rs", "topo")

secretome_aa %>%
  filter(accession %in% uni_id) %>%
  nest_by(accession) %>%
  pwalk(list(accession, data), .f =

function(accession, data, filename) {
  dir.create(paste0("~/Desktop/", accession))
  for(m in metrics) {
    if(TRUE) {
      dat <- data[[m]]
      dat[is.na(dat)] <- "None"
      if(!is.numeric(dat)) {
        dat <- all_desc_colors[[m]] %>% unname
      }
      filename <- paste0("~/Desktop/", accession, "/", m, ".defattr")
      data.table::fwrite(as.list(c(paste0("attribute: ", sub(pattern = "->", "_", m)),
                       paste0("match mode: 1-to-1"),
                       "recipient: residues",
                       paste0("\t", ":", 1:nrow(data), "\t", dat))),
             file = filename,
             sep = "\n")

    }
  }
}

)




to_plot <- c("pathogenicity", "relASA", "pLDDT", "conservation_og", "known_peptide", "sites")

dat <- secretome_aa %>% filter(accession == "P01042") %>% select(starts_with("CA_")) %>% as.matrix(.)

dat2 <- secretome_aa %>% filter(accession == "P01042") %>% select(all_of(to_plot)) %>% replace_na(replace = list(sites = "other"))

html_3dPlot(coordinates = dat,
            color = dat2,
            out_dir = "~/Desktop",
            file_name = "/test.html")





ps <- to_plot %>%
  filter(type == "interacting")

res_sub %>%
  mutate(map2_df(.x = chain1, .y = chain2, .f = \(x) {


  }))
  filter(accession %in% c("P01042", "O00230")) %>%
  nest_by(accession) %>%
  pmap(list(accession, data), .f =

         function(accession, data, filename) {
           dir.create(paste0("~/Desktop/", accession))
           metrics <- colnames(data)
           for(m in metrics) {
             if(TRUE) {
               dat <- data[[m]]
               dat[is.na(dat)] <- "None"
               filename <- paste0("~/Desktop/", accession, "/", m, ".defattr")
               data.table::fwrite(as.list(c(paste0("attribute: ", sub(pattern = "->", "_", m)),
                                            paste0("match mode: 1-to-1"),
                                            "recipient: residues",
                                            paste0("\t", ":", 1:nrow(data), "\t", dat))),
                                  file = filename,
                                  sep = "\n")

             }
           }
         }

  )





