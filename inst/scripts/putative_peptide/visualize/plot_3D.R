
residue_anno %>%
  pmap(list(gene, data = c("cons", "dssp", "af_missense", "af_xyz", "msa_mean")), .f = 
         
         function(gene, data, filename) {
           dir.create(paste0("~/Desktop/", gene))
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

to_plot <- c("pathogenicity", "relASA", "pLDDT", "conservation_og", "known_peptide", "sites")

dat <- secretome_aa %>% filter(accession == "P01042") %>% select(starts_with("CA_")) %>% as.matrix(.)

dat2 <- secretome_aa %>% filter(accession == "P01042") %>% select(all_of(to_plot)) %>% replace_na(replace = list(sites = "other"))

html_3dPlot(coordinates = dat, 
            color = dat2, 
            out_dir = "~/Desktop", 
            file_name = "/test.html")



