


library(bio3d)
s_localDir <- "~/peptide_alg/build_residue_db"

pdb_dir <- "~/peptide_alg/UP000005640_9606_HUMAN_v4/"

list.files(pdb_dir) %>% head

accession <- secretome %>% filter(gene == "CXCL14") %>% pull(accession)
features <- secretome %>% filter(gene == "CXCL14") %>% pull(features) %>% `[[`(1)

accession <- secretome[["accession"]][1]
features <- secretome[["features"]][1]


get_afdsb <- function(features, accession) {

pdb_file <- fs::path(pdb_dir, paste0("AF-", accession, "-F1-model_v4.pdb"))

if(!file.exists(pdb_file)) {
  return(features)
}
pdb <- bio3d::read.pdb(pdb_file)

# extract sulfur atoms from cysteines
sg <- subset(pdb$atom, resid == "CYS" & elety == "SG")

# distance matrix
d <- as.matrix(dist(cbind(sg$x, sg$y, sg$z)))

# find bonded pairs (exclude self)
hits <- which(d < 3 & d > 1.6, arr.ind = TRUE)

# unique pairs only
hits <- hits[hits[,1] < hits[,2], , drop = FALSE]

disulfides <- tibble(
  type = "afdsb",
  evidence = "afdb",
  start = sg$resno[hits[,1]],
  end = sg$resno[hits[,2]],
  source = "afdb",
  description = as.character(d[hits])
)

if(nrow(disulfides) > 0) {
features <- bind_rows(disulfides, features)
}

return(features)
}


start <- Sys.time()
#future::plan(strategy = future::multisession(workers = 6))
future::plan(strategy = future::sequential())


test <- furrr::future_map2(.x = secretome[["features"]],
                            .y = secretome[["accession"]],
                            get_afdsb)

end <- Sys.time()

end - start


secretome$features <- test





map_lgl(test, \(x) {
  "afdsb" %in% x$type
}) %>% sum



