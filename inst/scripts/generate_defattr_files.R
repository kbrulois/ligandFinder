


expand_by_residue <- function(x,
                              dat_to_expand = c("ref_index",
                                                "topo",
                                                "features_expanded",
                                                "modification",
                                                "domain",
                                                "SV",
                                                "DBC",
                                                "dssp",
                                                "af_missense",
                                                "cons",
                                                "alignment_AA",
                                                "aa_scores")
) {

  x <- x %>%
    mutate(features_expanded = map2(features, sequence_uni, expand_features)) %>%
    mutate(modification = map2(features, sequence_uni, assemble_res_mod)) %>%
    mutate(domain = map2(features, sequence_uni, assemble_dom)) %>%
    mutate(SV = map2(features, sequence_uni, assemble_SV)) %>%
    mutate(DBC = map2(features, sequence_uni, assemble_DBC)) %>%
    mutate(ref_index = map(sequence_uni, \(x) {tibble(index = 1:nchar(x),
                                                      AA = stringr::str_split(x, "", simplify = TRUE) %>% `c`)}))

  x <- x %>%
    mutate(topo = map2(sequence_uni, topo, \(x, y) tibble(AA = stringr::str_split(x, "", simplify = TRUE) %>% c,
                                                          topo = stringr::str_split(y, "", simplify = TRUE) %>% c))) %>%
    mutate(to_expand = pmap(pick(any_of(dat_to_expand)),
                            bind_cols, .name_repair = "minimal")) %>%
    mutate(to_expand = map(to_expand, \(x) x[, !duplicated(colnames(x))])) %>%
    mutate(to_expand = map(to_expand, \(x) x[, colnames(x) != ""]))

  to_return <- x %>%
    select(accession, to_expand) %>%
    unnest(to_expand)

  return(to_return)

}


.libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1')
library(dplyr)
library(purrr)
library(tidyr)
remotes::install_github("kbrulois/ligandFinder")
library(ligandFinder)

secretome <- readRDS("/oak/stanford/groups/ebutcher/kevin/secretome_latest.rds")



mets <- list(cons = c("blos_wt_all_n", "cons_rs", "blos_wt_mam", "blos_wt_all", "gran_wt_all"),
             af_missense = c("mean_afm", "min_afm"),
             dssp = c("relASA"),
             aa_scores = c("pep_xgb4c", "chem_xgb3c", "pep_nn4c", "chem_nn4c")
)

mets$aa_scores <- c(mets$aa_scores, paste0(mets$aa_scores, "_s6"))

all_mets <- do.call(`c`, mets) %>% unname

names(all_mets) <- rep(names(mets), sapply(mets, length))

dir_path <- "/oak/stanford/groups/ebutcher/kevin/defattr"

dir.create(dir_path)

input <- secretome %>%
  group_split(accession)

num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) %/% 1.5

#num_cores <- 60
future::plan(strategy = future::multicore(workers = num_cores))


furrr::future_walk(input, .f = function(x) {
  try({
            accession <- unique(x[["accession"]])

            dat <- expand_by_residue(x)
            dir.create(paste0(dir_path, "/", accession))
            for(m in all_mets) {
              if(TRUE) {
                dat2 <- dat[[m]]
                dat2[is.na(dat2)] <- "None"
                filename <- paste0(dir_path, "/", accession, "/", m, ".defattr")
                data.table::fwrite(as.list(c(paste0("attribute: ", sub(pattern = "->", "_", m)),
                                             paste0("match mode: 1-to-1"),
                                             "recipient: residues",
                                             paste0("\t", ":", 1:nrow(dat), "\t", dat2))),
                                   file = filename,
                                   sep = "\n")

              }
            }
  })
          }

  )


