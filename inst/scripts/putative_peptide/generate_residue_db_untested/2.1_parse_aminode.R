####parse aminode (~1h)

##download aminode: data https://drive.google.com/file/d/1RGGxmud7YL-ipBA6IF9-03TEZH9_0sx3/view

aminode <- readLines(paste0(s_localDir, "/raw/Aminode Bulk Data/Bulk_AlignedSequences.txt"))

entry_inds <- stringr::str_detect(aminode, "^>")

entries <- aminode[entry_inds]

entries <- dplyr::bind_rows(lapply(entries, \(x) {
  x <- strsplit(x, " ")[[1]]
  tibble(gene = x[1],
         ensemble = x[2])
}))

last_entry_og <- nrow(entries)

cons_dir <- paste0(s_localDir,"/processed/chunked_cons_species")

dir.create(cons_dir)

chunks <- 1:nrow(entries)

chunks2 <- split.default(chunks, cut(chunks, breaks = c(seq(0, length(chunks), by = 1000), max(chunks))))

for(c in seq_along(chunks2)) {
  message("processing ", c)
  x <- chunks2[[c]]
  if(x[length(x)] == last_entry_og) {
    inds <- c(which(stringr::str_detect(aminode, paste0(entries[["gene"]][x[1]]))),
              length(aminode))
    writeLines(aminode[inds[1]:inds[2]], paste0(cons_dir, "/", c, ".txt"))
  } else {
    inds <- c(which(stringr::str_detect(aminode, paste0(entries[["gene"]][x[1]]))),
              which(stringr::str_detect(aminode, entries[["gene"]][x[length(x)] + 1])))
    inds <- c(inds[1], inds[2] - 1)
    writeLines(aminode[inds[1]:inds[2]], paste0(cons_dir, "/", c, ".txt"))
  }
}




length(unique(entries[["gene"]])) == last_entry_og

start <- Sys.time()

future::plan(strategy = future::multisession(workers = 10))

cons_dat <- furrr::future_map(.x = list.files(cons_dir), .f = \(cf) {

  aminode <- readLines(paste0(cons_dir, "/", cf))

  entry_inds <- stringr::str_detect(aminode, "^>")

  entries <- aminode[entry_inds]

  entries <- tibble(gene = entries)

  last_entry <- nrow(entries)


  entries[["cons"]] <- lapply(1:nrow(entries), \(x) {
    if(x == last_entry) {
      inds <- c(which(stringr::str_detect(aminode, entries[["gene"]][x])) + 1,
                length(aminode))
      return(aminode[inds[1]:inds[2]])
    } else {
      inds <- c(which(stringr::str_detect(aminode, entries[["gene"]][x])),
                which(stringr::str_detect(aminode, entries[["gene"]][x + 1])))
      inds <- c(inds[1] + 1, inds[2] - 1)
      return(aminode[inds[1]:inds[2]])
    }
  })

  clean_cons <- function(x) {
    bind_rows(lapply(x, \(y) {
      z <- strsplit(y, " ")[[1]]
      ens <- strsplit(z[2], "\\|")[[1]]

      if(length(ens) == 3) {
        ens <- tibble(ensg = ens[1],
                      enst = ens[2],
                      ensp = ens[3])
      }


      bind_cols(tibble(species = z[1]),
                ens,
                sequence = z[3])
    }
    ))
  }

  entries %>%
    dplyr::mutate(gene = sub("^>", "", gene)) %>%
    dplyr::mutate(cons = map(.x = cons, .f = clean_cons))

})

end <- Sys.time()

end - start

cons_dat <- bind_rows(cons_dat)

saveRDS(cons_dat, paste0(s_localDir,"/processed/aminode_seq.rds"))




cons_dat <- readRDS(paste0(s_localDir,"/processed/aminode_seq.rds"))



all_species <- lapply(cons_dat[[2]], `[[`, "species") %>% do.call(c, .) %>% unique(.) %>% stringr::str_subset(., "^>", negate = TRUE)

test <- cons_dat %>%
  mutate(num_species = map_int(cons, ~nrow(.))) %>%
  filter(num_species < 5)








####assemble species data

species_dat <- tibble(aminode = lapply(cons_dat[[2]], `[[`, "species") %>%
                        do.call(c, .) %>%
                        unique(.) %>%
                        stringr::str_subset(., "^>", negate = TRUE)
)

species_dat <- species_dat %>%
                mutate(timetree_in = case_when(aminode == "Canis_familiaris" ~ "Canis_lupus",
                                               aminode == "Mustela_putorius_furo" ~ "Mustela_putorius",
                                               TRUE ~ aminode)) %>%
                arrange(aminode)


data.table::fwrite(list(test = species_dat[["timetree_in"]]), "~/AF2_analysis/species_stuff/all_species.txt", sep = "\n", col.names = FALSE)



tree <- ape::read.tree("~/AF2_analysis/species_stuff/all_species.nwk")
dist_mat <- as.matrix(cophenetic(tree))

evol_dat <- tibble(timetree_out = rownames(dist_mat),
                   timetree_in = rownames(dist_mat),
                   myo = dist_mat["Homo_sapiens", ],
                   rank = rank(myo, ties.method = "first")) %>%
            arrange(timetree_out)

setdiff(evol_dat$timetree_in, species_dat$timetree_in)

evol_dat[["timetree_in"]][evol_dat[["timetree_in"]] == "Carlito_syrichta"] <- "Tarsius_syrichta"
evol_dat[["timetree_in"]][evol_dat[["timetree_in"]] == "Notamacropus_eugenii"] <- "Macropus_eugenii"

setdiff(evol_dat$timetree_in, species_dat$timetree_in)

species_dat <- left_join(species_dat, evol_dat, by = "timetree_in")

library(taxize)

classifications <- classification(species_dat$aminode, db = "ncbi")

broad_group <- sapply(classifications, function(x) {
  if (!is.null(x)) {
    val <- x$name[x$rank == "class"]
    if (length(val) > 0) val else NA
  } else NA
})

broad_group[is.na(broad_group)] <- c("Actinopteri", "Lepidosauria")

class_dat <- tibble(aminode = names(broad_group),
                    class = unname(broad_group))

species_dat <- left_join(species_dat, class_dat, by = "aminode")

species_dat <- species_dat %>%
                mutate(aminode = factor(aminode, levels = .[["aminode"]][order(.[["rank"]], decreasing = TRUE)]))


p <- ggplot2::ggplot(data = species_dat) +
  ggplot2::geom_point(aes(x = myo, y = aminode, fill = class), pch = 21, stroke = 0.5, size = 4) +
  ggplot2::scale_fill_discrete(palette = ggsci::pal_simpsons(), name = "") +
  ylab("") +
  ggplot2::theme_bw()

svglite::svglite(filename = "~/AF2_analysis/species_stuff/all_species_aminode.svg", width = 7, height = 11)
p
dev.off()

data(BLOSUM62)

library(Biostrings)
BLOSUM62






uniprot <- readRDS('~/peptide_alg/uniprot.rds')

test <- uniprot %>%
  mutate(accession = paste0(">", accession)) %>%
  mutate(comb = map2(accession, sequence, c)) %>%
  {do.call(c, .[["comb"]])}

writeLines(text = test, con = "~/AF2_analysis/species_stuff/consurf_input.txt")












