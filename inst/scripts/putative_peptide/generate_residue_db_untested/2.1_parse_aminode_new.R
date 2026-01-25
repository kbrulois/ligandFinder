####parse aminode msa (~1h)

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
    inds <- c(which(stringr::str_detect(aminode, paste0(entries[["gene"]][x[1]])))[1],
              which(stringr::str_detect(aminode, entries[["gene"]][x[length(x)] + 1]))[1])
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







###combine aminode conservation and msa data

library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"

cons_dat <- readRDS("~/peptide_alg/aminode.rds")

cons_dat_s <- readRDS(paste0(s_localDir,"/processed/aminode_seq.rds"))

cons_dat <- cons_dat %>%
            select(-ensemble)

cons_dat_s <- cons_dat_s %>%
                dplyr::rename(alignment = cons)

cons_dat <- left_join(cons_dat, cons_dat_s, by = "gene")



####assemble species data. skip...already done...

if(FALSE) {

species_dat <- tibble(aminode = lapply(cons_dat[["alignment"]], `[[`, "species") %>%
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



tree <- ape::read.tree(system.file("extdata/all_species.nwk", package = "ligandFinder"))
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
  mutate(aminode = factor(aminode, levels = .[["aminode"]][order(.[["rank"]], decreasing = TRUE)])) %>%
  mutate(myo_sim = myo/max(myo)) %>%
  mutate(myo_dissim = 1 - myo_sim)


p <- ggplot2::ggplot(data = species_dat) +
  ggplot2::geom_point(aes(x = myo, y = aminode, fill = class), pch = 21, stroke = 0.5, size = 4) +
  ggplot2::scale_fill_discrete(palette = ggsci::pal_simpsons(), name = "") +
  ylab("") +
  ggplot2::theme_bw()

svglite::svglite(filename = "~/AF2_analysis/species_stuff/all_species_aminode.svg", width = 7, height = 11)
p
dev.off()

}

#####do custom conservation

species_dat <- readRDS(system.file("extdata/species_dat.rds", package = "ligandFinder"))

g_mat <- grantham::grantham_distances_matrix

rownames(g_mat) <- grantham::as_one_letter(rownames(g_mat))
colnames(g_mat) <- grantham::as_one_letter(colnames(g_mat))

g_mat <- abs(g_mat/max(g_mat) - 1)

g_mat <- rbind(g_mat, X = rep(NA, 20))
g_mat <- cbind(g_mat, X = c(rep(NA, 20), 1))

g_mat <- rbind(g_mat, `-` = rep(NA, 21))
g_mat <- cbind(g_mat, `-` = c(rep(NA, 21), 1))


library(msa)
data(BLOSUM62)

blos_mat <- BLOSUM62

rownames(blos_mat)[25] <- colnames(blos_mat)[25] <- "-"

blos_mat <- apply(blos_mat, 2, \(x) scales::rescale(x, to = c(0, 1)))

blos_mat["X",] <- NA
blos_mat[ , "X"] <- NA
blos_mat["-",] <- NA
blos_mat[ , "-"] <- NA

sim_mats <- list(blos = blos_mat,
                 gran = g_mat)



cons_dat <- cons_dat %>%
  mutate(alignment = map(alignment, ~left_join(., species_dat, by = join_by(species == aminode)) %>% arrange(myo_sim)))


cons_dat <- cons_dat %>%
  mutate(species_total = map_int(alignment, ~nrow(.))) %>%
  mutate(species_limit = map_dbl(alignment, ~max(.[["myo_sim"]])))


cons_dat[["cons"]] <- cons_dat %>%
  select(gene, cons) %>%
  unnest(cons) %>%
  mutate(cons_og = as.numeric(cons)) %>%
  select(-cons) %>%
  mutate(cons_rs = scales::rescale(cons_og, c(1,0))) %>%
  mutate(cons_lrs = scales::rescale(log(cons_og + 0.1), c(1,0))) %>%
  mutate(index = as.numeric(index)) %>%
  mutate(index_og = index) %>%
  mutate(index = row_number(), .by = gene) %>%
  nest(.by = gene) %>%
  pull(data)

cons_dat <- cons_dat %>%
  filter(gene %in% !!genes)

cons_dat <- cons_dat %>%
  mutate(alignment_AA = map(alignment, \(x) {

    tmp <- map(x[["sequence"]], \(y) {
      stringr::str_split(y, "")[[1]]
    }) %>%
      do.call(cbind, .) %>%
      as_tibble
    colnames(tmp) <- x[["species"]]
    tmp
  }))

cons_dat <- cons_dat %>%
  mutate(sim_mat = map(alignment_AA, \(x) {
    ref_seq <- x[["Homo_sapiens"]]
    Map(\(sm) {
      x %>%
        mutate(across(everything(), \(y) {
          sm[cbind(ref_seq, y)]
        }))
    }, sim_mats)

  }))

cons_dat <- cons_dat %>%
  mutate(cons = map2(cons, sim_mat,  do_cons_hu_ref))

cons_dat <- cons_dat %>%
  dplyr::slice(-17102) %>%
  mutate(cons = map2(cons, alignment_AA, do_cons))

saveRDS(cons_dat, fs::path(s_localDir, "processed/aminode.rds"))





























