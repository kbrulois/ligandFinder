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

cons_dir <- paste0(s_localDir,"/processed/chunked_cons")

dir.create(cons_dir)

chunks <- 1:nrow(entries)

chunks2 <- split.default(chunks, cut(chunks, breaks = c(seq(0, length(chunks), by = 1000), max(chunks))))

for(c in seq_along(chunks2)) {
  message("processing ", c)
  x <- chunks2[[c]]
  if(x[length(x)] == last_entry_og) {
    inds <- c(which(stringr::str_detect(aminode, paste0(entries[["gene"]][x[1]], " "))),
              length(aminode))
    writeLines(aminode[inds[1]:inds[2]], paste0(cons_dir, "/", c, ".txt"))
  } else {
    inds <- c(which(stringr::str_detect(aminode, paste0(entries[["gene"]][x[1]], " "))),
              which(stringr::str_detect(aminode, entries[["gene"]][x[length(x)] + 1])))
    inds <- c(inds[1], inds[2] - 1)
    writeLines(aminode[inds[1]:inds[2]], paste0(cons_dir, "/", c, ".txt"))
  }
}




length(unique(entries[["gene"]])) == last_entry

future::plan(strategy = future::multisession(workers = 10))

cons_dat <- furrr::future_map(.x = list.files(cons_dir), .f = \(cf) {

  aminode <- readLines(paste0(cons_dir, "/", cf))

  entry_inds <- stringr::str_detect(aminode, "^>")

  entries <- aminode[entry_inds]

  entries <- dplyr::bind_rows(lapply(entries, \(x) {
    x <- strsplit(x, " ")[[1]]
    tibble(gene = x[1],
           ensemble = x[2])
  }))

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
      tibble(AA = z[1],
             index = z[2],
             cons = z[3])
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

saveRDS(cons_dat, paste0(s_localDir,"/processed/aminode.rds"))

