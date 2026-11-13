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


library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"

cons_dat <- readRDS(paste0(s_localDir,"/processed/aminode_seq.rds"))





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












og_cons_dat <- readRDS("~/peptide_alg/aminode.rds")


g_mat <- grantham::grantham_distances_matrix

rownames(g_mat) <- grantham::as_one_letter(rownames(g_mat))
colnames(g_mat) <- grantham::as_one_letter(colnames(g_mat))

g_mat <- abs(g_mat/max(g_mat) - 1)

######addd X and dash

g_mat <- rbind(g_mat, X = rep(NA, 20))
g_mat <- cbind(g_mat, X = c(rep(NA, 20), 1))

g_mat <- rbind(g_mat, `-` = rep(NA, 21))
g_mat <- cbind(g_mat, `-` = c(rep(NA, 21), 1))

rescale_g_mat <- function(x, inflection_point = 0.85) {
  ifelse(x >= inflection_point, x,  x - 1)
}

g_mat_0.85 <- rescale_g_mat(g_mat)



## read sequences
filepath <- system.file("examples", "HemoglobinAA.fasta", package="msa")
mySeqs <- Biostrings::readAAStringSet(filepath)

## perform multiple alignment
myAlignment <- msa::msa(mySeqs)
library(msa)



## compute consensus scores using the BLOSUM62 matrix
data(BLOSUM62)

blos_mat <- BLOSUM62

rownames(blos_mat)[25] <- colnames(blos_mat)[25] <- "-"

blos_mat <- apply(blos_mat, 2, \(x) scales::rescale(x, to = c(0, 1)))

blos_mat["X",] <- NA
blos_mat[ , "X"] <- NA
blos_mat["-",] <- NA
blos_mat[ , "-"] <- NA

msaConservationScore(myAlignment, BLOSUM62)

## compute consensus scores using the BLOSUM62 matrix
## without scoring gap-gap pairs and using a different consensus sequence
msaConservationScore(myAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower")

## compute a consensus matrix first
conMat <- consensusMatrix(myAlignment)
data(PAM250)
msaConservationScore(conMat, PAM250, gapVsGap=0)





aln <- x %>%
  select(matches("V\\d")) %>%
  as.matrix

consensus <- Biostrings::consensusString(aln)
rel_sub <- colMeans(mat != strsplit(consensus, "")[[1]])

msaConservationScore.matrix <- function(x, substitutionMatrix, gapVsGap=NULL, …) {
  if (!is.matrix(substitutionMatrix) ||
      … ) stop("substitution matrix is not in proper format")
  …
  out <- drop(apply(x, 2, function(y) crossprod(y, substitutionMatrix %*% y)))
  names(out) <- unlist(strsplit(msaConsensusSequence(x, …), “”))
  out
}




top300 <- og_cons_dat %>%
            mutate(avg_cons = map_dbl(cons, ~mean(as.numeric(.[["cons"]]), na.rm = TRUE))) %>%
            filter(!is.na(avg_cons)) %>%
            slice_min(order_by = avg_cons, n = 300) %>%
            pull(gene)

sim1 <- cons_dat %>%
          filter(species_limit == 1 & species_total > 35) %>%
          pull(gene)

sim_low <- cons_dat %>%
          filter(species_limit < 0.3) %>%
          pull(gene)

intersect(top300, sim1)



cons_dat <- cons_dat %>%
              mutate(cons = map(cons, ~left_join(., species_dat, by = join_by(species == aminode)) %>% arrange(myo_sim)))


cons_dat <- cons_dat %>%
  mutate(species_total = map_int(cons, ~nrow(.))) %>%
  mutate(species_limit = map_dbl(cons, ~max(.[["myo_sim"]])))









genes <- c("CCL22", "CXCL12", "BRINP1", "CXCL3", "CEACAM4", "CCL18", "ACTC1", "RPS14", "SEC61G", "HOXC8", "DAD1", "CALM1")

names(genes) <- c(rep("goi", 3), rep("species_low", 3), rep("cons_high_species_high", 3), rep("cons_high_species_low", 3))


cons_test <- cons_dat %>%
              filter(gene %in% !!genes)


cons_test <- cons_test %>%
  mutate(cons = map(cons, \(x) {

    tmp <- map(x[["sequence"]], \(y) {
      stringr::str_split(y, "")[[1]]
    }) %>%
      do.call(rbind, .) %>%
      as_tibble

    bind_cols(x, tmp)
  }))


cons_test <- cons_test %>%
  mutate(cons_gran = map(cons, \(x) {

    ref <- x %>%
            filter(species == "Homo_sapiens") %>%
            select(matches("V\\d"))

    x %>%
      mutate(across(matches("V\\d"), \(y) {
        g_mat[cbind(y, rep(ref[[cur_column()]], length(y)))]
      })) %>%
      filter(species != "Homo_sapiens")

  })) %>%
  mutate(cons_blos = map(cons, \(x) {

    ref <- x %>%
      filter(species == "Homo_sapiens") %>%
      select(matches("V\\d"))

    cons_input <- x %>%
      mutate(across(matches("V\\d"), \(y) {
        blos_mat[cbind(y, rep(ref[[cur_column()]], length(y)))]
      })) %>%
      filter(species != "Homo_sapiens")

  })) %>%
  mutate(cons_score = pmap(list(cons, cons_gran, cons_blos), \(a, b, c) {

    ref <- a %>%
      filter(species == "Homo_sapiens") %>%
      select(matches("V\\d"))

    cons_gran_all <- b %>%
      summarise(across(matches("V\\d"), \(y) {
        wts <- case_when(is.na(y) ~ NA,
                         y < 0.85 ~ .data[["myo_dissim"]],
                         y >= 0.85 ~ .data[["myo_sim"]])
        weighted.mean(x = y, w = wts, na.rm = TRUE)
      })) %>%
      unlist

    cons_gran_mam <- b %>%
      filter(class == "Mammalia") %>%
      summarise(across(matches("V\\d"), \(y) {
        wts <- case_when(is.na(y) ~ NA,
                         y < 0.85 ~ .data[["myo_dissim"]],
                         y >= 0.85 ~ .data[["myo_sim"]])
        weighted.mean(x = y, w = wts, na.rm = TRUE)
      })) %>%
      unlist

    cons_blos_all <- c %>%
      summarise(across(matches("V\\d"), \(y) {
        wts <- case_when(is.na(y) ~ NA,
                         y < 0.85 ~ .data[["myo_dissim"]],
                         y >= 0.85 ~ .data[["myo_sim"]])
        weighted.mean(x = y, w = wts, na.rm = TRUE)
      })) %>%
      unlist

    cons_blos_mam <- c %>%
      filter(class == "Mammalia") %>%
      summarise(across(matches("V\\d"), \(y) {
        wts <- case_when(is.na(y) ~ NA,
                         y < 0.85 ~ .data[["myo_dissim"]],
                         y >= 0.85 ~ .data[["myo_sim"]])
        weighted.mean(x = y, w = wts, na.rm = TRUE)
      })) %>%
      unlist

    tibble(AA = unlist(ref),
           index = 1:ncol(ref),
           cons_gran_all = cons_gran_all,
           cons_gran_mam = cons_gran_mam,
           cons_blos_all = cons_blos_all,
           cons_blos_mam = cons_blos_mam)
  }))



og_cons_test <- og_cons_dat %>%
  unnest(cons) %>%
  mutate(cons_og = as.numeric(cons)) %>%
  select(-cons) %>%
  mutate(cons_rs = scales::rescale(cons_og, c(1,0))) %>%
  mutate(cons_lrs = scales::rescale(log(cons_og + 0.1), c(1,0))) %>%
  filter(gene %in% !!genes) %>%
  mutate(index = as.numeric(index))

cons_test1 <- cons_test %>%
  select(-cons, -cons_gran, -cons_blos) %>%
  unnest(cons_score)

cons_test1 <- left_join(cons_test1, og_cons_test, by = c("gene", "AA", "index"))

cons_test1 <- cons_test1 %>%
  group_by(gene) %>%
  mutate(max = max(index))

cons_test2 <- cons_test %>%
          mutate(alignment = pmap(list(cons, cons_gran, cons_blos), \(x, y, z) {

bind_cols(
x %>%
  filter(species != "Homo_sapiens") %>%
  mutate(V0 = rep("", max(row_number())), .before = "V1") %>%
  pivot_longer(cols = matches("V\\d+"), names_to = "index", values_to = "AA") %>%
  mutate(index = stringr::str_remove(index, "^V") %>% as.numeric),

  y %>%
  dplyr::rename(V0 = myo_sim) %>%
  pivot_longer(cols = matches("V\\d+"), names_to = "index", values_to = "value") %>%
  mutate(index = stringr::str_remove(index, "^V") %>% as.numeric) %>%
  select(value),

  z %>%
  dplyr::rename(V0 = myo_sim) %>%
  pivot_longer(cols = matches("V\\d+"), names_to = "index", values_to = "value2") %>%
  mutate(index = stringr::str_remove(index, "^V") %>% as.numeric) %>%
  select(value2)
)

          })) %>%
  select(gene, alignment) %>%
  unnest(alignment)  %>%
  dplyr::rename(metric = species)



cols_to_plot <- c("cons_gran_all", "cons_gran_mam", "cons_blos_all", "cons_blos_mam", "cons_rs", "cons_lrs", "score_nn7_s8", "pathogenicity", "conservation_og")


to_plot <- cons_test1 %>%
  pivot_longer(cols = any_of(cols_to_plot),
               names_to = "metric",
               values_to = "value")

to_plot <- bind_rows(to_plot,
                     cons_test2) %>%
            mutate(metric_type = if_else(metric %in% cols_to_plot, "", "blosum62\n--------\ngrantham")) %>%
  mutate(metric = factor(metric, levels = c(cols_to_plot, levels(species_dat[["aminode"]])))) %>%
  mutate(value2 = if_else(metric_type == "blosum62\n--------\ngrantham", value2, value))



library(ggplot2)
library(ggtext)

c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF",
  "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF",
  "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF"
)



y_axis_colors <- c("#709AE1FF", "#FD7446FF", "#46732EFF", "#370335FF", "#8A9197FF")

names(y_axis_colors) <- c("Mammalia", "Aves", "Lepidosauria", "Actinopteri", "Amphibia"
)

y_axis_colors <- to_plot %>%
  mutate(y_axis_colors = if_else(class %in% names(y_axis_colors), y_axis_colors[class], "black")) %>%
  pull(y_axis_colors, metric)

y_axis_colors <- y_axis_colors[which(!duplicated(names(y_axis_colors)))]


p <- Map(\(x) {

ggplot2::ggplot(to_plot %>% filter(gene == !!x)) +
  ggplot2::geom_tile(mapping = aes(x = index, y = metric, fill = value),
                     #color = "black",
                     #lwd = 0.1,
                     width = 1,
                     height = 0.5,
                     position = position_nudge(y = -0.25)) +
  ggplot2::geom_tile(mapping = aes(x = index, y = metric, fill = value2),
                     #color = "black",
                     #lwd = 0.1,
                     width = 1,
                     height = 0.5,
                     position = position_nudge(y = 0.25)) +
  ggplot2::geom_text(aes(x = index, y = metric, label = AA), size = 1.8, fontface = "bold", color = "black") %>%
  ggfx::with_outer_glow(., colour = "white", sigma = 0.8, expand = 3, blend_type = "add") %>%
  ggrastr::rasterise(., dev = "ragg", dpi = 400) +
  ggh4x::facet_nested(rows = vars(metric_type), scales = "free_y", switch = "y", space = "free_y") +
  scale_fill_viridis_c(option = "H", name = "") +
  scale_y_discrete(labels = function(labs) {
      purrr::map_chr(labs, ~ glue::glue("<span style='color:{y_axis_colors[.x]}'>{.x}</span>"))
    }) +
  scale_x_continuous(expand = grid::unit(0, "lines"),
                     breaks = seq(0, max_index[1], by = 10),
                     minor_breaks = seq(0, max_index[1], by = 5)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = ggtext::element_markdown(size = 9),
        panel.spacing = grid::unit(0, "lines"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.y.left = element_text(face = 2, angle = 0, hjust = 1, size = 5),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.justification = "left",
        legend.ticks = element_line(color = "black", linewidth = 0.1),
        legend.key.height = unit(0.2, "cm"))
}, unique(to_plot[["gene"]]))


plot_dir <- "~/AF2_analysis/tmp"
dir.create(plot_dir)
unlink(plot_dir)

future::plan(strategy = future::sequential())

furrr::future_map(names(p), \(x) {
       num_res <- max(grep("V\\d+", colnames(cons_test %>% filter(gene == !!x) %>% pull(cons) %>% `[[`(1))))
       svglite::svglite(filename = paste0(plot_dir, "/", x, ".svg"), width = num_res/9.35, height = 6)
       print(p[[x]])
       dev.off()
}
)




ligandFinder::html_slide_show(svg_directory = plot_dir,
                output_file = "~/AF2_analysis/conservation.html",
                frames = genes,
                categories = names(genes),
                title = "conservation",
                columns = 1)












htmlwidgets::saveWidget(plotly::ggplotly(p[[1]]),
                        "~/Desktop/interactive_plot.html",
                        selfcontained = TRUE)









uniprot <- readRDS('~/peptide_alg/uniprot.rds')

test <- uniprot %>%
  mutate(accession = paste0(">", accession)) %>%
  mutate(comb = map2(accession, sequence, c)) %>%
  {do.call(c, .[["comb"]])}

writeLines(text = test, con = "~/AF2_analysis/species_stuff/consurf_input.txt")












