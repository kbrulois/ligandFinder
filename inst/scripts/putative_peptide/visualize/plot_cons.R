





library(tidyverse)
s_localDir <- "~/peptide_alg/build_residue_db"

cons_dat <- readRDS(fs::path(s_localDir, "aminode.rds"))





####Get weird genes
top300 <- cons_dat %>%
  mutate(avg_cons = map_dbl(cons, ~mean(as.numeric(.[["cons_og"]]), na.rm = TRUE))) %>%
  filter(!is.na(avg_cons)) %>%
  slice_min(order_by = avg_cons, n = 300) %>%
  pull(gene)

mid_cons <- cons_dat %>%
  mutate(avg_cons = map_dbl(cons, ~mean(as.numeric(.[["cons_og"]]), na.rm = TRUE))) %>%
  filter(!is.na(avg_cons) & avg_cons > 0.1 & avg_cons < 0.2) %>%
  pull(gene)

sim1 <- cons_dat %>%
  filter(species_limit == 1 & species_total > 45) %>%
  pull(gene)

sim_low <- cons_dat %>%
  filter(species_limit < 0.3) %>%
  pull(gene)

intersect(mid_cons, sim1)




####visualize conservation

genes <- c("CCL22", "CXCL12", "BRINP1", "SCGB3A2", "TLR3", "IL1B", "BUB1", "CXCL3", "CEACAM4", "CCL18", "ACTC1", "RPS14", "SEC61G", "HOXC8", "DAD1", "CALM1")
names(genes) <- c(rep("goi", 4), rep("cons_mid_species_high", 3), rep("species_low", 3), rep("cons_high_species_high", 3), rep("cons_high_species_low", 3))

genes <- c("brinp2", "cort", "cxcl14", "cxcl17", "pomc") %>% toupper
names(genes) <- rep("goi", 5)


html_file_name <- "~/AF2_analysis/conservation5.html"

plot_dir <- "~/AF2_analysis/tmp"

plot_title <- "conservation"

cons_sub <- cons_dat %>%
  filter(gene %in% !!genes)


to_plot <- cons_sub %>%
  mutate(cons = map2(cons, species_limit, \(x, y) {
    x %>%
      dplyr::rename(blos_lrs_all = cons_lrs,
                    blos_nlrs_all = cons_rs) %>%
      mutate(across(!all_of(c("AA", "index", "cons_og", "index_og")), \(z) z * y, .names = "{.col}_n"))

  })) %>%
  select(gene, cons) %>%
  unnest(cons)

to_plot <- to_plot %>%
  group_by(gene) %>%
  mutate(max = max(index))

transfrom_aligments <- function(x, vals_to = "AA") {
  x %>%
    mutate(index = row_number()) %>%
    pivot_longer(cols = -index, names_to = "metric", values_to = vals_to)
}

to_plot_aln <- cons_sub %>%
  mutate(alignment_AA = map(alignment_AA, ~transfrom_aligments(.))) %>%
  mutate(sim_mat = map(sim_mat, \(x) {

    tmp <- map(names(x), ~transfrom_aligments(x[[.]], vals_to = .))

    tmp <- purrr::reduce(tmp, function(x, y) {
      bind_cols(x, y %>% select(!any_of(names(x))))
    })


    bind_rows(
      species_dat %>%
        mutate(metric = as.character(aminode)) %>%
        mutate(index = 0) %>%
        mutate(blos = myo_sim,
               gran = myo_sim) %>%
        select(index, metric, blos, gran),
      tmp)

  })) %>%
  mutate(alignment_final = map2(alignment_AA, sim_mat, ~right_join(.x, .y, by = c("index", "metric")))) %>%
  select(gene, alignment_final) %>%
  unnest(alignment_final)



cols_to_plot <- c("cons_gran_all", "cons_gran_mam", "cons_blos_all", "cons_blos_mam", "cons_rs", "cons_lrs", "score_nn7_s8", "pathogenicity", "conservation_og")

cols_to_plot <- colnames(to_plot)
cols_to_plot <- cols_to_plot[!cols_to_plot %in% c("gene", "AA", "index", "index_og", "cons_og", "max")]

cols_to_plot <- cols_to_plot[gtools::mixedorder(cols_to_plot, decreasing = TRUE)]

cols_to_plot <- c(stringr::str_subset(cols_to_plot, "_n$", negate = TRUE),
                  stringr::str_subset(cols_to_plot, "_n$", negate = FALSE))

to_plot <- to_plot %>%
  pivot_longer(cols = any_of(cols_to_plot),
               names_to = "metric",
               values_to = "value")


to_plot <- bind_rows(to_plot,
                     to_plot_aln) %>%
  mutate(metric_type = if_else(metric %in% cols_to_plot, "", "blosum62\n-------------\ngrantham")) %>%
  mutate(metric = factor(metric, levels = c(cols_to_plot, levels(species_dat[["aminode"]]))))


library(ggplot2)
library(ggtext)

c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF", "#FD7446FF",
  "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF",
  "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF"
)



y_axis_colors <- tibble(color = c("#709AE1FF", "#FD7446FF", "#46732EFF", "#370335FF", "#8A9197FF"),
                        class = c("Mammalia", "Aves", "Lepidosauria", "Actinopteri", "Amphibia")
)

y_axis_colors <- left_join(species_dat, y_axis_colors, by = "class") %>%
  mutate(metric = as.character(aminode)) %>%
  select(color, metric)

y_axis_colors <- bind_rows(tibble(metric = cols_to_plot,
                                  color = "black"),
                           y_axis_colors)

y_axis_colors <- setNames(y_axis_colors[["color"]],
                          y_axis_colors[["metric"]])


p <- Map(\(x) {

  df <- to_plot %>% filter(gene == !!x)
  max_index <- max(df$index, na.rm = TRUE)

  ggplot2::ggplot(data = df) +
    ggplot2::geom_tile(data = df %>% filter(metric_type == ""),
                       mapping = aes(x = index, y = metric, fill = value),
                       #color = "black",
                       #lwd = 0.1
    ) +
    ggplot2::geom_tile(data = df %>% filter(metric_type != ""),
                       mapping = aes(x = index, y = metric, fill = gran),
                       #color = "black",
                       #lwd = 0.1,
                       width = 1,
                       height = 0.5,
                       position = position_nudge(y = -0.25)) +
    ggplot2::geom_tile(data = df %>% filter(metric_type != ""),
                       mapping = aes(x = index, y = metric, fill = blos),
                       #color = "black",
                       #lwd = 0.1,
                       width = 1,
                       height = 0.5,
                       position = position_nudge(y = 0.25)) +
    ggplot2::geom_text(aes(x = index, y = metric, label = AA), size = 1.8, fontface = "bold", color = "black") %>%
    ggfx::with_outer_glow(., colour = "white", sigma = 0.8, expand = 3, blend_type = "add") %>%
    ggrastr::rasterise(., dev = "ragg", dpi = 300) +
    ggh4x::facet_nested(rows = vars(metric_type), scales = "free_y", switch = "y", space = "free_y") +
    scale_fill_viridis_c(option = "H", name = "") +
    scale_y_discrete(labels = function(labs) {
      purrr::map_chr(labs, ~ glue::glue("<span style='color:{y_axis_colors[.x]}'>{.x}</span>"))
    }) +
    scale_x_continuous(expand = grid::unit(0, "lines"),
                       breaks = seq(0, max_index, by = 10),
                       minor_breaks = seq(0, max_index, by = 5)) +
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


dir.create(plot_dir)
unlink(plot_dir)

#future::plan(strategy = future::multisession(workers = 10))
future::plan(strategy = future::sequential())

options(future.globals.maxSize = Inf)

furrr::future_map(names(p), \(x) {
  message("saving ", x)
  num_res <- cons_sub %>% filter(gene == !!x) %>% pull(cons) %>% `[[`(1) %>% nrow
  svglite::svglite(filename = paste0(plot_dir, "/", x, ".svg"), width = num_res/9.35, height = 13)
  print(p[[x]])
  dev.off()
}
)


ligandFinder::html_slide_show(svg_directory = plot_dir,
                              output_file = html_file_name,
                              frames = genes,
                              categories = names(genes),
                              title = plot_title,
                              columns = 1)



























######output data for batch consurf...


uniprot <- readRDS('~/peptide_alg/uniprot.rds')

test <- uniprot %>%
  mutate(accession = paste0(">", accession)) %>%
  mutate(comb = map2(accession, sequence, c)) %>%
  {do.call(c, .[["comb"]])}

writeLines(text = test, con = "~/AF2_analysis/species_stuff/consurf_input.txt")






####Misc code...Try msa package
## read sequences
filepath <- system.file("examples", "HemoglobinAA.fasta", package="msa")
mySeqs <- Biostrings::readAAStringSet(filepath)

## perform multiple alignment
myAlignment <- msa::msa(mySeqs)

library(msa)

msa:::msaConservationScore(myAlignment, BLOSUM62)

## compute consensus scores using the BLOSUM62 matrix
data(BLOSUM62)


msaConservationScore(myAlignment, BLOSUM62)

## compute consensus scores using the BLOSUM62 matrix
## without scoring gap-gap pairs and using a different consensus sequence
msaConservationScore(myAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower")

## compute a consensus matrix first
conMat <- consensusMatrix(myAlignment)
data(PAM250)
msaConservationScore(conMat, PAM250, gapVsGap=0)






