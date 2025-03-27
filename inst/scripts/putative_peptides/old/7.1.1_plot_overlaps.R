

x_cutoff <- 80
y_cutoff <- 80

x_var <- "percent_ol_gtp"
y_var <- "percent_ol_phs_hsr"


fetch_goi_dif <- \(x) {
  to_ret <- secretome_roi[["dterm"]][secretome_roi[["roi_name"]] == x]
  if(length(to_ret) == 0) {
    to_ret <- NA
  }
  if(length(to_ret) > 1) {
    to_ret <- max(to_ret)
  }
  return(to_ret)
}

to_plot <- secretome_roi %>%
              filter(!is.na(!!sym(x_var)) & 
                       roi_type %in% c("phs_hsr", "phs_hsr_N", "phs_hsr_C") & 
                       roi_length > 3) %>%
              rowwise() %>%
              mutate(dterm = fetch_goi_dif(overlap_region)) %>%
              mutate(across(c(x_var, y_var), .fns = ~jitter(., factor = 3)))






to_plot2 <- to_plot %>%
  mutate(
    quadrant = case_when(
      !!sym(x_var) >= x_cutoff & !!sym(y_var) >= y_cutoff ~ "good",
      !!sym(x_var) < x_cutoff & !!sym(y_var) >= y_cutoff ~ "too small",
      !!sym(x_var) < x_cutoff & !!sym(y_var) < y_cutoff ~ "too big",
      !!sym(x_var) >= x_cutoff & !!sym(y_var) < y_cutoff ~ "bad"
    )
  ) %>%
  group_by(quadrant) %>%
  summarise(count = n(),
         percent = count/ nrow(to_plot) * 100)

to_plot3 <- tibble(xmin = c(0, x_cutoff, x_cutoff, 0), 
       xmax = c(x_cutoff, 100, 100, x_cutoff),
       ymin = c(0, y_cutoff, 0, y_cutoff),
       ymax = c(y_cutoff, 100, y_cutoff, 100))

fill_scale <- c(bad = "#8A9197FF", good = "#FD7446FF", `too big` = "#FED439FF", `too small` = "#709AE1FF")

to_plot2 <- bind_cols(to_plot2, to_plot3)

fill_vars <- c("roi_length", "score_nn8_s8_entire","dterm")



plot_dir <- "~/peptide_alg/plots_tmp12"
dir.create(plot_dir)

all_p <- lapply(fill_vars, \(x) {

p <- ggplot2::ggplot(data = to_plot2) +
  ggplot2::geom_rect(aes(xmin = xmin,
                         xmax = xmax,
                         ymin = ymin,
                         ymax = ymax,
                         fill = quadrant)) +
  scale_fill_manual(values = fill_scale, labels = paste0(to_plot2[["quadrant"]], " (", 
                                                         to_plot2[["count"]], "; ", 
                                                         round(to_plot2[["percent"]], 0), "%)")) +
  ggnewscale::new_scale_fill() +
  ggplot2::geom_point(data = to_plot,
                       mapping = aes(x = percent_ol_gtp, 
                                     y = percent_ol_phs_hsr,
                                     fill = !!sym(x)),
                       size = 3,
                       pch = 21) +
  ggplot2::geom_hline(yintercept = 80) +
  ggplot2::geom_vline(xintercept = 80) +
  scale_fill_viridis_c(option = "H") + #name = "N verus C score\ndifference\nof corresponding\nknown peptide") +
  ggtitle("Putative Peptides that overlap with known peptides") +
  xlab("Overlap (Percent of Known Peptide)") +
  ylab("Overlap (Percent of Putative Peptide)") +
  #scale_shape_manual(values = c(phs_hsr = 21, phs_hsr_N = 22, phs_hsr_C = 23)) +
  ggplot2::theme_bw()
  
ggsave(filename = paste0(plot_dir, "/", x, ".svg"),
       plot = p, svglite::svglite, width = 6, height = 5)

})





files <- list.files(plot_dir)

files <- sub(".svg$", "", files)

html_slide_show(svg_directory = plot_dir,
                output_file = paste0("~/peptide_alg/", "peptide_overlaps2", ".html"),
                frames = files,
                categories = NULL,
                title = "peptide_eval",
                columns = 1)









