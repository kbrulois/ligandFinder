## Shared chart tokens for the window-model figures.
##
## Values are the dataviz reference palette's light-surface instance. Only
## categorical slots 1-3 are documented as clearing the all-pairs separation
## floors, and dot forms (beeswarm, scatter) are held to the all-pairs list, so
## these figures never use more than three series.
##
## Labels in these figures must be ASCII: the png device renders "+/-", the
## middle dot and en/em dashes as ".." and "...".

LF_VIZ <- list(
  surface = "#fcfcfb",   # chart surface
  ink     = "#0b0b0b",   # text-primary
  ink2    = "#52514e",   # text-secondary
  grid    = "#e6e5e1",   # recessive grid
  s1      = "#2a78d6",   # categorical slot 1, blue
  s2      = "#eb6834",   # categorical slot 2, orange
  s3      = "#1baf7a"    # categorical slot 3, aqua
)

#' Minimal theme for the window-model figures
#'
#' Recessive grid on the y only (these are dot plots over a categorical x),
#' text in ink tokens rather than series colour.
#' @param base_size passed to [ggplot2::theme_minimal()].
#' @export
lf_viz_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background    = ggplot2::element_rect(fill = LF_VIZ$surface, colour = NA),
      panel.background   = ggplot2::element_rect(fill = LF_VIZ$surface, colour = NA),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = LF_VIZ$grid, linewidth = 0.3),
      axis.text          = ggplot2::element_text(colour = LF_VIZ$ink2),
      axis.title         = ggplot2::element_text(colour = LF_VIZ$ink2),
      plot.title         = ggplot2::element_text(colour = LF_VIZ$ink, face = "bold", size = 12),
      plot.subtitle      = ggplot2::element_text(colour = LF_VIZ$ink2, size = 9.5, lineheight = 1.15),
      plot.caption       = ggplot2::element_text(colour = LF_VIZ$ink2, size = 8, hjust = 0))
}
