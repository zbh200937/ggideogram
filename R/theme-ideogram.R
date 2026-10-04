# Standard ggplot2 theme preset for the dimensionless renderer.

#' Theme for chromosome ideograms
#'
#' `theme_ideogram()` is a small [ggplot2::theme_void()] preset. It does not
#' remove guides, set a page size, or alter coordinate scaling. A later
#' `theme()` call therefore overrides it exactly as for any other ggplot.
#' The default gap between the plot and its legends is 6 mm. Circular plots
#' can reserve additional space for chromosome text; set `legend.box.spacing`
#' in a later `theme()` call to choose an explicit gap.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param plot_margin A [ggplot2::margin()] object.
#' @param ... Additional theme settings passed to [ggplot2::theme()].
#'
#' @return A ggplot2 theme.
#' @export
theme_ideogram <- function(
    base_size = 11,
    base_family = "",
    plot_margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
    ...) {
  spacing <- grid::unit(6, "mm")
  attr(spacing, "ideogram_auto_spacing") <- TRUE
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.margin = plot_margin,
      panel.spacing = grid::unit(0, "pt"),
      legend.box.spacing = spacing
    ) +
    ggplot2::theme(...)
}
