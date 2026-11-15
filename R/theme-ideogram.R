# Standard ggplot2 theme preset for the dimensionless renderer.

#' Theme for chromosome ideograms
#'
#' `theme_ideogram()` is a small [ggplot2::theme_void()] preset. It does not
#' remove guides, set a page size, or alter coordinate scaling. A later
#' `theme()` call therefore overrides it exactly as for any other ggplot.
#' The default gap between the plot and its legends is 6 mm. Plots can reserve
#' additional space for text below the panel; set `legend.box.spacing`
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
  if (missing(plot_margin)) {
    attr(plot_margin, "ideogram_auto_margin") <- TRUE
  }
  auto_margin <- function(t = 0, r = 0, b = 0, l = 0) {
    m <- ggplot2::margin(t, r, b, l)
    attr(m, "ideogram_auto_margin") <- TRUE
    m
  }
  spacing <- grid::unit(6, "mm")
  attr(spacing, "ideogram_auto_spacing") <- TRUE
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.margin = plot_margin,
      panel.spacing = grid::unit(0, "pt"),
      legend.box.spacing = spacing,
      plot.title = ggplot2::element_text(hjust = 0, margin = auto_margin(t = 5.5)),
      plot.subtitle = ggplot2::element_text(hjust = 0, margin = auto_margin(t = 5.5)),
      plot.caption = ggplot2::element_text(hjust = 1, margin = auto_margin(t = 5.5)),
      axis.title.x = ggplot2::element_text(vjust = 1, margin = auto_margin(t = 2.75)),
      axis.title.x.top = ggplot2::element_text(vjust = 0, margin = auto_margin(b = 2.75)),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 1, margin = auto_margin(r = 2.75)),
      axis.title.y.right = ggplot2::element_text(angle = -90, vjust = 1, margin = auto_margin(l = 2.75))
    ) +
    ggplot2::theme(...)
}
