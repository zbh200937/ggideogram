# Standard save boundary; no page, pixel or DPI compatibility renderer.

#' Save an ideogram with ggplot2
#'
#' This convenience wrapper validates the plot and delegates to
#' [ggplot2::ggsave()]. Width and height are ordinary device choices and are
#' never inferred from a historical canvas. Omitting both uses ggplot2's own
#' device defaults.
#'
#' @param filename Output filename.
#' @param plot A [ggideogram()] plot or any standard ggplot composition.
#' @param width,height,units,dpi Standard [ggplot2::ggsave()] arguments.
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return `filename`, invisibly.
#' @export
save_ideogram <- function(
    filename,
    plot = ggplot2::last_plot(),
    width = NULL,
    height = NULL,
    units = "in",
    dpi = 300,
    ...) {
  if (!inherits(plot, "ggplot")) {
    stopf("`plot` must be a ggplot object or compatible ggplot composition.")
  }
  arguments <- list(
    filename = filename,
    plot = plot,
    units = units,
    dpi = dpi,
    ...
  )
  if (!is.null(width)) arguments$width <- width
  if (!is.null(height)) arguments$height <- height
  do.call(ggplot2::ggsave, arguments)
  invisible(filename)
}
