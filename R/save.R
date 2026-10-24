# Standard ggplot saving with explicit output dimensions.

#' Save an ideogram with ggplot2
#'
#' This convenience wrapper delegates plotting and validation to
#' [ggplot2::ggsave()]. Width and height are ordinary device choices.
#' Omitting both uses ggplot2's own
#' device defaults.
#'
#' @param filename Output filename.
#' @param plot A [ggideogram()] plot or any object supported by
#'   [ggplot2::ggsave()], including grobs and compositions with a
#'   [grid::grid.draw()] method.
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
