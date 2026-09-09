# Standard text and segment geoms in chromosome coordinates.

StatChrInterval <- ggplot2::ggproto(
  "StatChrInterval", ggplot2::Stat,
  required_aes = c("chr", "start", "end"),
  compute_panel = function(
      data, scales, placement = "overlay",
      side = "right", gap = 0.25, track = NULL,
      na.rm = FALSE) {
    if (any(data$start > data$end, na.rm = TRUE)) {
      stopf("Chromosome intervals require `start <= end`.")
    }
    data$ideogram_chr <- as.character(data$chr)
    data$ideogram_position <- data$start
    data$ideogram_end_position <- data$end
    data$ideogram_interval <- TRUE
    data$ideogram_interval_placement <- placement
    if (placement == "beside") {
      data$ideogram_side <- side
      data$ideogram_gap <- gap
      if (!is.null(track)) data$ideogram_marker_track <- track
    }
    data$x <- 0
    data$y <- 0
    data$xend <- 0
    data$yend <- 0
    data
  }
)

#' Draw chromosome intervals with a standard segment geom
#'
#' `geom_chr_interval()` maps genomic `start`/`end` to the chromosome long
#' axis and delegates drawing, physical line width, aesthetics, scales and
#' guides to [ggplot2::GeomSegment]. Use `placement = "beside"` for an
#' interval track or the default `"overlay"` for a line on the body axis.
#'
#' @param mapping Aesthetic mapping with `chr`, `start`, and `end`.
#' @param data Interval data frame.
#' @param stat Chromosome interval Stat.
#' @param position Standard ggplot2 position adjustment.
#' @param ... Standard segment parameters and constant aesthetics.
#' @param placement Draw on the body axis or beside it.
#' @param side,gap Beside placement in chromosome-width units.
#' @param track Optional declared track used for beside placement.
#' @param na.rm,show.legend,inherit.aes Standard layer arguments.
#'
#' @return A ggplot2 layer using `GeomSegment`.
#' @export
geom_chr_interval <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrInterval,
    position = "identity",
    ...,
    placement = c("overlay", "beside"),
    side = c("right", "left"),
    gap = 0.25,
    track = NULL,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE) {
  placement <- match.arg(placement)
  side <- match.arg(side)
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  if (!is.null(track) && placement != "beside") {
    stopf("`track` requires `placement = \"beside\"`.")
  }
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomSegment,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      placement = placement,
      side = side,
      gap = gap,
      track = track,
      na.rm = na.rm
    ), list(...))
  )
}

#' @rdname geom_chr_interval
#' @export
geom_chr_region <- function(...) geom_chr_interval(...)

#' Draw text at chromosome positions with the standard text geom
#'
#' This is the text counterpart of [geom_chr_marker()]. It changes only the
#' chromosome coordinates; text size, family, alignment, aesthetics, scales
#' and guides remain ordinary ggplot2 behaviour. Add [geom_chr_link()] when a
#' displaced label should retain an explicit locus leader.
#'
#' @inheritParams geom_chr_marker
#' @param gap Distance from the chromosome edge to the text anchor, in
#'   chromosome-body-width units.
#'
#' @return A ggplot2 layer using `GeomText`.
#' @export
geom_chr_text <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrMarker,
    position = "identity",
    ...,
    side = c("right", "left"),
    gap = 0.4,
    track = NULL,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE) {
  side <- match.arg(side)
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomText,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      side = side,
      gap = gap,
      track = track,
      component = "marker",
      na.rm = na.rm
    ), list(...))
  )
}
