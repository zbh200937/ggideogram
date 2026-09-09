# Native point and segment marker layers for the dimensionless renderer.

StatChrMarker <- ggplot2::ggproto(
  "StatChrMarker", ggplot2::Stat,
  required_aes = c("chr", "position"),
  compute_panel = function(data, scales, side = "right", gap = 0.25,
                           track = NULL, component = "marker",
                           na.rm = FALSE) {
    data$ideogram_chr <- as.character(data$chr)
    data$ideogram_anchor_position <- data$position
    data$ideogram_position <- data$position
    data$ideogram_side <- side
    data$ideogram_gap <- gap
    if (!is.null(track)) data$ideogram_marker_track <- track
    data$ideogram_link <- identical(component, "link")
    data$x <- 0
    data$y <- 0
    if (identical(component, "link")) {
      data$xend <- 0
      data$yend <- 0
    }
    data
  }
)

#' Draw native ggplot2 markers beside chromosomes
#'
#' `geom_chr_marker()` is a thin chromosome-coordinate adapter around
#' [ggplot2::GeomPoint]. Point glyphs, physical size, scales and guides are all
#' supplied by ggplot2; the package does not construct marker polygons or
#' convert canvas pixels.
#'
#' Map `chr` and `position` to the chromosome identifier and genomic point.
#' All ordinary point aesthetics remain available, including `shape`, `size`,
#' `colour`, `fill`, `alpha` and `stroke`.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()]. It
#'   must contain `chr` and `position` unless inherited from the plot.
#' @param data Marker data frame. When `NULL`, the plot data is used.
#' @param stat Statistical transformation. The default creates the semantic
#'   chromosome coordinates consumed by the ideogram coordinate system.
#' @param position Position adjustment. Use [position_chr_repel()] for bounded
#'   one-dimensional avoidance along each chromosome.
#' @param ... Other arguments passed to [ggplot2::layer()], including standard
#'   point aesthetics set to constants.
#' @param side Side of the chromosome body, `"left"` or `"right"`. For a
#'   horizontal ideogram these names follow the chromosome's local transverse
#'   axis rather than page left/right.
#' @param gap Distance from the chromosome edge to the point centre, in
#'   chromosome-body-width units. Without a declared track, the marker reserves
#'   a symmetric lane extending the same `gap` beyond its centre; bp axes are
#'   placed outside this lane. Use a wider track or gap for large point sizes.
#' @param track Optional identifier declared by [track_layout()]. When supplied,
#'   the point is centred in that track and `side`/`gap` only remain fallback
#'   values for plots without a marker track.
#' @param na.rm If `FALSE`, missing required aesthetics are removed with a
#'   warning.
#' @param show.legend Should this layer be included in standard ggplot2 guides?
#' @param inherit.aes If `FALSE`, the default, aesthetics are not inherited
#'   from the base plot.
#'
#' @return A ggplot2 layer whose geom inherits directly from
#'   [ggplot2::GeomPoint].
#' @examples
#' markers <- data.frame(
#'   Chr = c("Chr1", "Chr1", "Chr2"),
#'   Pos = c(20, 60, 40),
#'   Type = c("circle", "box", "circle")
#' )
#' chromosomes <- data.frame(
#'   Chr = c("Chr1", "Chr2"), Start = 0, End = 100
#' )
#' ggideogram(chromosomes) +
#'   geom_chr_marker(
#'     data = markers,
#'     ggplot2::aes(chr = Chr, position = Pos, shape = Type, fill = Type),
#'     size = 2
#'   ) +
#'   ggplot2::scale_shape_manual(values = c(circle = 21, box = 22))
#' @export
geom_chr_marker <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrMarker,
    position = "identity",
    ...,
    side = c("right", "left"),
    gap = 0.25,
    track = NULL,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE) {
  side <- match.arg(side)
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  layer <- ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      side = side, gap = gap, track = track,
      component = "marker", na.rm = na.rm
    ), list(...))
  )
  ggplot2::ggproto("LayerChrMarker", layer,
    ideogram_marker_spec = list(side = side, gap = gap, track = track))
}

#' Draw native leader links from chromosomes to markers
#'
#' `geom_chr_link()` uses [ggplot2::GeomSegment]. Its start remains attached to
#' the true genomic position at the chromosome edge; its end follows the marker
#' display position, including an optional [position_chr_repel()] adjustment.
#'
#' @inheritParams geom_chr_marker
#' @param show.legend Should this segment layer contribute to ggplot2 guides?
#'   The default is `FALSE` so a matching point layer owns the marker key.
#'
#' @return A ggplot2 layer whose geom inherits directly from
#'   [ggplot2::GeomSegment].
#' @examples
#' chromosomes <- data.frame(Chr = "Chr1", Start = 0, End = 100)
#' markers <- data.frame(Chr = "Chr1", Pos = c(49, 50, 51))
#' repel <- position_chr_repel(8, units = "bp")
#' ggideogram(chromosomes) +
#'   geom_chr_link(
#'     data = markers, ggplot2::aes(chr = Chr, position = Pos),
#'     position = repel
#'   ) +
#'   geom_chr_marker(
#'     data = markers, ggplot2::aes(chr = Chr, position = Pos),
#'     position = repel, size = 1.5
#'   )
#' @export
geom_chr_link <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrMarker,
    position = "identity",
    ...,
    side = c("right", "left"),
    gap = 0.25,
    track = NULL,
    na.rm = FALSE,
    show.legend = FALSE,
    inherit.aes = FALSE) {
  side <- match.arg(side)
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomSegment,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      side = side, gap = gap, track = track,
      component = "link", na.rm = na.rm
    ), list(...))
  )
}

validate_optional_marker_track <- function(track) {
  if (!is.null(track) &&
      (!is.character(track) || length(track) != 1L ||
       is.na(track) || !nzchar(track))) {
    stopf("`track` must be NULL or one non-empty track identifier.")
  }
  invisible(track)
}

PositionChrRepel <- ggplot2::ggproto(
  "PositionChrRepel", ggplot2::Position,
  required_aes = c(
    "ideogram_chr", "ideogram_position", "ideogram_anchor_position"
  ),
  setup_params = function(self, data) {
    list(min_distance = self$min_distance, units = self$units)
  },
  compute_panel = function(data, params, scales) {
    required <- c(
      "ideogram_chr", "ideogram_position", "ideogram_anchor_position"
    )
    missing <- setdiff(required, names(data))
    if (length(missing)) {
      stopf(paste0(
        "`position_chr_repel()` requires a chromosome-semantic layer; ",
        "missing %s."), paste0("`", missing, "`", collapse = ", "))
    }
    data$ideogram_repel_distance <- params$min_distance
    data$ideogram_repel_units <- params$units
    data
  }
)

#' Repel marker display positions along each chromosome
#'
#' This position adjustment keeps `position` as the true genomic anchor and
#' solves a bounded one-dimensional spacing problem independently for every
#' chromosome and side. [geom_chr_link()] makes any displacement explicit.
#'
#' The requested separation is a data/layout distance, not a conversion from
#' point millimetres: the physical point size is only known when the plot is
#' drawn. If all points cannot fit at the requested distance, they are spread
#' evenly over the chromosome while their true anchors remain unchanged.
#'
#' @param min_distance Positive minimum separation.
#' @param units `"bp"` interprets `min_distance` as genomic base pairs;
#'   `"body_width"` interprets it in normalized chromosome-width units.
#'
#' @return A chromosome-aware ggplot2 Position object.
#' @export
position_chr_repel <- function(min_distance,
                               units = c("bp", "body_width")) {
  check_positive_layout(min_distance, "min_distance")
  units <- match.arg(units)
  ggplot2::ggproto(
    NULL, PositionChrRepel,
    min_distance = min_distance,
    units = units
  )
}

#' @method ggplot_add LayerChrMarker
#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.LayerChrMarker <- function(object, plot, ...) {
  plot <- NextMethod()
  layout <- plot$coordinates$layout
  if (!inherits(layout, "ideogram_layout_v2")) return(plot)
  spec <- object$ideogram_marker_spec
  if (!is.null(spec$track)) return(plot)
  layout$marker_extent <- layout$marker_extent %||% c(left = 0, right = 0)
  # An untracked marker is centred in a symmetric lane: `gap` from the
  # body edge to its centre and the same distance beyond the centre.
  extent <- 2 * spec$gap
  layout$marker_extent[[spec$side]] <- max(
    layout$marker_extent[[spec$side]], extent)
  point <- offset_chr_points(layout, layout$chrom$.chr,
    list(x = layout$chrom$.axis_start_x, y = layout$chrom$.axis_start_y),
    side = spec$side, distance = layout$chromosome_width / 2 + extent)
  layout$bounds$x <- range(layout$bounds$x, point$x)
  layout$bounds$y <- range(layout$bounds$y, point$y)
  plot <- update_plot_ideogram_layout(plot, layout)
  axis_layers <- list()
  for (index in seq_along(plot$layers)) {
    layer <- plot$layers[[index]]
    if (is.null(layer$ideogram_axis_spec)) next
    replacement <- do.call(chromosome_axis_layers,
      c(list(layout = layout), layer$ideogram_axis_spec))
    plot$layers[[index]] <- replacement[[layer$ideogram_axis_part]]
    axis_layers <- c(axis_layers, list(plot$layers[[index]]))
  }
  reserve_chromosome_axis_space(plot, axis_layers)
}
