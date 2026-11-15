# Native point and segment marker layers for the dimensionless renderer.

validate_and_filter_chr_marker <- function(data, ideogram_layout) {
  if (!nrow(data) || !all(c("chr", "position") %in% names(data))) return(data)
  chr <- validate_chr(data$chr, "chr")
  pos <- validate_coordinate(data$position, "position")
  if (!is.null(ideogram_layout)) {
    source <- ideogram_layout$data$karyotype
    i <- match(chr, source$.chr)
    if (anyNA(i)) {
      stopf("Unknown chromosome%s: %s.",
            if (sum(is.na(i)) > 1L) "s" else "",
            format_chr_rows(unique(chr[is.na(i)])))
    }
    if (any(pos < source$.start[i] | pos > source$.end[i])) {
      stopf("Genomic position is outside its chromosome bounds.")
    }
    vis <- match(chr, ideogram_layout$chrom$.chr)
    keep <- !is.na(vis) & pos >= ideogram_layout$chrom$.start[vis] &
      pos <= ideogram_layout$chrom$.end[vis]
    data <- data[keep, , drop = FALSE]
    chr <- chr[keep]
    pos <- pos[keep]
  }
  data$chr <- chr
  data$position <- pos
  data
}

StatChrMarker <- ggplot2::ggproto(
  "StatChrMarker", ggplot2::Stat,
  required_aes = c("chr", "position"),
  compute_layer = function(self, data, params, layout) {
    params$ideogram_layout <- layout$coord$layout
    data <- validate_and_filter_chr_marker(data, params$ideogram_layout)
    ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_layer(data, params, layout)
  },
  compute_panel = function(data, scales, side = "right", gap = 0.25,
                           track = NULL, component = "marker", track_position = 0.5,
                           na.rm = FALSE, ideogram_layout = NULL) {
    if (!nrow(data)) return(data.frame())
    data <- validate_and_filter_chr_marker(data, ideogram_layout)
    if (!nrow(data)) return(data.frame())
    chr <- data$chr
    pos <- data$position
    data$ideogram_chr <- chr
    data$ideogram_anchor_position <- pos
    data$ideogram_position <- pos
    data$ideogram_side <- side
    data$ideogram_gap <- gap
    if (!is.null(track)) {
      data$ideogram_marker_track <- track
      data$ideogram_track_position <- track_position
    }
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

#' @noRd
geom_chr_marker <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrMarker,
    position = "identity",
    ...,
    side = c("right", "left", "inner", "outer"),
    gap = 0.25,
    track = NULL,
    track_position = 0.5,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE,
    .geom = ggplot2::GeomPoint) {
  side <- canonical_chr_side(match.arg(side))
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  check_track_position(track_position)
  layer <- ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = .geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      side = side, gap = gap, track = track, track_position = track_position,
      component = "marker", na.rm = na.rm
    ), list(...))
  )
  layer$position <- locus_position(layer$position)
  ggplot2::ggproto("LayerChrMarker", layer,
    ideogram_marker_spec = list(side = side, gap = gap, track = track))
}

#' @noRd
geom_chr_link <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrMarker,
    position = "identity",
    ...,
    side = c("right", "left", "inner", "outer"),
    gap = 0.25,
    track = NULL,
    track_position = 0.5,
    na.rm = FALSE,
    show.legend = FALSE,
    inherit.aes = FALSE) {
  side <- canonical_chr_side(match.arg(side))
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  check_track_position(track_position)
  layer <- ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomSegment,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      side = side, gap = gap, track = track, track_position = track_position,
      component = "link", na.rm = na.rm
    ), list(...))
  )
  layer$position <- circular_segment_position(layer$position)
  layer
}

validate_optional_marker_track <- function(track) {
  if (!is.null(track) &&
      (!is.character(track) || length(track) != 1L ||
       is.na(track) || !nzchar(track))) {
    stopf("`track` must be NULL or one non-empty track identifier.")
  }
  invisible(track)
}

check_track_position <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stopf("`track_position` must be one number from 0 (near edge) to 1 (far edge).")
  }
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
#' chromosome and side. geom_locus(geom = "link") makes displacement explicit.
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
  if (extent <= layout$marker_extent[[spec$side]]) return(plot)
  base <- layout$base_spec
  base$tracks <- base$tracks %||% track_layout()
  space <- layout$marker_extent
  space[[spec$side]] <- extent
  attr(base$tracks, "marker_extent") <- space
  rebuild_ideogram(plot, base)
}
