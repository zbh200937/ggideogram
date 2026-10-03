# Standard text and segment geoms in chromosome coordinates.

StatChrInterval <- ggplot2::ggproto(
  "StatChrInterval", ggplot2::Stat,
  required_aes = c("chr", "start", "end"),
  compute_layer = function(self, data, params, layout) {
    params$ideogram_layout <- layout$coord$layout
    ggplot2::ggproto_parent(ggplot2::Stat, self)$compute_layer(data, params, layout)
  },
  compute_panel = function(
      data, scales, placement = "overlay",
      side = "right", gap = 0.25, track = NULL, track_position = 0.5,
      na.rm = FALSE, ideogram_layout = NULL) {
    first <- validate_coordinate(data$start, "start")
    last <- validate_coordinate(data$end, "end")
    if (any(first < 1 | first > last | first != floor(first) | last != floor(last))) {
      stopf("Chromosome intervals require positive integer `start <= end` in one-based closed coordinates.")
    }
    chr <- validate_chr(data$chr, "chr")
    source <- ideogram_layout$data$karyotype
    i <- match(chr, source$.chr)
    if (anyNA(i) || any(first - 1 < source$.start[i] | last > source$.end[i])) {
      stopf("Chromosome intervals are outside the source chromosome bounds.")
    }
    i <- match(chr, ideogram_layout$chrom$.chr)
    low <- pmax(first - 1, ideogram_layout$chrom$.start[i])
    high <- pmin(last, ideogram_layout$chrom$.end[i])
    keep <- !is.na(i) & low < high
    if (!any(keep)) return(data.frame())
    data <- data[keep, , drop = FALSE]
    data$ideogram_chr <- chr[keep]
    data$ideogram_position <- low[keep]
    data$ideogram_end_position <- high[keep]
    data$ideogram_interval <- TRUE
    data$ideogram_interval_placement <- placement
    if (placement == "beside") {
      data$ideogram_side <- side
      data$ideogram_gap <- gap
      if (!is.null(track)) {
        data$ideogram_marker_track <- track
        data$ideogram_track_position <- track_position
      }
    }
    data$x <- 0
    data$y <- 0
    data$xend <- 0
    data$yend <- 0
    data
  }
)

#' @noRd
geom_chr_interval <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrInterval,
    position = "identity",
    ...,
    placement = c("overlay", "beside"),
    side = c("right", "left", "inner", "outer"),
    gap = 0.25,
    track = NULL,
    track_position = 0.5,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE) {
  placement <- match.arg(placement)
  side <- canonical_chr_side(match.arg(side))
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  check_track_position(track_position)
  if (!is.null(track) && placement != "beside") {
    stopf("`track` requires `placement = \"beside\"`.")
  }
  layer <- ggplot2::layer(
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
      track = track, track_position = track_position,
      na.rm = na.rm
    ), list(...))
  )
  layer$position <- circular_segment_position(layer$position)
  layer
}

#' @noRd
geom_chr_region <- function(...) geom_chr_interval(...)

#' @noRd
geom_chr_text <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrMarker,
    position = "identity",
    ...,
    side = c("right", "left", "inner", "outer"),
    gap = 0.4,
    track = NULL,
    track_position = 0.5,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE) {
  side <- canonical_chr_side(match.arg(side))
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  check_track_position(track_position)
  layer <- ggplot2::layer(
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
      track = track, track_position = track_position,
      component = "marker",
      na.rm = na.rm
    ), list(...))
  )
  layer$position <- locus_position(layer$position)
  layer
}
