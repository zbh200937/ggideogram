#' Select a genomic window while retaining its source karyotype
#'
#' A view preserves the full source data and changes only the visible genomic
#' range. Axes show original base-pair positions. Truncated chromosome ends
#' are flat; real chromosome ends remain rounded. Use [view_chr_data()] to
#' select observations for marker, interval and track layers before adding them.
#'
#' @param data Karyotype data frame or `ideogram_data`.
#' @param chr One chromosome identifier.
#' @param start,end Window boundaries in the original base-pair coordinates.
#'
#' @return An `ideogram_data` object carrying an explicit `view` specification.
#' @examples
#' k <- data.frame(Chr = "Chr1", Start = 0, End = 1000)
#' view <- chr_view(k, "Chr1", 200, 500)
#' ggideogram(view, orientation = "horizontal", axis = TRUE)
#' @export
chr_view <- function(data, chr, start, end) {
  semantic <- as_ideogram_data(data)
  chr <- validate_chr(chr, "chr")
  if (length(chr) != 1L || !chr %in% semantic$karyotype$.chr) {
    stopf("`chr` must select one chromosome in the source karyotype.")
  }
  start <- validate_coordinate(start, "start")
  end <- validate_coordinate(end, "end")
  if (length(start) != 1L || length(end) != 1L || start >= end) {
    stopf("A view requires one interval with start < end.")
  }
  source <- semantic$karyotype[semantic$karyotype$.chr == chr, ]
  if (start < source$.start || end > source$.end) {
    stopf("The view is outside the source chromosome bounds.")
  }
  semantic$view <- list(chr = chr, start = start, end = end)
  semantic
}

#' Select points or clip intervals to a chromosome view
#'
#' Point observations are selected by their original position. One-based closed
#' intervals are intersected with the view's continuous boundaries and retain
#' their source bounds in `.source_start` and `.source_end`. For fractional
#' view boundaries, the returned closed interval includes every intersecting
#' base; chromosome layers clip its geometry to the exact view boundary.
#' No interpolation or aggregation is performed: line tracks
#' join only the selected observations. Data on other chromosomes are excluded.
#'
#' @param data Data frame of observations.
#' @param view An object returned by [chr_view()].
#' @param chr Column containing chromosome identifiers.
#' @param position Point-coordinate column, or `NULL` for interval data.
#' @param start,end One-based closed interval-coordinate columns when
#'   `position = NULL`.
#'
#' @return A data frame suitable for existing chromosome layers. The input
#'   object is not modified.
#' @export
view_chr_data <- function(data, view, chr = "Chr", position = NULL,
                          start = "Start", end = "End") {
  check_projection_data(data)
  if (!inherits(view, "ideogram_data") || is.null(view$view)) {
    stopf("`view` must be created by chr_view().")
  }
  window <- view$view
  ids <- validate_chr(projection_column(data, chr, "chr"), "chr")
  result <- data[ids == window$chr, , drop = FALSE]
  source <- view$karyotype[view$karyotype$.chr == window$chr, ]
  if (!is.null(position)) {
    p <- validate_coordinate(projection_column(result, position, "position"), "position")
    if (any(p < source$.start | p > source$.end)) {
      stopf("Point observations are outside the source chromosome bounds.")
    }
    result <- result[p >= window$start & p <= window$end, , drop = FALSE]
  } else {
    lo <- validate_coordinate(projection_column(result, start, "start"), "start")
    hi <- validate_coordinate(projection_column(result, end, "end"), "end")
    if (any(lo < 1 | lo > hi | lo != floor(lo) | hi != floor(hi) |
            lo - 1 < source$.start | hi > source$.end)) {
      stopf("Intervals require closed integer coordinates inside the source chromosome bounds.")
    }
    keep <- hi > window$start & lo - 1 < window$end
    if (!".source_start" %in% names(result)) result$.source_start <- lo
    if (!".source_end" %in% names(result)) result$.source_end <- hi
    result[[start]] <- pmax(lo, floor(window$start) + 1)
    result[[end]] <- pmin(hi, ceiling(window$end))
    result <- result[keep, , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}
