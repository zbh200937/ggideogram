# Start-position counts from positional GFF columns.

#' Count features per window from a GFF file
#'
#' Produces a `Chr` / `Start` / `End` / `Value` data frame suitable for
#' [geom_track()] with ordinary ggplot2 line or column layers.
#'
#' @param input Path to a GFF file, or a data frame of one already read (columns
#'   are used positionally, as in GFF: 1 = sequence name, 3 = feature type,
#'   4 = start, 5 = end when present). Selected features must have positive
#'   integer, 1-based closed coordinates inside the source chromosome bounds.
#'   Unknown chromosomes and invalid coordinates are errors.
#' @param karyotype Path to a tab-separated karyotype file with a header, or a
#'   data frame with `Chr` and `End`, or an [as_ideogram_data()] result.
#'   Optional `Start` is the 0-based left
#'   boundary and defaults to zero.
#' @param feature Feature type to count, matched against column 3 of the GFF.
#' @param window Positive integer window size in base pairs.
#' @return A data frame with columns `Chr`, `Start`, `End`, `Value`, in
#'   karyotype order. Output intervals are 1-based closed. Each feature is
#'   counted once by its start position, in the window with `Start <= p <= End`.
#'   Windows begin at each chromosome's `Start + 1`; the last window is
#'   truncated to its `End`. Chromosomes without selected features retain
#'   zero counts. The `window_summary` attribute follows [bin_genome()] and
#'   records method `"count"` and `count_position = "start"`.
#' @export
GFFex <- function(input, karyotype, feature = "gene", window = 1000000) {
  gff <- if (is.data.frame(input)) {
    input
  } else {
    utils::read.table(input, stringsAsFactors = FALSE, header = FALSE,
                      comment.char = "#", sep = "\t", quote = "")
  }
  if (ncol(gff) < 4) stopf("`input` needs at least 4 GFF columns, got %d.", ncol(gff))

  if (!is.character(feature) || length(feature) != 1L || is.na(feature) || !nzchar(feature)) {
    stopf("`feature` must be one non-empty feature type.")
  }
  keep <- which(as.character(gff[[3]]) == feature)
  positions <- data.frame(Chr = as.character(gff[[1]][keep]), Start = gff[[4]][keep])
  if (ncol(gff) >= 5L) positions$End <- gff[[5]][keep]
  bin_genome(positions, karyotype, window = window, count_position = "start")
}
