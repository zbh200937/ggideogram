# Dimensionless layout and projection for ideogram_data.
#
# The body width is the unit of the local coordinate system.  All other
# geometry is derived from public layout arguments and biological lengths.

#' Compute an ideogram layout
#'
#' `ideogram_layout()` converts either a karyotype data frame or an
#' `ideogram_data` object into the single dimensionless layout used by every
#' chromosome, marker, track and component layer.
#'
#' @param karyotype A karyotype data frame or `ideogram_data` object.
#' @param ... Arguments dispatched to the class method.
#'
#' @return An `ideogram_layout_v2` object.
#' @export
ideogram_layout <- function(karyotype, ...) {
  UseMethod("ideogram_layout")
}

#' @rdname ideogram_layout
#' @export
ideogram_layout.data.frame <- function(karyotype, ...) {
  ideogram_layout(as_ideogram_data(karyotype), ...)
}

#' Dimensionless ideogram layout
#'
#' This `ideogram_layout()` method is the package's layout core. It uses
#' chromosome-width units and native ggplot2 orientation.
#' Within each row, chromosomes align at their input start (normally 0 bp).
#'
#' @param karyotype An [as_ideogram_data()] result.
#' @param ncol Number of chromosomes per linear row. `NULL` uses one row and
#'   is required for a circular sector ring.
#' @param chromosome_width Width of a chromosome body in layout units.  The
#'   default `1` defines the package's normalized body-width unit.
#' @param chromosome_gap Clear gap between adjacent body columns, in body-width
#'   units. Circular gaps use `gap_angle`.
#' @param max_chr_length Linear display length of the longest chromosome, in
#'   body-width units. Circular geometry uses `radius`.
#' @param row_gap Gap between linear chromosome rows, in body-width units.
#' @param orientation `"vertical"`, `"horizontal"`, or `"circular"`.
#' @param radius Circular chromosome-center radius in body-width units.
#' @param start_angle Circular starting angle, clockwise from the top, in degrees.
#'   `NULL` keeps the closing gap above and to the right: its left boundary is
#'   the vertical radius at 12 o'clock, for either arrangement direction.
#' @param gap_angle Angular gap after each circular chromosome: one angle or
#'   a vector named by chromosome keys. `NULL` uses 2 degrees between chromosomes
#'   and a 20-degree closing gap for track names and numerical scales.
#' @param opening_angle One non-negative angle in degrees controlling the closing
#'   opening independently of the other chromosome gaps. `NULL` uses the closing
#'   value from `gap_angle`, or 20 degrees when both are `NULL`.
#' @param clockwise Arrange circular sectors clockwise when TRUE.
#' @param scale_length `"global"` preserves length comparisons between
#'   chromosomes. `"per_chr"` gives every chromosome the same display length
#'   and is therefore an explicit non-comparable mode.
#' @param curve_points Number of profile samples used by chromosome-shaped
#'   components, including rounded-body clipping masks.
#' @param tracks `NULL` or a declarative [track_layout()]. Declared tracks are
#'   included in chromosome pitch and plot bounds before any layers are added.
#' @param order_by `"input"`, `"genome"` or `"homolog"`. Grouped layouts use
#'   explicitly mapped genome and homolog metadata. Missing homologs leave slots.
#' @param genome_order,homolog_order Optional complete group orders.
#' @param reverse_chr Chromosome identifiers whose display direction is reversed.
#'   For multiple genomes, supply [chr_key()] values. Source bp remain unchanged.
#' @param ... Reserved for future layout components; currently must be empty.
#'
#' @return An `ideogram_layout_v2` object.  Its `chrom` table contains only
#'   dimensionless body coordinates and biological metadata.
#' @rdname ideogram_layout
#' @export
ideogram_layout.ideogram_data <- function(
    karyotype,
    ncol = NULL,
    chromosome_width = 1,
    chromosome_gap = 3,
    max_chr_length = 40,
    row_gap = 2,
    orientation = c("vertical", "horizontal", "circular"),
    scale_length = c("global", "per_chr"),
    curve_points = 32,
    tracks = NULL,
    order_by = c("input", "genome", "homolog"),
    genome_order = NULL,
    homolog_order = NULL,
    reverse_chr = character(),
    radius = 20,
    start_angle = NULL,
    gap_angle = NULL,
    opening_angle = NULL,
    clockwise = TRUE,
    ...) {
  check_unused_args(...)
  orientation <- match.arg(orientation)
  scale_length <- match.arg(scale_length)
  order_by <- match.arg(order_by)
  check_positive_layout(chromosome_width, "chromosome_width")
  check_nonnegative_layout(chromosome_gap, "chromosome_gap")
  check_positive_layout(max_chr_length, "max_chr_length")
  check_nonnegative_layout(row_gap, "row_gap")
  curve_points <- validate_curve_points(curve_points)
  track_geometry <- resolve_track_geometry(tracks, chromosome_width)

  k <- karyotype$karyotype
  if (!is.null(karyotype$chr_selection)) k <- k[match(karyotype$chr_selection, k$.chr), , drop = FALSE]
  if (!is.null(karyotype$view)) {
    window <- karyotype$view
    k <- k[k$.chr == window$chr, , drop = FALSE]
    k$.source_start <- k$.start
    k$.source_end <- k$.end
    k$.start <- window$start
    k$.end <- window$end
  }
  if (orientation == "circular") {
    if (!is.null(ncol)) stopf("Circular layouts use one sector ring; combine separate plots for multiple rings.")
    return(circular_ideogram_layout(karyotype, k, tracks, track_geometry,
      chromosome_width, radius, start_angle, gap_angle, opening_angle, clockwise,
      scale_length, order_by, genome_order, homolog_order, reverse_chr, curve_points))
  }
  chromosome_count <- nrow(k)
  slots <- chromosome_layout_slots(k, ncol, order_by, genome_order, homolog_order)
  k <- k[slots$index, , drop = FALSE]
  per_row <- slots$ncol
  row <- slots$row
  column <- slots$col
  row_count <- max(row) + 1L

  biological_length <- k$.end - k$.start
  if (scale_length == "global") {
    units_per_bp <- max_chr_length / max(biological_length)
    display_length <- biological_length * units_per_bp
  } else {
    units_per_bp <- max_chr_length / biological_length
    display_length <- rep(max_chr_length, chromosome_count)
  }

  row_height <- vapply(seq_len(row_count) - 1L, function(row_id) {
    if (any(row == row_id)) max(display_length[row == row_id]) else max(display_length)
  }, numeric(1))
  total_long <- sum(row_height) + row_gap * max(0L, row_count - 1L)
  row_top <- total_long - c(0, utils::head(cumsum(row_height + row_gap), -1L))

  cell_width <- track_geometry$left + chromosome_width + track_geometry$right
  pitch <- cell_width + chromosome_gap
  cross_center <- column * pitch + track_geometry$left + chromosome_width / 2
  body_top <- row_top[row + 1L]
  body_bottom <- body_top - display_length
  cross_extent <- per_row * cell_width +
    max(0L, per_row - 1L) * chromosome_gap

  if (orientation == "vertical") {
    axis_start_x <- cross_center
    axis_start_y <- body_top
    axis_end_x <- cross_center
    axis_end_y <- body_bottom
    x_min <- cross_center - chromosome_width / 2
    x_max <- cross_center + chromosome_width / 2
    y_min <- body_bottom
    y_max <- body_top
    bounds <- list(x = c(0, cross_extent), y = c(0, total_long))
  } else {
    axis_start_x <- total_long - body_top
    axis_start_y <- cross_extent - cross_center
    axis_end_x <- total_long - body_bottom
    axis_end_y <- cross_extent - cross_center
    x_min <- axis_start_x
    x_max <- axis_end_x
    y_min <- axis_start_y - chromosome_width / 2
    y_max <- axis_start_y + chromosome_width / 2
    bounds <- list(x = c(0, total_long), y = c(0, cross_extent))
  }

  chrom <- data.frame(
    .chr = k$.chr,
    .start = k$.start,
    .end = k$.end,
    .centromere_start = k$.centromere_start,
    .centromere_end = k$.centromere_end,
    .length = biological_length,
    .display_length = display_length,
    .units_per_bp = units_per_bp,
    .row = row,
    .col = column,
    .axis_start_x = axis_start_x,
    .axis_start_y = axis_start_y,
    .axis_end_x = axis_end_x,
    .axis_end_y = axis_end_y,
    .x_min = pmin(x_min, x_max),
    .x_max = pmax(x_min, x_max),
    .y_min = pmin(y_min, y_max),
    .y_max = pmax(y_min, y_max),
    stringsAsFactors = FALSE
  )
  if (!is.null(karyotype$view)) {
    chrom$.source_start <- k$.source_start
    chrom$.source_end <- k$.source_end
    chrom$.cut_start <- chrom$.start > chrom$.source_start
    chrom$.cut_end <- chrom$.end < chrom$.source_end
  }

  semantic_fields <- intersect(c(".chr_name", ".genome", ".assembly",
                                 ".homolog_group", ".display_name"), names(k))
  chrom[semantic_fields] <- k[semantic_fields]
  chrom$.name_x <- chrom$.axis_start_x
  chrom$.name_y <- chrom$.axis_start_y
  reverse_chr <- validate_chr(reverse_chr, "reverse_chr")
  unknown <- setdiff(reverse_chr, chrom$.chr)
  if (length(unknown)) stopf("Unknown reverse chromosome: %s.", format_chr_rows(unknown))
  reverse <- chrom$.chr %in% reverse_chr
  chrom$.direction <- ifelse(reverse, -1, 1)
  for (axis in c("x", "y")) {
    start <- paste0(".axis_start_", axis)
    end <- paste0(".axis_end_", axis)
    original <- chrom[[start]][reverse]
    chrom[[start]][reverse] <- chrom[[end]][reverse]
    chrom[[end]][reverse] <- original
  }

  centromere_start <- project_positions_raw(
    chrom, k$.chr, k$.centromere_start, allow_na = TRUE)
  centromere_end <- project_positions_raw(
    chrom, k$.chr, k$.centromere_end, allow_na = TRUE)
  chrom$.centromere_start_x <- centromere_start$x
  chrom$.centromere_start_y <- centromere_start$y
  chrom$.centromere_end_x <- centromere_end$x
  chrom$.centromere_end_y <- centromere_end$y

  structure(
    list(
      data = karyotype,
      chrom = chrom,
      index = stats::setNames(seq_len(chromosome_count), chrom$.chr),
      ncol = per_row,
      nrow = row_count,
      row_height = row_height,
      chromosome_width = chromosome_width,
      chromosome_gap = chromosome_gap,
      row_gap = row_gap,
      max_chr_length = max_chr_length,
      curve_points = curve_points,
      orientation = orientation,
      scale_length = scale_length,
      order_by = order_by,
      length_comparable = identical(scale_length, "global"),
      track_layout = tracks %||% track_layout(),
      tracks = track_geometry$table,
      track_ranges = initialize_track_ranges(
        track_geometry$table, chrom$.chr),
      track_extent = c(left = track_geometry$left,
                       right = track_geometry$right),
      attachments = list(chr = list(), locus = list()),
      bounds = bounds,
      units = "chromosome_width"
    ),
    class = c("ideogram_layout_v2", "ideogram_layout")
  )
}

#' Project genomic points into an ideogram layout
#'
#' These low-level functions are the stable bridge for ggplot2 and third-party
#' extensions.  They append projected fields to the supplied data and never
#' create a grob.
#' In circular layouts, projected x is unfolded centerline arc length and y is
#' radius. These are data coordinates consumed by the ideogram Coord, which
#' places native layers in the Cartesian panel. Original bp fields are retained.
#'
#' @param layout An `ideogram_layout_v2` object or a plot from [ggideogram()].
#' @param data A data frame.
#' @param chr,position,start,end Column names in `data`.
#'
#' @return `data` with stable prefixed projection columns.
#' @export
project_chr_point <- function(layout, data, chr, position) {
  layout <- chr_projection_layout(layout)
  check_layout_v2(layout)
  check_projection_data(data)
  chr_value <- projection_column(data, chr, "chr")
  position_value <- projection_column(data, position, "position")
  point <- project_positions_checked(layout, chr_value, position_value)
  data$.chr <- as.character(chr_value)
  data$.position <- as.numeric(position_value)
  data$.x <- point$x
  data$.y <- point$y
  data$.chr_index <- point$index
  data$.row <- layout$chrom$.row[point$index]
  data$.col <- layout$chrom$.col[point$index]
  data
}

#' @rdname project_chr_point
#' @export
project_chr_interval <- function(layout, data, chr, start, end) {
  layout <- chr_projection_layout(layout)
  check_layout_v2(layout)
  check_projection_data(data)
  chr_value <- projection_column(data, chr, "chr")
  start_value <- projection_column(data, start, "start")
  end_value <- projection_column(data, end, "end")
  if (!is.numeric(start_value) || !is.numeric(end_value) ||
      anyNA(start_value) || anyNA(end_value) ||
      any(!is.finite(start_value)) || any(!is.finite(end_value))) {
    stopf("Projected interval `start` and `end` must be finite numeric coordinates.")
  }
  if (any(start_value > end_value)) {
    stopf("Projected intervals require `start <= end`.")
  }
  first <- project_positions_checked(layout, chr_value, start_value)
  last <- project_positions_checked(layout, chr_value, end_value)
  data$.chr <- as.character(chr_value)
  data$.start <- as.numeric(start_value)
  data$.end <- as.numeric(end_value)
  data$.x_start <- first$x
  data$.y_start <- first$y
  data$.x_end <- last$x
  data$.y_end <- last$y
  data$.chr_index <- first$index
  data$.row <- layout$chrom$.row[first$index]
  data$.col <- layout$chrom$.col[first$index]
  data
}

#' Invert projected ideogram points
#'
#' Projects an x/y position orthogonally onto the named chromosome's long axis
#' and returns its genomic coordinate.  Transverse track offsets therefore do
#' not alter the recovered genomic position.
#' Circular input uses the same unfolded arc-length/radius coordinates returned
#' by [project_chr_point()] and [project_chr_track()].
#'
#' @param layout An `ideogram_layout_v2` object or a plot from [ggideogram()].
#' @param data A data frame.
#' @param chr,x,y Column names in `data`.
#' @param tolerance Allowed numerical distance beyond a chromosome endpoint,
#'   in layout units.
#'
#' @return `data` with `.chr` and `.position` appended.
#' @export
unproject_chr_point <- function(layout, data, chr, x, y,
                                tolerance = sqrt(.Machine$double.eps)) {
  layout <- chr_projection_layout(layout)
  check_layout_v2(layout)
  check_projection_data(data)
  check_nonnegative_layout(tolerance, "tolerance")
  chr_value <- as.character(projection_column(data, chr, "chr"))
  x_value <- projection_column(data, x, "x")
  y_value <- projection_column(data, y, "y")
  if (!is.numeric(x_value) || !is.numeric(y_value) || anyNA(x_value) ||
      anyNA(y_value) || any(!is.finite(x_value)) || any(!is.finite(y_value))) {
    stopf("Projected `x` and `y` must be finite numeric coordinates.")
  }
  index <- match_layout_chr(layout, chr_value)
  g <- layout$chrom[index, , drop = FALSE]
  dx <- g$.axis_end_x - g$.axis_start_x
  dy <- g$.axis_end_y - g$.axis_start_y
  length2 <- dx^2 + dy^2
  fraction <- ((x_value - g$.axis_start_x) * dx +
                 (y_value - g$.axis_start_y) * dy) / length2
  axis_tolerance <- tolerance / sqrt(length2)
  outside <- fraction < -axis_tolerance | fraction > 1 + axis_tolerance
  if (any(outside)) {
    stopf("Projected point lies beyond the chromosome axis for %s.",
          format_chr_rows(chr_value[outside]))
  }
  fraction <- pmin(pmax(fraction, 0), 1)
  data$.chr <- chr_value
  data$.position <- g$.start + fraction * (g$.end - g$.start)
  data
}

#' @export
print.ideogram_layout_v2 <- function(x, ...) {
  cat("<ideogram_layout_v2>\n")
  cat("  chromosomes : ", nrow(x$chrom), "\n", sep = "")
  cat("  rows/columns: ", x$nrow, "/", x$ncol, "\n", sep = "")
  cat("  orientation : ", x$orientation, "\n", sep = "")
  cat("  length scale: ", x$scale_length,
      if (x$length_comparable) " (comparable)" else " (not comparable)",
      "\n", sep = "")
  cat("  units       : ", x$units, "\n", sep = "")
  invisible(x)
}

validate_layout_ncol <- function(ncol, chromosome_count) {
  if (is.null(ncol)) return(chromosome_count)
  if (!is.numeric(ncol) || length(ncol) != 1L || is.na(ncol) ||
      !is.finite(ncol) || ncol < 1 || ncol != floor(ncol)) {
    stopf("`ncol` must be NULL or one positive whole number.")
  }
  min(as.integer(ncol), chromosome_count)
}

check_positive_layout <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stopf("`%s` must be one positive finite number.", what)
  }
  invisible(x)
}

check_nonnegative_layout <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0) {
    stopf("`%s` must be one non-negative finite number.", what)
  }
  invisible(x)
}

check_layout_v2 <- function(layout) {
  if (!inherits(layout, "ideogram_layout_v2")) {
    stopf("`layout` must be a dimensionless `ideogram_layout_v2` object.")
  }
  invisible(layout)
}

# Unit normal pointing to semantic chromosome side "right".  Horizontal
# layouts are a clockwise rotation of vertical layouts, so their transverse
# normal has the opposite handedness of the usual left normal of the tangent.
# Keeping this rule in one place makes bodies, tracks, markers, axes and insets
# rotate together instead of allowing horizontal tracks to escape the bounds.
chromosome_right_normal <- function(layout, g) {
  dx <- g$.axis_end_x - g$.axis_start_x
  dy <- g$.axis_end_y - g$.axis_start_y
  length <- sqrt(dx^2 + dy^2)
  handedness <- (if (identical(layout$orientation, "horizontal")) -1 else 1) *
    (g$.direction %||% 1)
  list(
    nx = handedness * (-dy / length),
    ny = handedness * (dx / length),
    length = length
  )
}

check_projection_data <- function(data) {
  if (!is.data.frame(data)) {
    stopf("`data` must be a data frame, not %s.", class(data)[1])
  }
  invisible(data)
}

projection_column <- function(data, column, what) {
  if (!is.character(column) || length(column) != 1L || !nzchar(column)) {
    stopf("`%s` must name one column in `data`.", what)
  }
  if (!column %in% names(data)) {
    stopf("`data` has no `%s` column named `%s`.", what, column)
  }
  data[[column]]
}

match_layout_chr <- function(layout, chr) {
  chr <- as.character(chr)
  if (anyNA(chr) || any(!nzchar(chr))) {
    stopf("Projected chromosome identifiers must be non-missing and non-empty.")
  }
  index <- unname(layout$index[chr])
  if (anyNA(index)) {
    stopf("Unknown chromosome%s: %s.",
          if (length(unique(chr[is.na(index)])) > 1L) "s" else "",
          format_chr_rows(chr[is.na(index)]))
  }
  index
}

project_positions_checked <- function(layout, chr, position) {
  if (!is.numeric(position) || anyNA(position) || any(!is.finite(position))) {
    stopf("Projected `position` must contain finite numeric coordinates.")
  }
  index <- match_layout_chr(layout, chr)
  g <- layout$chrom[index, , drop = FALSE]
  outside <- position < g$.start | position > g$.end
  if (any(outside)) {
    detail <- unique(sprintf("%s:%s", as.character(chr[outside]),
                             format(position[outside], trim = TRUE)))
    stopf(paste0(
      "Genomic position is outside its chromosome: %s.\n",
      "  Require chromosome start <= position <= chromosome end."),
      paste0("`", detail, "`", collapse = ", "))
  }
  point <- project_positions_raw(layout$chrom, chr, position)
  point$index <- index
  point
}

project_positions_raw <- function(chrom, chr, position, allow_na = FALSE) {
  index <- match(as.character(chr), chrom$.chr)
  fraction <- (position - chrom$.start[index]) /
    (chrom$.end[index] - chrom$.start[index])
  x <- chrom$.axis_start_x[index] +
    fraction * (chrom$.axis_end_x[index] - chrom$.axis_start_x[index])
  y <- chrom$.axis_start_y[index] +
    fraction * (chrom$.axis_end_y[index] - chrom$.axis_start_y[index])
  if (allow_na) {
    x[is.na(position)] <- NA_real_
    y[is.na(position)] <- NA_real_
  }
  list(x = x, y = y)
}
