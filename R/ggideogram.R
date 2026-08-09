#' Build a grammar-of-graphics chromosome ideogram
#'
#' `ggideogram()` creates semantic
#' chromosome data, computes a dimensionless layout, and returns an ordinary
#' ggplot object with a fixed coordinate ratio. It has no output side effect,
#' page canvas, DPI conversion or secondary compatibility renderer.
#'
#' The constructor provides chromosome bodies, cytobands, names and local
#' base-pair axes. Marker, declared tracks and complete inset plots can all be
#' added with exported ggplot2-style components.
#'
#' @param data A karyotype data frame or an [as_ideogram_data()] object.
#' @param mapping,centromere,cytoband,cytoband_mapping Passed to
#'   [as_ideogram_data()] for a data-frame input.
#' @param ncol,chromosome_width,chromosome_gap,max_chr_length,row_gap,orientation,scale_length Passed to the dimensionless
#'   [ideogram_layout()] method.
#' @param tracks `NULL` or a named [track_layout()] reserved before rendering.
#' @param fill,colour,linewidth Chromosome body style. Sizes use ordinary
#'   ggplot2 units.
#' @param curve_points Number of vertices used to approximate each rounded cap.
#' @param show_names Draw chromosome names.
#' @param name_size,name_gap,name_colour Name size in ggplot2 millimetres, gap
#'   in em, and colour.
#' @param axis `FALSE`, `TRUE`, or a chromosome vector selecting local base-pair
#'   axes.
#' @param axis_side Side of the chromosome axis.
#' @param axis_breaks Explicit breaks, a function receiving `c(start, end)`, or
#'   `NULL` for round automatic breaks.
#' @param axis_n Target automatic break intervals.
#' @param axis_units Label unit.
#' @param axis_gap Clear axis-to-body gap in body-width units. Tick marks start
#'   at the axis spine and extend away from the chromosome, so they do not eat
#'   into this clearance.
#' @param axis_tick_length Tick length in millimetres.
#' @param axis_label_gap Axis-label clearance beyond the tick tip, in em.
#' @param axis_size,axis_colour,axis_linewidth Axis text size, colour and line
#'   width in ordinary ggplot2 units.
#' @param cytoband_scheme,cytoband_palette,cytoband_bleach Cytoband palette
#'   controls passed to [cytoband_colours()].
#' @param padding Coordinate padding in body-width units.
#' @param base_family Base font family.
#'
#' @return A standard ggplot object whose coordinate object owns the single
#'   `ideogram_layout_v2` source of truth.
#' @examples
#' data(human_karyotype, package = "ggideogram")
#' ggideogram(human_karyotype, ncol = 12)
#' @export
ggideogram <- function(
    data,
    mapping = NULL,
    centromere = NULL,
    cytoband = NULL,
    cytoband_mapping = NULL,
    ncol = NULL,
    chromosome_width = 1,
    chromosome_gap = 1,
    max_chr_length = 40,
    row_gap = 2,
    orientation = c("vertical", "horizontal"),
    scale_length = c("global", "per_chr"),
    tracks = NULL,
    fill = "#F7F7F7",
    colour = "#4D4D4D",
    linewidth = 0.4,
    curve_points = 32,
    show_names = TRUE,
    name_size = 3.2,
    name_gap = 0.5,
    name_colour = "#202020",
    axis = FALSE,
    axis_side = c("left", "right"),
    axis_breaks = NULL,
    axis_n = 6,
    axis_units = c("auto", "bp", "kb", "Mb", "Gb"),
    axis_gap = 0.3,
    axis_tick_length = 1.5,
    axis_label_gap = 0.25,
    axis_size = 2.4,
    axis_colour = "#666666",
    axis_linewidth = 0.3,
    cytoband_scheme = c("circos", "biovizbase", "only.centromeres"),
    cytoband_palette = NULL,
    cytoband_bleach = 0,
    padding = 0.5,
    base_family = "") {
  orientation <- match.arg(orientation)
  scale_length <- match.arg(scale_length)
  axis_side <- match.arg(axis_side)
  axis_units <- match.arg(axis_units)
  cytoband_scheme <- match.arg(cytoband_scheme)
  curve_points <- validate_curve_points(curve_points)
  check_positive_layout(name_size, "name_size")
  check_nonnegative_layout(name_gap, "name_gap")
  check_nonnegative_layout(axis_gap, "axis_gap")
  check_nonnegative_layout(axis_tick_length, "axis_tick_length")
  check_nonnegative_layout(axis_label_gap, "axis_label_gap")
  check_positive_layout(axis_size, "axis_size")
  check_nonnegative_layout(linewidth, "linewidth")
  check_nonnegative_layout(axis_linewidth, "axis_linewidth")

  if (!inherits(data, "ideogram_data") && is_cytoband(data)) {
    if (!is.null(cytoband)) {
      stopf(paste0(
        "`data` already looks like a cytoband table; do not also supply ",
        "`cytoband`."))
    }
    cytoband <- cytoband_alias(data)
    data <- cytoband_karyotype(cytoband)
  }

  semantic <- if (inherits(data, "ideogram_data")) {
    supplied <- list(mapping = mapping, centromere = centromere,
                     cytoband = cytoband,
                     cytoband_mapping = cytoband_mapping)
    supplied <- names(supplied)[!vapply(supplied, is.null, logical(1))]
    if (length(supplied)) {
      stopf("An `ideogram_data` input cannot also supply %s.",
            paste0("`", supplied, "`", collapse = ", "))
    }
    data
  } else {
    as_ideogram_data(
      data, mapping = mapping, centromere = centromere,
      cytoband = cytoband, cytoband_mapping = cytoband_mapping)
  }

  layout <- ideogram_layout(
    semantic, ncol = ncol,
    chromosome_width = chromosome_width,
    chromosome_gap = chromosome_gap,
    max_chr_length = max_chr_length,
    row_gap = row_gap,
    orientation = orientation,
    scale_length = scale_length,
    curve_points = curve_points,
    tracks = tracks
  )

  layers <- chromosome_body_layers(
    layout, fill = fill, colour = colour, linewidth = linewidth,
    curve_points = curve_points,
    cytoband_scheme = cytoband_scheme,
    cytoband_palette = cytoband_palette,
    cytoband_bleach = cytoband_bleach
  )
  if (isTRUE(show_names)) {
    layers <- c(layers, list(chromosome_name_layer(
      layout, size = name_size, gap = name_gap,
      family = base_family, colour = name_colour)))
  }
  layers <- c(layers, chromosome_axis_layers(
    layout, chr = axis, side = axis_side,
    breaks = axis_breaks, n = axis_n, units = axis_units,
    labels = TRUE, gap = axis_gap, tick_length = axis_tick_length,
    label_gap = axis_label_gap, size = axis_size,
    family = base_family, colour = axis_colour,
    linewidth = axis_linewidth
  ))

  ggplot2::ggplot() + layers +
    coord_ideogram_next(layout, padding = padding, clip = "off") +
    theme_ideogram(base_family = base_family)
}
