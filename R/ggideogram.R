#' Build a grammar-of-graphics chromosome ideogram
#'
#' `ggideogram()` creates semantic
#' chromosome data, computes a dimensionless layout, and returns an ordinary
#' ggplot object with a fixed coordinate ratio.
#'
#' The constructor provides chromosome bodies, cytobands, names and local
#' base-pair axes. Marker, declared tracks and complete inset plots can all be
#' added with exported ggplot2-style components.
#'
#' @param data A karyotype data frame or an [as_ideogram_data()] object.
#' @param mapping,centromere,cytoband,cytoband_mapping Passed to
#'   [as_ideogram_data()] for a data-frame input.
#' @param ncol,chromosome_width,max_chr_length,row_gap,orientation,scale_length,order_by,genome_order,homolog_order,reverse_chr,radius,clockwise Passed to the dimensionless
#'   [ideogram_layout()] method.
#' @param chromosome_gap Clear space between chromosome columns, in body-width
#'   units. `NULL` uses 8 when constructing a vertical plot with a bp axis,
#'   and 3 otherwise. An explicit value is retained.
#' @param start_angle Circular starting angle, clockwise from the top, in degrees.
#'   `NULL` keeps the closing gap above and to the right: its left boundary is
#'   the vertical radius at 12 o'clock, for either arrangement direction.
#' @param gap_angle Angular gap after each circular chromosome: one angle or
#'   a vector named by chromosome keys. `NULL` uses 2 degrees between chromosomes
#'   and a 20-degree closing gap for track names and numerical scales.
#' @param opening_angle One non-negative angle in degrees controlling the closing
#'   opening independently of the other chromosome gaps. `NULL` uses the closing
#'   value from `gap_angle`, or 20 degrees when both are `NULL`.
#' @param tracks `NULL` or a named list of [geom_track()] declarations.
#' @param fill,colour,linewidth,alpha,linetype Chromosome body style. Sizes use ordinary
#'   ggplot2 units.
#' @param curve_points Number of samples for chromosome boundaries, including
#'   rounded caps in linear layouts and arcs in circular layouts.
#' @param show_names Draw chromosome names.
#' @param name_size,name_gap,name_colour Name size in ggplot2 millimetres, gap
#'   in em, and colour.
#' @param name_position Name placement: `"auto"`, `"start"`, `"end"`, or
#'   `"middle"`. Endpoints refer to display direction. Circular names share an
#'   outer label ring.
#' @param axis `FALSE`, `TRUE`, or a chromosome vector selecting local base-pair
#'   axes.
#' @param axis_side Side of the chromosome axis. Circular `"inner"`/`"outer"`
#'   correspond to `"left"`/`"right"`; its default is outer.
#' @param axis_breaks Explicit breaks, a function receiving `c(start, end)`, or
#'   `NULL` for round automatic breaks. Circular automatic breaks omit the
#'   terminal endpoint to keep adjacent seam labels separate.
#' @param axis_n Target automatic break intervals.
#' @param axis_units Label unit.
#' @param axis_gap Clear gap beyond the outermost same-side track or marker lane
#'   (or the body when neither is present), in body-width units. Tick marks start
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
#' @param chr Optional chromosome identifiers to select, in display order.
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
    chromosome_gap = NULL,
    max_chr_length = 40,
    row_gap = 2,
    orientation = c("vertical", "horizontal", "circular"),
    scale_length = c("global", "per_chr"),
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
    fill = ideogram_colour("body"),
    colour = ideogram_colour("outline"),
    linewidth = ideogram_linewidth("body"),
    alpha = NA,
    linetype = 1,
    curve_points = 32,
    show_names = TRUE,
    name_size = ideogram_text_size("chromosome"),
    name_gap = 0.5,
    name_colour = ideogram_colour("text"),
    name_position = c("auto", "start", "end", "middle"),
    axis = FALSE,
    axis_side = c("left", "right", "inner", "outer"),
    axis_breaks = NULL,
    axis_n = 6,
    axis_units = c("auto", "bp", "kb", "Mb", "Gb"),
    axis_gap = 0.3,
    axis_tick_length = 1.5,
    axis_label_gap = 0.25,
    axis_size = ideogram_text_size("bp"),
    axis_colour = ideogram_colour("axis"),
    axis_linewidth = ideogram_linewidth("axis"),
    cytoband_scheme = c("circos", "biovizbase", "only.centromeres"),
    cytoband_palette = NULL,
    cytoband_bleach = 0,
    padding = 0.5,
    base_family = "",
    chr = NULL) {
  orientation <- match.arg(orientation)
  scale_length <- match.arg(scale_length)
  order_by <- match.arg(order_by)
  name_position <- match.arg(name_position)
  tracks <- normalise_chr_tracks(tracks)
  outer_axis <- orientation == "circular" && missing(axis_side)
  axis_side <- canonical_chr_side(match.arg(axis_side))
  if (outer_axis) axis_side <- "right"
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

  if (!is.null(chr)) {
    chr <- validate_chr(chr, "chr")
    if (!length(chr) || anyDuplicated(chr) || any(!chr %in% semantic$karyotype$.chr)) {
      stopf("`chr` must select unique chromosomes in the source karyotype.")
    }
    if (!is.null(semantic$view) && !semantic$view$chr %in% chr) stopf("`chr` does not include the view chromosome.")
    semantic$chr_selection <- chr
  }

  base_spec <- mget(names(formals(ggideogram)), envir = environment())
  base_spec$data <- semantic
  base_spec[c("mapping", "centromere", "cytoband", "cytoband_mapping")] <- rep(list(NULL), 4)

  layout <- ideogram_layout(
    semantic, ncol = ncol,
    chromosome_width = chromosome_width,
    chromosome_gap = chromosome_gap %||%
      if (orientation == "vertical" && !identical(axis, FALSE)) 8 else 3,
    max_chr_length = max_chr_length,
    row_gap = row_gap,
    orientation = orientation,
    scale_length = scale_length,
    curve_points = curve_points,
    tracks = tracks,
    order_by = order_by, genome_order = genome_order,
    homolog_order = homolog_order, reverse_chr = reverse_chr,
    radius = radius, start_angle = start_angle, gap_angle = gap_angle,
    opening_angle = opening_angle, clockwise = clockwise
  )
  layout$base_spec <- base_spec

  layers <- chromosome_body_layers(
    layout, fill = fill, colour = colour, linewidth = linewidth,
    curve_points = curve_points,
    cytoband_scheme = cytoband_scheme,
    cytoband_palette = cytoband_palette,
    cytoband_bleach = cytoband_bleach, alpha = alpha, linetype = linetype
  )
  if (isTRUE(show_names)) {
    layers <- c(layers, list(chromosome_name_layer(
      layout, size = name_size, gap = name_gap,
      family = base_family, colour = name_colour, position = name_position)))
  }
  axis_layers <- chromosome_axis_layers(
    layout, chr = axis, side = axis_side,
    breaks = axis_breaks, n = axis_n, units = axis_units,
    labels = TRUE, gap = axis_gap, tick_length = axis_tick_length,
    label_gap = axis_label_gap, size = axis_size,
    family = base_family, colour = axis_colour,
    linewidth = axis_linewidth
  )

  for (i in seq_along(layers)) layers[[i]]$ideogram_base <- TRUE
  for (i in seq_along(axis_layers)) axis_layers[[i]]$ideogram_base <- TRUE

  plot <- ggplot2::ggplot() + layers + axis_layers +
    coord_ideogram_next(layout, padding = padding, clip = "off") +
    theme_ideogram(base_family = base_family)
  if (is_circular_layout(layout)) plot <- plot + ggplot2::guides(x = "none", y = "none")
  plot <- add_chr_track_contents(plot, tracks)
  plot$layers <- lapply(plot$layers, record_ideogram_aesthetics)
  reserve_chromosome_axis_space(plot, axis_layers)
}
