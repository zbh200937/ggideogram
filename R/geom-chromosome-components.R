# Public additive components backed by the shared dimensionless layout.

new_chromosome_component <- function(kind, params) {
  structure(
    list(kind = kind, params = params),
    class = "ggideogram_chromosome_component"
  )
}

#' Add chromosome bodies to an ideogram
#'
#' This public component uses the same layout and geometry as
#' [ggideogram()]. It is useful when a base plot was constructed without the
#' desired styling; it never creates a second coordinate system.
#'
#' @param fill,colour,linewidth Standard chromosome body aesthetics.
#' @param curve_points Vertices used to approximate each rounded cap.
#' @param cytoband Include cytobands stored in the plot's semantic data.
#' @param cytoband_scheme,cytoband_palette,cytoband_bleach Cytoband styling.
#'
#' @return An additive ggplot component.
#' @export
geom_chromosome <- function(
    fill = "#F7F7F7",
    colour = "#4D4D4D",
    linewidth = 0.4,
    curve_points = 32,
    cytoband = TRUE,
    cytoband_scheme = c("circos", "biovizbase", "only.centromeres"),
    cytoband_palette = NULL,
    cytoband_bleach = 0) {
  check_nonnegative_layout(linewidth, "linewidth")
  curve_points <- validate_curve_points(curve_points)
  if (!is.logical(cytoband) || length(cytoband) != 1L || is.na(cytoband)) {
    stopf("`cytoband` must be `TRUE` or `FALSE`.")
  }
  new_chromosome_component("body", list(
    fill = fill,
    colour = colour,
    linewidth = linewidth,
    curve_points = curve_points,
    cytoband = cytoband,
    scheme = match.arg(cytoband_scheme),
    palette = cytoband_palette,
    bleach = cytoband_bleach
  ))
}

#' Add cytobands already attached to an ideogram
#'
#' Supply the band table through `ggideogram(cytoband = ...)`; this component
#' only controls its rendering and therefore cannot drift from the host layout.
#'
#' @param scheme,palette,bleach Passed to [cytoband_colours()].
#' @param curve_points Vertices used to approximate rounded boundaries.
#'
#' @return An additive ggplot component.
#' @export
geom_chr_cytoband <- function(
    scheme = c("circos", "biovizbase", "only.centromeres"),
    palette = NULL,
    bleach = 0,
    curve_points = 32) {
  new_chromosome_component("cytoband", list(
    scheme = match.arg(scheme),
    palette = palette,
    bleach = bleach,
    curve_points = validate_curve_points(curve_points)
  ))
}

#' @rdname geom_chr_cytoband
#' @param ... Passed to [geom_chr_cytoband()].
#' @export
geom_cytoband <- function(...) geom_chr_cytoband(...)

#' Add chromosome names
#'
#' @param size,gap,colour,family Standard text settings. `gap` is typographic
#'   spacing in em and remains independent of panel scaling.
#'
#' @return An additive ggplot component.
#' @export
geom_chr_name <- function(
    size = 3.2,
    gap = 0.5,
    colour = "#202020",
    family = "") {
  check_positive_layout(size, "size")
  check_nonnegative_layout(gap, "gap")
  new_chromosome_component("name", list(
    size = size, gap = gap, colour = colour, family = family
  ))
}

#' Add local base-pair axes
#'
#' The axis uses the chromosome's own bp projection. Tick length and text size
#' are physical ggplot2 units; `gap` is a clear distance from the chromosome in
#' body-width units, and ticks extend only outwards.
#'
#' @param chr `TRUE` for every chromosome or a chromosome vector.
#' @param side Axis side in the chromosome's local orientation.
#' @param breaks Explicit breaks, a function of `c(start, end)`, or `NULL`.
#' @param n Target number of automatic intervals.
#' @param units Label unit.
#' @param labels Draw labels as well as ticks.
#' @param gap Clear chromosome-to-axis distance in body-width units.
#' @param tick_length Tick length in millimetres.
#' @param label_gap Typographic label clearance beyond the tick tip, in em.
#' @param size,family,colour,linewidth Standard ggplot2 appearance settings.
#'
#' @return An additive ggplot component.
#' @export
geom_chr_axis <- function(
    chr = TRUE,
    side = c("left", "right"),
    breaks = NULL,
    n = 6,
    units = c("auto", "bp", "kb", "Mb", "Gb"),
    labels = TRUE,
    gap = 0.3,
    tick_length = 1.5,
    label_gap = 0.25,
    size = 2.4,
    family = "",
    colour = "#666666",
    linewidth = 0.3) {
  if (!isTRUE(chr) && (!is.character(chr) || !length(chr))) {
    stopf("`chr` must be `TRUE` or a non-empty chromosome vector.")
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 1) {
    stopf("`n` must be one positive number.")
  }
  if (!is.logical(labels) || length(labels) != 1L || is.na(labels)) {
    stopf("`labels` must be `TRUE` or `FALSE`.")
  }
  check_nonnegative_layout(gap, "gap")
  check_nonnegative_layout(tick_length, "tick_length")
  check_nonnegative_layout(label_gap, "label_gap")
  check_positive_layout(size, "size")
  check_nonnegative_layout(linewidth, "linewidth")
  new_chromosome_component("axis", list(
    chr = chr,
    side = match.arg(side),
    breaks = breaks,
    n = n,
    units = match.arg(units),
    labels = labels,
    gap = gap,
    tick_length = tick_length,
    label_gap = label_gap,
    size = size,
    family = family,
    colour = colour,
    linewidth = linewidth
  ))
}

#' @method ggplot_add ggideogram_chromosome_component
#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.ggideogram_chromosome_component <- function(
    object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  object_name <- ggplot_add_object_name(
    ..., fallback = "chromosome_component")
  parameters <- object$params
  layers <- switch(
    object$kind,
    body = {
      result <- list(chromosome_fill_layer(
        layout, parameters$fill, parameters$curve_points))
      if (parameters$cytoband) {
        bands <- chromosome_cytoband_layer(
          layout, parameters$curve_points, parameters$scheme,
          parameters$palette, parameters$bleach)
        if (!is.null(bands)) result <- c(result, list(bands))
      }
      c(result, list(chromosome_outline_layer(
        layout, parameters$colour, parameters$linewidth,
        parameters$curve_points)))
    },
    cytoband = {
      layer <- chromosome_cytoband_layer(
        layout, parameters$curve_points, parameters$scheme,
        parameters$palette, parameters$bleach)
      if (is.null(layer)) {
        stopf(paste0(
          "This ideogram has no cytoband data. Supply it through ",
          "`ggideogram(cytoband = ...)`."))
      }
      list(layer)
    },
    name = list(chromosome_name_layer(
      layout, parameters$size, parameters$gap,
      parameters$family, parameters$colour)),
    axis = chromosome_axis_layers(
      layout,
      chr = parameters$chr,
      side = parameters$side,
      breaks = parameters$breaks,
      n = parameters$n,
      units = parameters$units,
      labels = parameters$labels,
      gap = parameters$gap,
      tick_length = parameters$tick_length,
      label_gap = parameters$label_gap,
      size = parameters$size,
      family = parameters$family,
      colour = parameters$colour,
      linewidth = parameters$linewidth
    ),
    stopf("Unknown chromosome component `%s`.", object$kind)
  )
  for (index in seq_along(layers)) {
    plot <- ggplot2::ggplot_add(
      layers[[index]], plot, paste0(object_name, "[[", index, "]]"))
  }
  plot
}
