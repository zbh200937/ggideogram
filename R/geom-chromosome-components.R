# Internal additive components backed by the shared dimensionless layout.

new_chromosome_component <- function(kind, params) {
  structure(
    list(kind = kind, params = params),
    class = "ggideogram_chromosome_component"
  )
}

#' Internal chromosome body component
#' @noRd
geom_chromosome <- function(
    fill = ideogram_colour("body"),
    colour = ideogram_colour("outline"),
    linewidth = ideogram_linewidth("body"),
    alpha = NA,
    linetype = 1,
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
    alpha = alpha,
    linetype = linetype,
    curve_points = curve_points,
    cytoband = cytoband,
    scheme = match.arg(cytoband_scheme),
    palette = cytoband_palette,
    bleach = cytoband_bleach
  ))
}

#' Internal cytoband component
#' @noRd
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

#' Internal cytoband alias
#' @noRd
geom_cytoband <- function(...) geom_chr_cytoband(...)

#' Internal chromosome name component
#' @noRd
geom_chr_name <- function(
    size = ideogram_text_size("chromosome"),
    gap = 0.5,
    colour = ideogram_colour("text"),
    family = "",
    position = c("auto", "start", "end", "middle")) {
  check_positive_layout(size, "size")
  check_nonnegative_layout(gap, "gap")
  new_chromosome_component("name", list(
    size = size, gap = gap, colour = colour, family = family,
    position = match.arg(position)
  ))
}

#' Internal chromosome bp-axis component
#' @noRd
geom_chr_axis <- function(
    chr = TRUE,
    side = c("left", "right", "inner", "outer"),
    breaks = NULL,
    n = 6,
    units = c("auto", "bp", "kb", "Mb", "Gb"),
    labels = TRUE,
    gap = 0.3,
    tick_length = 1.5,
    label_gap = 0.25,
    size = ideogram_text_size("bp"),
    family = "",
    colour = ideogram_colour("axis"),
    linewidth = ideogram_linewidth("axis")) {
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
    side = canonical_chr_side(match.arg(side)),
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
        layout, parameters$fill, parameters$curve_points, parameters$alpha))
      if (parameters$cytoband) {
        bands <- chromosome_cytoband_layer(
          layout, parameters$curve_points, parameters$scheme,
          parameters$palette, parameters$bleach, parameters$alpha)
        if (!is.null(bands)) result <- c(result, list(bands))
      }
      c(result, list(chromosome_outline_layer(
        layout, parameters$colour, parameters$linewidth,
        parameters$curve_points, parameters$alpha, parameters$linetype)))
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
      parameters$family, parameters$colour, parameters$position)),
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
  if (object$kind == "axis") {
    plot <- reserve_chromosome_axis_space(plot, layers)
  }
  if (object$kind == "name") plot <- reserve_chromosome_axis_space(plot, list())
  plot
}
