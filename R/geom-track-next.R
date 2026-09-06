# Standard ggplot2 geoms adapted to declared chromosome tracks.

StatChrTrack <- ggplot2::ggproto(
  "StatChrTrack", ggplot2::Stat,
  required_aes = c("chr", "position", "value"),
  compute_panel = function(data, scales, track, chromosome_separation,
                           na.rm = FALSE) {
    data$ideogram_chr <- as.character(data$chr)
    data$ideogram_track <- track
    chromosome_code <- match(
      data$ideogram_chr, unique(data$ideogram_chr)) - 1L
    data$ideogram_position_shift <- chromosome_code * chromosome_separation
    data$x <- data$position + data$ideogram_position_shift
    data$y <- data$value
    data$group <- as.integer(interaction(
      data$ideogram_chr, data$group,
      drop = TRUE, lex.order = TRUE
    ))
    data
  }
)

#' Adapt an ordinary ggplot2 geom to a chromosome track
#'
#' `geom_chr_track()` keeps the requested geom unchanged and supplies it with
#' standard `x`/`y` coordinates through a chromosome-aware Stat and Coord. The
#' target track must have been declared in `ggideogram(tracks = ...)`.
#'
#' The minimum generic protocol is an identity-stat geom that consumes `x` and
#' `y` and passes its data to `coord$transform()`. This covers ordinary point,
#' line/path, column/rect, tile and text geoms as well as compatible third-party
#' geoms. Components whose standard geom internally reconstructs data are
#' provided by dedicated wrappers such as [geom_track_ribbon()] and
#' [geom_track_boxplot()].
#'
#' @param mapping Aesthetic mappings containing `chr`, `position` and `value`,
#'   plus aesthetics consumed by `geom`.
#' @param data Track data frame.
#' @param geom A ggplot2 geom constructor such as [ggplot2::geom_line()].
#' @param track One identifier declared by [track_layout()].
#' @param stat Currently `"identity"`. The adapter supplies its own semantic
#'   identity Stat; use a dedicated wrapper for statistical geoms.
#' @param position A standard ggplot2 position adjustment.
#' @param ... Parameters and constant aesthetics passed unchanged to `geom`.
#' @param na.rm Passed to the geom layer.
#' @param show.legend Should the standard ggplot2 guide include this layer?
#' @param inherit.aes Whether to inherit plot aesthetics. The ideogram base plot
#'   has no data mapping, so the default is `FALSE`.
#'
#' @return An additive ggplot component. After addition, the resulting layer
#'   retains the ordinary geom class supplied in `geom`.
#' @examples
#' chromosomes <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
#' signal <- data.frame(
#'   Chr = rep(c("A", "B"), each = 5),
#'   Pos = rep(seq(10, 90, length.out = 5), 2),
#'   Value = c(1:5, 5:1)
#' )
#' tracks <- track_layout(
#'   signal = track(side = "right", width = 1, limits = c(0, 5))
#' )
#' ggideogram(chromosomes, tracks = tracks) +
#'   geom_chr_track(
#'     data = signal,
#'     geom = ggplot2::geom_line,
#'     track = "signal",
#'     mapping = ggplot2::aes(chr = Chr, position = Pos, value = Value),
#'     linewidth = 0.4
#'   )
#' @export
geom_chr_track <- function(
    mapping = NULL,
    data = NULL,
    geom,
    track,
    stat = "identity",
    position = "identity",
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = FALSE) {
  if (!is.function(geom)) {
    stopf("`geom` must be a ggplot2 geom constructor function.")
  }
  if (!identical(stat, "identity")) {
    stopf(paste0(
      "The generic chromosome-track protocol currently requires ",
      "`stat = \"identity\"`; use a dedicated `geom_track_*()` wrapper ",
      "for statistical geoms."))
  }
  new_track_component(
    kind = "dynamic", mapping = mapping, data = data,
    geom = geom, track = track, position = position,
    params = list(...), na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, include_zero = FALSE
  )
}

new_track_component <- function(
    kind, mapping, data, geom, track, position, params,
    na.rm, show.legend, inherit.aes, include_zero = FALSE) {
  structure(
    list(
      kind = kind, mapping = mapping, data = data, geom = geom,
      track = track, position = position, params = params,
      na.rm = na.rm, show.legend = show.legend,
      inherit.aes = inherit.aes, include_zero = include_zero
    ),
    class = "ggideogram_track_component"
  )
}

#' @method ggplot_add ggideogram_track_component
#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.ggideogram_track_component <- function(object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  object_name <- ggplot_add_object_name(..., fallback = "chromosome_track")
  if (is.null(object$data) || !is.data.frame(object$data)) {
    stopf("A chromosome track component currently requires data-frame `data`.")
  }
  if (is.null(object$mapping)) {
    stopf("A chromosome track component requires an aesthetic `mapping`.")
  }
  track_table_row(layout, object$track)

  if (object$kind == "dynamic") {
    result <- add_dynamic_track_component(object, plot, layout, object_name)
  } else {
    result <- add_preprojected_track_component(
      object, plot, layout, object_name)
  }
  result
}

ideogram_plot_layout <- function(plot) {
  layout <- plot$coordinates$layout
  if (!inherits(layout, "ideogram_layout_v2")) {
    stopf("Chromosome track components can only be added to `ggideogram()`.")
  }
  layout
}

update_plot_ideogram_layout <- function(plot, layout) {
  coordinate <- plot$coordinates
  plot$coordinates <- coord_ideogram_next(
    layout,
    padding = coordinate$padding %||% 0.5,
    clip = coordinate$clip %||% "off"
  )
  plot
}

add_dynamic_track_component <- function(object, plot, layout, object_name) {
  fields <- mapped_fields(
    object$data, object$mapping, c("chr", "position", "value"),
    what = "mapping")
  validate_track_semantics(layout, fields$chr, fields$position, fields$value)
  chr <- as.character(fields$chr)
  value <- as.numeric(fields$value)
  tile_range <- track_tile_value_range(object, value)
  if (!is.null(tile_range)) {
    chr <- rep(chr, 3L)
    value <- c(value, tile_range$lower, tile_range$upper)
  }
  stack_range <- track_position_value_range(
    object$position, chr, fields$position, value)
  if (!is.null(stack_range)) {
    chr <- c(chr, stack_range$chr)
    value <- c(value, stack_range$value)
  }
  if (isTRUE(object$include_zero)) {
    present_chr <- unique(chr)
    chr <- c(chr, present_chr)
    value <- c(value, rep(0, length(present_chr)))
  }
  layout <- register_track_values(layout, object$track, chr, value)
  plot <- update_plot_ideogram_layout(plot, layout)

  geom <- object$geom
  if (!is.null(object$tile_clip)) {
    geom <- resolve_track_tile_geom(layout, object$track, object$tile_clip)
  }
  arguments <- c(
    list(
      mapping = object$mapping,
      data = object$data,
      stat = StatChrTrack,
      position = object$position,
      track = object$track,
      chromosome_separation = track_chromosome_separation(
        layout, object$data, object$mapping, object$params),
      na.rm = object$na.rm,
      show.legend = object$show.legend,
      inherit.aes = object$inherit.aes
    ),
    object$params
  )
  layer <- do.call(geom, arguments)
  ggplot2::ggplot_add(layer, plot, object_name)
}

resolve_track_tile_geom <- function(layout, track_id, clip) {
  spec <- track_table_row(layout, track_id)
  clipped <- identical(clip, "on") ||
    (identical(clip, "auto") && identical(spec$side, "overlay"))
  if (clipped && !identical(spec$side, "overlay")) {
    stopf(paste0(
      "Chromosome-body tile clipping is available only for an overlay ",
      "track; track `%s` is on the `%s` side."), track_id, spec$side)
  }
  if (clipped) geom_ideogram_tile else ggplot2::geom_tile
}

geom_ideogram_tile <- function(
    mapping = NULL, data = NULL, stat = "identity", position = "identity",
    ..., na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomIdeogramTile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang::list2(na.rm = na.rm, ...)
  )
}

GeomIdeogramTile <- ggplot2::ggproto(
  "GeomIdeogramTile", ggplot2::GeomTile,
  draw_panel = function(
      data, panel_params, coord, lineend = "butt", linejoin = "mitre",
      na.rm = FALSE) {
    if (!nrow(data)) return(grid::nullGrob())
    if (!inherits(coord$layout, "ideogram_layout_v2") ||
        !"ideogram_chr" %in% names(data)) {
      return(ggplot2::GeomTile$draw_panel(
        data, panel_params, coord,
        lineend = lineend, linejoin = linejoin))
    }

    layout <- coord$layout
    chromosome <- as.character(data$ideogram_chr)
    children <- lapply(split(seq_len(nrow(data)), chromosome), function(rows) {
      chr <- chromosome[rows[1]]
      g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
      silhouette <- chromosome_section_polygon(
        g, layout$chromosome_width, layout$curve_points)
      silhouette <- coord$transform(silhouette, panel_params)
      mask <- grid::polygonGrob(
        x = silhouette$x, y = silhouette$y,
        default.units = "native",
        gp = grid::gpar(fill = "white", col = NA)
      )
      tiles <- ggplot2::GeomTile$draw_panel(
        data[rows, , drop = FALSE], panel_params, coord,
        lineend = lineend, linejoin = linejoin)
      grid::grobTree(
        tiles,
        vp = grid::viewport(mask = grid::as.mask(mask))
      )
    })
    do.call(grid::grobTree, children)
  }
)

track_tile_value_range <- function(object, value) {
  if (!identical(object$geom, ggplot2::geom_tile)) return(NULL)
  height <- if ("height" %in% names(object$mapping)) {
    mapped <- rlang::eval_tidy(object$mapping$height, data = object$data)
    recycle_semantic(mapped, nrow(object$data), "mapping$height")
  } else if (!is.null(object$params$height)) {
    rep(object$params$height, nrow(object$data))
  } else {
    rep(ggplot2::resolution(value, zero = FALSE) * 0.9, length(value))
  }
  if (!is.numeric(height) || anyNA(height) || any(!is.finite(height)) ||
      any(height < 0)) {
    stopf("Tile `height` must contain finite non-negative numbers.")
  }
  list(lower = value - height / 2, upper = value + height / 2)
}

track_chromosome_separation <- function(layout, data, mapping, params) {
  genomic_span <- max(layout$chrom$.end) - min(layout$chrom$.start)
  width <- if ("width" %in% names(mapping)) {
    rlang::eval_tidy(mapping$width, data = data)
  } else {
    params$width %||% 0
  }
  width <- suppressWarnings(max(abs(width), na.rm = TRUE))
  if (!is.finite(width)) width <- 0
  genomic_span + max(width, genomic_span * sqrt(.Machine$double.eps))
}

track_position_value_range <- function(position, chr, genomic_position, value) {
  is_fill <- identical(position, "fill") || inherits(position, "PositionFill")
  is_stack <- identical(position, "stack") ||
    inherits(position, "PositionStack")
  if (!is_fill && !is_stack) return(NULL)
  if (is_fill) {
    return(data.frame(chr = unique(as.character(chr)), value = 1))
  }
  key <- interaction(chr, genomic_position, drop = TRUE, lex.order = TRUE)
  positive <- tapply(pmax(value, 0, na.rm = FALSE), key, sum, na.rm = TRUE)
  negative <- tapply(pmin(value, 0, na.rm = FALSE), key, sum, na.rm = TRUE)
  first <- match(levels(key), key)
  data.frame(
    chr = rep(as.character(chr[first]), each = 2L),
    value = as.vector(rbind(negative, positive)),
    stringsAsFactors = FALSE
  )
}

validate_track_semantics <- function(layout, chr, position, value) {
  chr <- as.character(chr)
  match_layout_chr(layout, chr)
  if (!is.numeric(position) || anyNA(position) ||
      any(!is.finite(position))) {
    stopf("Track `position` must contain finite numeric coordinates.")
  }
  project_positions_checked(layout, chr, position)
  if (!is.numeric(value)) {
    stopf("Track `value` must be numeric.")
  }
  finite <- !is.na(value)
  if (any(!is.finite(value[finite]))) {
    stopf("Track `value` must contain only finite numbers or missing values.")
  }
  invisible(TRUE)
}

require_explicit_track_limits <- function(layout, track_id, component) {
  spec <- track_table_row(layout, track_id)
  if (is.null(spec$limits[[1]])) {
    stopf(paste0(
      "`%s` requires explicit `limits` in `track(%s = ...)` because its ",
      "standard ggplot2 Stat/Geom reconstructs coordinates before drawing."),
      component, "limits")
  }
  spec
}

add_preprojected_track_component <- function(
    object, plot, layout, object_name) {
  require_explicit_track_limits(layout, object$track, object$kind)
  if (object$kind %in% c("area", "boxplot", "violin")) {
    fields <- mapped_fields(
      object$data, object$mapping, c("chr", "position", "value"),
      what = "mapping")
    range_value <- as.numeric(fields$value)
    if (object$kind == "area") range_value <- c(range_value, 0)
  } else if (object$kind == "ribbon") {
    fields <- mapped_fields(
      object$data, object$mapping, c("chr", "position", "ymin", "ymax"),
      what = "mapping")
    if (!is.numeric(fields$ymin) || !is.numeric(fields$ymax) ||
        any(fields$ymin > fields$ymax, na.rm = TRUE)) {
      stopf("Ribbon `ymin` and `ymax` must be numeric with `ymin <= ymax`.")
    }
    range_value <- c(fields$ymin, fields$ymax)
  } else {
    stopf("Unknown preprojected track component `%s`.", object$kind)
  }

  primary_value <- if (object$kind == "ribbon") fields$ymin else fields$value
  validate_track_semantics(
    layout, fields$chr, fields$position, primary_value)
  finite_range <- range_value[!is.na(range_value)]
  layout <- register_track_values(
    layout, object$track,
    rep(as.character(fields$chr),
        length.out = length(range_value)),
    range_value
  )
  plot <- update_plot_ideogram_layout(plot, layout)

  layer <- if (object$kind %in% c("area", "ribbon")) {
    build_projected_ribbon_layer(object, layout, fields)
  } else {
    build_projected_distribution_layer(object, layout, fields)
  }
  ggplot2::ggplot_add(layer, plot, object_name)
}

track_mapping_without <- function(mapping, fields) {
  mapping[setdiff(names(mapping), fields)]
}

combine_track_mapping <- function(mapping, coordinates) {
  mapping[names(coordinates)] <- coordinates
  class(mapping) <- class(coordinates)
  mapping
}

track_group_values <- function(data, mapping, chr, position = NULL) {
  group <- if ("group" %in% names(mapping)) {
    value <- rlang::eval_tidy(mapping$group, data = data)
    recycle_semantic(value, nrow(data), "mapping$group")
  } else if (!is.null(position)) {
    position
  } else {
    rep(1L, nrow(data))
  }
  groups <- list(as.character(chr), group)
  if (!"group" %in% names(mapping)) {
    aesthetics <- setdiff(names(mapping),
      c("chr", "position", "value", "x", "y", "ymin", "ymax"))
    for (aesthetic in aesthetics) {
      value <- rlang::eval_tidy(mapping[[aesthetic]], data = data)
      if (is.factor(value) || is.character(value) || is.logical(value)) {
        groups[[length(groups) + 1L]] <- recycle_semantic(
          value, nrow(data), paste0("mapping$", aesthetic))
      }
    }
  }
  do.call(interaction, c(groups, list(drop = TRUE, lex.order = TRUE)))
}

build_projected_distribution_layer <- function(object, layout, fields) {
  projected <- project_track_values_raw(
    layout, object$track, fields$chr, fields$position, fields$value)
  data <- object$data
  data$.track_x <- projected$x
  data$.track_y <- projected$y
  data$.track_group <- track_group_values(
    data, object$mapping, fields$chr, fields$position)
  mapping <- track_mapping_without(
    object$mapping, c("chr", "position", "value", "group", "x", "y"))
  coordinates <- ggplot2::aes(
    x = .data$.track_x, y = .data$.track_y,
    group = .data$.track_group
  )
  mapping <- combine_track_mapping(mapping, coordinates)
  orientation <- if (layout$orientation == "vertical") "y" else "x"
  arguments <- c(
    list(
      mapping = mapping, data = data, position = object$position,
      orientation = orientation, na.rm = object$na.rm,
      show.legend = object$show.legend,
      inherit.aes = FALSE
    ),
    object$params
  )
  do.call(object$geom, arguments)
}

build_projected_ribbon_layer <- function(object, layout, fields) {
  lower <- if (object$kind == "area") rep(0, length(fields$value)) else fields$ymin
  upper <- if (object$kind == "area") fields$value else fields$ymax
  first <- project_track_values_raw(
    layout, object$track, fields$chr, fields$position, lower)
  last <- project_track_values_raw(
    layout, object$track, fields$chr, fields$position, upper)
  data <- object$data
  if (layout$orientation == "vertical") {
    data$.track_position <- first$y
    data$.track_min <- pmin(first$x, last$x)
    data$.track_max <- pmax(first$x, last$x)
    coordinates <- ggplot2::aes(
      y = .data$.track_position,
      xmin = .data$.track_min,
      xmax = .data$.track_max,
      group = .data$.track_group
    )
    orientation <- "y"
  } else {
    data$.track_position <- first$x
    data$.track_min <- pmin(first$y, last$y)
    data$.track_max <- pmax(first$y, last$y)
    coordinates <- ggplot2::aes(
      x = .data$.track_position,
      ymin = .data$.track_min,
      ymax = .data$.track_max,
      group = .data$.track_group
    )
    orientation <- "x"
  }
  data$.track_group <- track_group_values(
    data, object$mapping, fields$chr, position = NULL)
  mapping <- track_mapping_without(
    object$mapping,
    c("chr", "position", "value", "ymin", "ymax", "group",
      "x", "y", "xmin", "xmax"))
  mapping <- combine_track_mapping(mapping, coordinates)
  arguments <- c(
    list(
      mapping = mapping, data = data, stat = "identity",
      position = object$position, orientation = orientation,
      na.rm = object$na.rm, show.legend = object$show.legend,
      inherit.aes = FALSE
    ),
    object$params
  )
  do.call(ggplot2::geom_ribbon, arguments)
}

#' Common chromosome track geom wrappers
#'
#' These functions are small presets over [geom_chr_track()] and retain
#' ordinary ggplot2 aesthetic, scale and guide semantics. An overlay
#' `geom_track_tile()` uses a [ggplot2::GeomTile] subclass only to mask its
#' standard tile output to the chromosome silhouette. All wrappers require
#' mappings for `chr`, `position` and `value`.
#'
#' @inheritParams geom_chr_track
#' @return An additive chromosome-track component.
#' @name geom_track_geoms
NULL

#' @rdname geom_track_geoms
#' @export
geom_track_point <- function(mapping = NULL, data = NULL, track,
                             position = "identity", ...,
                             na.rm = FALSE, show.legend = NA,
                             inherit.aes = FALSE) {
  geom_chr_track(
    mapping, data, ggplot2::geom_point, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
}

#' @rdname geom_track_geoms
#' @export
geom_track_line <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  geom_chr_track(
    mapping, data, ggplot2::geom_line, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
}

#' @rdname geom_track_geoms
#' @export
geom_track_col <- function(mapping = NULL, data = NULL, track,
                           position = "stack", ...,
                           na.rm = FALSE, show.legend = NA,
                           inherit.aes = FALSE) {
  component <- geom_chr_track(
    mapping, data, ggplot2::geom_col, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
  component$include_zero <- TRUE
  component
}

#' @rdname geom_track_geoms
#' @param clip Tile clipping mode. `"auto"` clips an overlay tile layer to the
#'   chromosome silhouette and leaves beside tracks unchanged. `"on"` requires
#'   an overlay track; `"off"` preserves unrestricted rectangular tiles.
#' @export
geom_track_tile <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            clip = c("auto", "on", "off"),
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  component <- geom_chr_track(
    mapping, data, ggplot2::geom_tile, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
  component$tile_clip <- match.arg(clip)
  component
}

#' @rdname geom_track_geoms
#' @export
geom_track_text <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  geom_chr_track(
    mapping, data, ggplot2::geom_text, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
}

#' Draw filled area and ribbon chromosome tracks
#'
#' These wrappers project explicit track limits first and then delegate drawing
#' to the standard [ggplot2::GeomRibbon]. `geom_track_area()` uses zero as its
#' baseline; `geom_track_ribbon()` maps standard `ymin` and `ymax` aesthetics.
#'
#' @inheritParams geom_chr_track
#' @return An additive chromosome-track component.
#' @name geom_track_ranges
NULL

#' @rdname geom_track_ranges
#' @export
geom_track_area <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  new_track_component(
    "area", mapping, data, ggplot2::geom_ribbon, track, position,
    list(...), na.rm, show.legend, inherit.aes, include_zero = TRUE)
}

#' @rdname geom_track_ranges
#' @export
geom_track_ribbon <- function(mapping = NULL, data = NULL, track,
                              position = "identity", ...,
                              na.rm = FALSE, show.legend = NA,
                              inherit.aes = FALSE) {
  new_track_component(
    "ribbon", mapping, data, ggplot2::geom_ribbon, track, position,
    list(...), na.rm, show.legend, inherit.aes)
}

#' Draw distribution geoms in chromosome tracks
#'
#' These wrappers project observations through an explicitly limited track and
#' then use the unmodified standard [ggplot2::GeomBoxplot] or
#' [ggplot2::GeomViolin] with their standard ggplot2 Stats. Map `chr`,
#' `position` and `value`; `position` identifies the locus/bin at which each
#' distribution is centred. Map `group` when multiple distributions share one
#' locus.
#'
#' @inheritParams geom_chr_track
#' @return An additive chromosome-track component.
#' @name geom_track_distributions
NULL

#' @rdname geom_track_distributions
#' @export
geom_track_boxplot <- function(mapping = NULL, data = NULL, track,
                               position = "dodge2", ...,
                               na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE) {
  new_track_component(
    "boxplot", mapping, data, ggplot2::geom_boxplot, track, position,
    list(...), na.rm, show.legend, inherit.aes)
}

#' @rdname geom_track_distributions
#' @export
geom_track_violin <- function(mapping = NULL, data = NULL, track,
                              position = "dodge", ...,
                              na.rm = FALSE, show.legend = NA,
                              inherit.aes = FALSE) {
  new_track_component(
    "violin", mapping, data, ggplot2::geom_violin, track, position,
    list(...), na.rm, show.legend, inherit.aes)
}

#' Add a local value axis to a chromosome track
#'
#' The axis is drawn at a chromosome endpoint and uses the track's resolved raw
#' value limits and transformation. Its spine is a standard segment. Tick
#' length and text size are physical ggplot2 sizes, and label clearance starts
#' at the tick tip rather than the spine.
#'
#' Add an automatic-range axis after at least one data component has registered
#' that range. An explicitly limited track can receive its axis in any order.
#'
#' @param track One declared track identifier.
#' @param chr `TRUE` for every chromosome or a character vector selecting
#'   chromosomes.
#' @param position Chromosome endpoint at which the transverse axis is drawn.
#' @param breaks `NULL` for pretty breaks, a numeric vector, or a function of
#'   the raw limits.
#' @param n Target number of automatic intervals.
#' @param labels A labelling function or a character vector matching breaks.
#' @param tick_length Tick length in millimetres.
#' @param label_gap Typographic gap beyond the tick tip, in em.
#' @param size,colour,linewidth Standard ggplot2 text, colour and line sizes.
#' @param family Font family.
#'
#' @return An additive track-axis component.
#' @export
geom_track_axis <- function(
    track,
    chr = TRUE,
    position = c("end", "start"),
    breaks = NULL,
    n = 4,
    labels = scales::label_number(),
    tick_length = 1.5,
    label_gap = 0.25,
    size = 2.4,
    colour = "#666666",
    linewidth = 0.3,
    family = "") {
  position <- match.arg(position)
  if (!isTRUE(chr) && (!is.character(chr) || !length(chr))) {
    stopf("`chr` must be `TRUE` or a non-empty chromosome vector.")
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 1) {
    stopf("`n` must be one positive number.")
  }
  check_nonnegative_layout(tick_length, "tick_length")
  check_nonnegative_layout(label_gap, "label_gap")
  check_positive_layout(size, "size")
  check_nonnegative_layout(linewidth, "linewidth")
  structure(
    list(
      track = track, chr = chr, position = position, breaks = breaks,
      n = n, labels = labels, tick_length = tick_length,
      label_gap = label_gap, size = size, colour = colour,
      linewidth = linewidth, family = family
    ),
    class = "ggideogram_track_axis_component"
  )
}

#' @method ggplot_add ggideogram_track_axis_component
#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.ggideogram_track_axis_component <- function(
    object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  object_name <- ggplot_add_object_name(..., fallback = "track_axis")
  track_table_row(layout, object$track)
  chromosomes <- if (isTRUE(object$chr)) {
    layout$chrom$.chr
  } else {
    as.character(object$chr)
  }
  unknown <- setdiff(chromosomes, layout$chrom$.chr)
  if (length(unknown)) {
    stopf("Unknown track-axis chromosome%s: %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  layers <- build_track_axis_layers(layout, object, chromosomes)
  for (index in seq_along(layers)) {
    plot <- ggplot2::ggplot_add(
      layers[[index]], plot, paste0(object_name, "[[", index, "]]"))
  }
  plot
}

# Resolve axes against the final coordinate layout, after all track layers have
# registered their values. Keep the ordinary segment/tick/text Geoms intact.
StatTrackAxis <- ggplot2::ggproto(
  "StatTrackAxis", ggplot2::StatIdentity,
  extra_params = c("na.rm", "axis_spec", "chromosomes", "axis_part"),
  compute_layer = function(self, data, params, layout) {
    layers <- build_track_axis_layers(
      layout$coord$layout, params$axis_spec, params$chromosomes)
    if (!length(layers)) return(data[FALSE, , drop = FALSE])
    result <- layers[[params$axis_part]]$data
    panels <- unique(data$PANEL)
    do.call(rbind, lapply(panels, function(panel) {
      result$PANEL <- panel
      result$group <- -1L
      result
    }))
  }
)

build_track_axis_layers <- function(layout, object, chromosomes) {
  spine <- list()
  ticks <- list()
  text <- list()
  for (index in seq_along(chromosomes)) {
    chr <- chromosomes[index]
    spec <- track_table_row(layout, object$track)
    key <- track_key_for_rows(spec, chr)
    limits <- layout$track_ranges[[key]]$limits
    if (is.null(limits)) {
      stopf(paste0(
        "Track `%s` has no resolved range for chromosome `%s`.\n",
        "  Add its data layer before an automatic-range axis, or declare ",
        "explicit track limits."), object$track, chr)
    }
    breaks <- if (is.function(object$breaks)) {
      object$breaks(limits)
    } else if (is.null(object$breaks)) {
      pretty(limits, n = object$n)
    } else {
      object$breaks
    }
    if (!is.numeric(breaks) || anyNA(breaks) || any(!is.finite(breaks))) {
      stopf("Track-axis `breaks` must resolve to finite numbers.")
    }
    breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]
    if (!length(breaks)) next
    label <- if (is.function(object$labels)) {
      object$labels(breaks)
    } else {
      object$labels
    }
    if (length(label) != length(breaks)) {
      stopf("Track-axis `labels` must match the number of breaks.")
    }

    g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
    genomic_position <- if (object$position == "end") g$.end else g$.start
    ends <- project_track_values_raw(
      layout, object$track, rep(chr, 2), rep(genomic_position, 2), limits)
    at <- project_track_values_raw(
      layout, object$track, rep(chr, length(breaks)),
      rep(genomic_position, length(breaks)), breaks)
    dx <- g$.axis_end_x - g$.axis_start_x
    dy <- g$.axis_end_y - g$.axis_start_y
    axis_length <- sqrt(dx^2 + dy^2)
    direction <- if (object$position == "end") 1 else -1
    spine[[index]] <- data.frame(
      x = ends$x[1], y = ends$y[1],
      xend = ends$x[2], yend = ends$y[2]
    )
    ticks[[index]] <- data.frame(
      x = at$x, y = at$y,
      nx = rep(direction * dx / axis_length, length(breaks)),
      ny = rep(direction * dy / axis_length, length(breaks))
    )
    text[[index]] <- data.frame(
      x = at$x, y = at$y, label = as.character(label),
      nx = rep(direction * dx / axis_length, length(breaks)),
      ny = rep(direction * dy / axis_length, length(breaks))
    )
  }
  spine <- do.call(rbind, spine)
  ticks <- do.call(rbind, ticks)
  text <- do.call(rbind, text)
  if (is.null(spine) || !nrow(spine)) return(list())
  layers <- list(
    ggplot2::geom_segment(
      data = spine,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      colour = object$colour, linewidth = object$linewidth,
      inherit.aes = FALSE, show.legend = FALSE
    ),
    ggplot2::layer(
      data = ticks,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity", geom = GeomIdeogramTick,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = object$tick_length, colour = object$colour,
        linewidth = object$linewidth, na.rm = FALSE)
    ),
    ggplot2::layer(
      data = text,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label,
        nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity",
      geom = GeomIdeogramAxisText,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = object$tick_length,
        label_gap = object$label_gap,
        size = object$size, colour = object$colour,
        family = object$family, na.rm = FALSE)
    )
  )
  for (part in seq_along(layers)) {
    layers[[part]]$stat <- StatTrackAxis
    layers[[part]]$stat_params <- list(
      na.rm = FALSE, axis_spec = object,
      chromosomes = chromosomes, axis_part = part)
  }
  layers
}
