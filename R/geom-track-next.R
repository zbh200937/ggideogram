# Standard ggplot2 geoms adapted to declared chromosome tracks.

StatChrTrack <- ggplot2::ggproto(
  "StatChrTrack", ggplot2::Stat,
  required_aes = c("chr", "position"),
  # Value mapping is validated before building. Leave missing observations
  # for the native Geom so an interior NA breaks a path rather than joining it.
  optional_aes = "value",
  compute_panel = function(data, scales, track, chromosome_separation,
                           circular = FALSE,
                           track_geometry = FALSE,
                           track_clip = "auto",
                           na.rm = FALSE) {
    data$ideogram_chr <- as.character(data$chr)
    data$ideogram_track <- track
    chromosome_code <- match(
      data$ideogram_chr, unique(data$ideogram_chr)) - as.integer(!circular)
    data$ideogram_position_shift <- chromosome_code * chromosome_separation
    data$x <- data$position + data$ideogram_position_shift
    data$y <- data$value
    data$group <- as.integer(interaction(
      data$ideogram_chr, data$group,
      drop = TRUE, lex.order = TRUE
    ))
    if (track_geometry) data$ideogram_track_geometry <- TRUE
    data$ideogram_track_clip <- track_clip
    data
  }
)

#' @noRd
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
  trained <- register_dynamic_track_range(object, layout)
  layout <- trained$layout
  coordinates <- trained$coordinates
  coordinate_names <- trained$coordinate_names
  plot <- update_plot_ideogram_layout(plot, layout)

  geom <- object$geom
  if (!is.null(object$tile_clip)) {
    geom <- resolve_track_tile_geom(layout, object$track, object$tile_clip)
  }
  mapping <- object$mapping
  data <- object$data
  params <- object$params
  position_fields <- intersect(c("xend", "xmin", "xmax"), coordinate_names)
  for (aesthetic in intersect(position_fields, names(params))) {
    column <- paste0(".chr_track_", aesthetic)
    data[[column]] <- coordinates[[aesthetic]]
    mapping[[aesthetic]] <- rlang::new_quosure(rlang::sym(column))
    params[[aesthetic]] <- NULL
  }
  arguments <- c(
    list(
      mapping = mapping,
      data = data,
      stat = StatChrTrack,
      position = object$position,
      track = object$track,
      chromosome_separation = track_chromosome_separation(
        layout, object$data, object$mapping, object$params),
      circular = is_circular_layout(layout),
      na.rm = object$na.rm,
      show.legend = object$show.legend,
      inherit.aes = object$inherit.aes
    ),
    params
  )
  if ("stat" %in% names(formals(geom))) {
    layer <- do.call(geom, arguments)
  } else {
    # Some native constructors fix their Stat internally (geom_col in 3.5).
    # Construct their ordinary layer, then install the semantic identity Stat.
    semantic_params <- arguments[c("track", "chromosome_separation", "circular", "na.rm")]
    arguments[c("stat", "track", "chromosome_separation", "circular")] <- NULL
    arguments["mapping"] <- list(NULL)
    layer <- do.call(geom, arguments)
    layer$mapping <- mapping
    layer$stat <- StatChrTrack
    layer$stat_params <- semantic_params
  }
  # Native rectangles become polygon vertices under a nonlinear Coord. Their
  # generated edges may extend beyond a valid observation's source bp.
  layer$stat_params$track_geometry <- inherits(layer$geom, "GeomRect")
  if (inherits(layer$stat, "StatChrTrack")) layer$stat_params$track_clip <-
    object$tile_clip %||% layout$base_spec$tracks[[object$track]]$clip %||% "auto"
  if (length(position_fields)) {
    layer$position <- track_coordinate_position(layer$position, position_fields)
  }
  if ("label" %in% layer$geom$required_aes) {
    # Text geoms may transform anonymous copies of x/y for native nudges.
    layer$position <- track_label_position(layer$position)
  }
  ggplot2::ggplot_add(layer, plot, object_name)
}

shift_track_coordinates <- function(data, fields) {
  for (field in intersect(fields, names(data))) {
    data[[field]] <- data[[field]] + data$ideogram_position_shift
  }
  data
}

track_coordinate_position <- function(position, fields) {
  ggplot2::ggproto("PositionChrTrack", position,
    setup_params = function(self, data) {
      ggplot2::ggproto_parent(position, self)$setup_params(
        shift_track_coordinates(data, fields))
    },
    setup_data = function(self, data, params) {
      ggplot2::ggproto_parent(position, self)$setup_data(
        shift_track_coordinates(data, fields), params)
    })
}

track_label_position <- function(position) {
  ggplot2::ggproto("PositionChrTrackLabel", position,
    compute_layer = function(self, data, params, layout) {
      data <- ggplot2::ggproto_parent(position, self)$compute_layer(
        data, params, layout)
      if (!nrow(data)) return(data)
      original <- data
      data <- transform_ideogram_track(layout$coord$layout, data)
      if (all(c("x_orig", "y_orig") %in% names(original))) {
        original$x <- original$x_orig
        original$y <- original$y_orig
        original <- transform_ideogram_track(layout$coord$layout, original)
        data$x_orig <- original$x
        data$y_orig <- original$y
      }
      data$ideogram_projected <- TRUE
      data
    })
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

    if (is_circular_layout(coord$layout)) data$ideogram_track_geometry <- TRUE
    mask_chr_geometry(data, panel_params, coord, function(rows) {
      ggplot2::GeomTile$draw_panel(rows, panel_params, coord,
        lineend = lineend, linejoin = linejoin)
    })
  }
)

mask_chr_geometry <- function(data, panel_params, coord, draw) {
  layout <- coord$layout
  chromosome <- as.character(data$ideogram_chr)
  children <- lapply(split(seq_len(nrow(data)), chromosome), function(rows) {
    g <- layout$chrom[layout$index[[chromosome[rows[1]]]], , drop = FALSE]
    silhouette <- chromosome_section_polygon(g, layout$chromosome_width, layout$curve_points)
    silhouette <- coord$transform(silhouette, panel_params)
    mask <- grid::polygonGrob(x = silhouette$x, y = silhouette$y,
      default.units = "native", gp = grid::gpar(fill = "white", col = NA))
    grid::grobTree(draw(data[rows, , drop = FALSE]),
      vp = grid::viewport(mask = grid::as.mask(mask)))
  })
  do.call(grid::grobTree, children)
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

track_initial_mapping <- function(mapping) {
  expression <- rlang::get_expr(mapping)
  if (rlang::is_call(expression, "stage")) {
    expression <- match.call(ggplot2::stage, expression)$start
    mapping <- if (rlang::is_quosure(mapping)) {
      rlang::quo_set_expr(mapping, expression)
    } else expression
  }
  delayed <- function(expression) {
    if (rlang::is_missing(expression) || !rlang::is_call(expression)) return(FALSE)
    if (rlang::is_call(expression, c("after_stat", "after_scale"))) return(TRUE)
    any(vapply(rlang::call_args(expression), delayed, logical(1)))
  }
  if (is.null(expression) || delayed(expression)) return(NULL)
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
      # Only initial aesthetics contribute to groups; later stages remain in
      # the original mapping for the native Stat and Geom to evaluate.
      initial <- track_initial_mapping(mapping[[aesthetic]])
      if (is.null(initial)) next
      value <- rlang::eval_tidy(initial, data = data)
      if (is.factor(value) || is.character(value) || is.logical(value)) {
        groups[[length(groups) + 1L]] <- recycle_semantic(
          value, nrow(data), paste0("mapping$", aesthetic))
      }
    }
  }
  do.call(interaction, c(groups, list(drop = TRUE, lex.order = TRUE)))
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
  geom <- if (!is.null(object$native_geom)) object$geom else ggplot2::geom_ribbon
  do.call(geom, arguments)
}

#' @noRd
NULL

#' @noRd
geom_track_point <- function(mapping = NULL, data = NULL, track,
                             position = "identity", ...,
                             na.rm = FALSE, show.legend = NA,
                             inherit.aes = FALSE) {
  geom_chr_track(
    mapping, data, ggplot2::geom_point, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
}

#' @noRd
geom_track_line <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  geom_chr_track(
    mapping, data, ggplot2::geom_line, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
}

#' @noRd
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

#' @noRd
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

#' @noRd
geom_track_text <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  geom_chr_track(
    mapping, data, ggplot2::geom_text, track,
    position = position, ..., na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes)
}

#' @noRd
NULL

#' @noRd
geom_track_area <- function(mapping = NULL, data = NULL, track,
                            position = "identity", ...,
                            na.rm = FALSE, show.legend = NA,
                            inherit.aes = FALSE) {
  new_track_component(
    "area", mapping, data, ggplot2::geom_ribbon, track, position,
    list(...), na.rm, show.legend, inherit.aes, include_zero = TRUE)
}

#' @noRd
geom_track_ribbon <- function(mapping = NULL, data = NULL, track,
                              position = "identity", ...,
                              na.rm = FALSE, show.legend = NA,
                              inherit.aes = FALSE) {
  new_track_component(
    "ribbon", mapping, data, ggplot2::geom_ribbon, track, position,
    list(...), na.rm, show.legend, inherit.aes)
}

#' @noRd
NULL

#' @noRd
geom_track_boxplot <- function(mapping = NULL, data = NULL, track,
                               position = "dodge2", ...,
                               na.rm = FALSE, show.legend = NA,
                               inherit.aes = FALSE) {
  new_track_component(
    "boxplot", mapping, data, ggplot2::geom_boxplot, track, position,
    list(...), na.rm, show.legend, inherit.aes)
}

#' @noRd
geom_track_violin <- function(mapping = NULL, data = NULL, track,
                              position = "dodge", ...,
                              na.rm = FALSE, show.legend = NA,
                              inherit.aes = FALSE) {
  new_track_component(
    "violin", mapping, data, ggplot2::geom_violin, track, position,
    list(...), na.rm, show.legend, inherit.aes)
}
