# Standard ggplot2 geoms adapted to declared chromosome tracks.

StatChrTrack <- ggplot2::ggproto(
  "StatChrTrack", ggplot2::Stat,
  required_aes = c("chr", "position", "value"),
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

register_dynamic_track_range <- function(object, layout) {
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
  stack_range <- if (!is.null(object$range_prototype)) {
    native_range_track_values(object$range_prototype, object$data, object$mapping)
  } else track_position_value_range(object$position, chr, fields$position, value)
  if (!is.null(stack_range)) {
    chr <- c(chr, stack_range$chr)
    value <- c(value, stack_range$value)
  }
  if (isTRUE(object$include_zero)) {
    present_chr <- unique(chr)
    chr <- c(chr, present_chr)
    value <- c(value, rep(0, length(present_chr)))
  }
  coordinate_names <- intersect(c("xend", "xmin", "xmax", "yend", "ymin", "ymax"),
    union(names(object$mapping), names(object$params)))
  coordinates <- lapply(stats::setNames(coordinate_names, coordinate_names), function(aesthetic) {
    input <- if (aesthetic %in% names(object$params)) object$params[[aesthetic]] else
      rlang::eval_tidy(object$mapping[[aesthetic]], data = object$data)
    recycle_semantic(input, nrow(object$data), paste0("track ", aesthetic))
  })
  for (aesthetic in intersect(c("xend", "xmin", "xmax"), coordinate_names)) {
    validate_track_semantics(layout, fields$chr, coordinates[[aesthetic]], fields$value)
  }
  for (aesthetic in intersect(c("yend", "ymin", "ymax"), coordinate_names)) {
    validate_track_semantics(layout, fields$chr, fields$position, coordinates[[aesthetic]])
    chr <- c(chr, as.character(fields$chr))
    value <- c(value, coordinates[[aesthetic]])
  }
  layout <- register_track_values(layout, object$track, chr, value)
  list(layout = layout, coordinates = coordinates, coordinate_names = coordinate_names)
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

track_tile_value_range <- function(object, value) {
  if (!identical(object$geom, ggplot2::geom_tile) &&
      !inherits(object$native_geom, "GeomTile")) return(NULL)
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
  separation <- genomic_span + max(width, genomic_span * sqrt(.Machine$double.eps))
  if (is_circular_layout(layout)) separation <- max(separation, layout$circular$semantic_separation)
  separation
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
  if (is.null(spec$limits[[1]]) && !track_id %in% (layout$scope_fixed_tracks %||% character())) {
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

distribution_track_layer <- function(object, fields) {
  data <- object$data
  data$.track_x <- fields$position
  data$.track_y <- fields$value
  data$.track_group <- track_group_values(
    data, object$mapping, fields$chr, fields$position)
  mapping <- track_mapping_without(
    object$mapping, c("chr", "position", "value", "group", "x", "y"))
  coordinates <- ggplot2::aes(
    x = .data$.track_x, y = .data$.track_y,
    group = .data$.track_group
  )
  mapping <- combine_track_mapping(mapping, coordinates)
  arguments <- c(
    list(
      mapping = mapping, data = data, position = object$position,
      orientation = "x", na.rm = object$na.rm,
      show.legend = object$show.legend,
      inherit.aes = FALSE
    ),
    object$params
  )
  layer <- do.call(object$geom, arguments)
  layer$data$.track_chr <- as.character(fields$chr)
  layer$mapping$chr <- ggplot2::aes(chr = .data$.track_chr)$chr
  native <- layer$stat
  layer$stat <- ggplot2::ggproto(NULL, native,
    required_aes = union(native$required_aes, "chr"))
  layer
}

distribution_track_values <- function(object, fields) {
  layer <- distribution_track_layer(object, fields)
  data <- ggplot2::ggplot_build(ggplot2::ggplot() + layer)$data[[1]]
  columns <- intersect(c("y", "ymin", "ymax",
    if (isTRUE(object$params$notch)) c("notchlower", "notchupper")), names(data))
  values <- unlist(data[columns], use.names = FALSE)
  list(chr = rep(as.character(data$chr), length(columns)), value = values)
}

build_projected_distribution_layer <- function(object, layout, fields) {
  layer <- distribution_track_layer(object, fields)
  layer$position <- chromosome_distribution_position(
    layer$position, object$track, isTRUE(object$params$notch))
  native <- layer$geom
  layer$geom <- ggplot2::ggproto(NULL, native,
    parameters = function(extra = FALSE) native$parameters(extra),
    draw_panel = function(data, panel_params, coord, ..., flipped_aes = FALSE) {
      native$draw_panel(data, panel_params, coord, ...,
        flipped_aes = coord$layout$orientation == "vertical")
    })
  layer
}

chromosome_distribution_position <- function(position, track, notch) {
  ggplot2::ggproto("PositionChrDistribution", position,
    setup_params = function(data) list(),
    setup_data = function(data, params) data,
    compute_layer = function(self, data, params, layout) {
      if (!nrow(data)) return(data)
      groups <- split(data, data$chr)
      computed <- do.call(rbind, lapply(groups, function(rows) {
        settings <- position$setup_params(rows)
        rows <- position$setup_data(rows, settings)
        position$compute_layer(rows, settings, layout)
      }))
      project_distribution_data(computed, layout$coord$layout, track, notch)
    })
}

project_distribution_data <- function(data, layout, track, notch) {
  if (!nrow(data)) return(data)
  raw <- data
  vertical <- layout$orientation == "vertical"
  project <- function(position, value) project_track_coordinate_pair(
    layout, track, as.character(raw$chr), position, value, check_position = FALSE)
  # Stat and Position use bp/value units; only their finished geometry is mapped.
  for (field in intersect(c("x", "xmin", "xmax"), names(raw))) {
    point <- project(raw[[field]], raw$middle %||% raw$y)
    data[[field]] <- if (vertical) point$y else point$x
  }
  values <- c("y", "ymin", "ymax", "lower", "middle", "upper", "ymin_final", "ymax_final",
    if (notch) c("notchlower", "notchupper"))
  for (field in intersect(values, names(raw))) {
    point <- project(raw$x, raw[[field]])
    data[[field]] <- if (vertical) point$x else point$y
  }
  if ("outliers" %in% names(raw)) {
    data$outliers <- lapply(seq_len(nrow(raw)), function(i) {
      values <- raw$outliers[[i]]
      point <- project_track_coordinate_pair(layout, track,
        rep(as.character(raw$chr[i]), length(values)), rep(raw$x[i], length(values)),
        values, check_position = FALSE)
      if (vertical) point$x else point$y
    })
  }
  # Reversed chromosomes still need ordered bounds for native rectangle geoms.
  for (bounds in list(c("xmin", "xmax"), c("ymin", "ymax"))) {
    if (all(bounds %in% names(data))) {
      lower <- pmin(data[[bounds[1]]], data[[bounds[2]]])
      data[[bounds[2]]] <- pmax(data[[bounds[1]]], data[[bounds[2]]])
      data[[bounds[1]]] <- lower
    }
  }
  data <- ggplot2::flip_data(data, vertical)
  data$flipped_aes <- vertical
  data$ideogram_projected <- TRUE
  data
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

#' @noRd
geom_track_axis <- function(
    track,
    chr = TRUE,
    position = c("auto", "end", "start", "gap"),
    breaks = NULL,
    n = 1,
    labels = scales::label_number(),
    tick_length = 1.5,
    label_gap = 0.25,
    size = ideogram_text_size("value"),
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
  if (length(layers) && inherits(layers[[3]]$geom, "GeomCircularGapText")) {
    layout$circular_gap_axes[[object$track]] <- object
    plot <- update_plot_ideogram_layout(plot, layout)
  }
  for (index in seq_along(layers)) {
    plot <- ggplot2::ggplot_add(
      layers[[index]], plot, paste0(object_name, "[[", index, "]]"))
  }
  reserve_chromosome_axis_space(plot, layers)
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
  spec <- track_table_row(layout, object$track)
  position <- object$position
  if (position == "auto") position <- if (is_circular_layout(layout) &&
      spec$value_scale != "per_chr" && isTRUE(object$chr)) "gap" else "end"
  in_gap <- position == "gap"
  if (in_gap && !is_circular_layout(layout)) stopf("Track-axis `position = \"gap\"` requires a circular layout.")
  if (in_gap && spec$value_scale == "per_chr" && length(chromosomes) > 1L) {
    stopf("A gap axis with `value_scale = \"per_chr\"` must select one chromosome.")
  }
  if (in_gap) chromosomes <- if (isTRUE(object$chr))
    layout$chrom$.chr[if (layout$circular$clockwise) nrow(layout$chrom) else 1L] else chromosomes[1L]
  spine <- list()
  ticks <- list()
  text <- list()
  for (index in seq_along(chromosomes)) {
    chr <- chromosomes[index]
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
    genomic_position <- if (position == "end") g$.end else g$.start
    ends <- project_track_values_raw(
      layout, object$track, rep(chr, 2), rep(genomic_position, 2), limits)
    at <- project_track_values_raw(
      layout, object$track, rep(chr, length(breaks)),
      rep(genomic_position, length(breaks)), breaks)
    dx <- g$.axis_end_x - g$.axis_start_x
    dy <- g$.axis_end_y - g$.axis_start_y
    axis_length <- sqrt(dx^2 + dy^2)
    direction <- if (position == "end") 1 else -1
    if (in_gap) {
      arc <- layout$circular$label_gap$left_arc
      ends <- circular_xy(layout, rep(arc, 2), ends$y)
      at <- circular_xy(layout, rep(arc, length(breaks)), at$y)
      theta <- layout$circular$label_gap$left_angle * pi / 180
      dx <- cos(theta); dy <- -sin(theta); axis_length <- direction <- 1
      if (abs(dx) < 1e-12) dx <- 0
      if (abs(dy) < 1e-12) dy <- 0
    }
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
  if (in_gap) {
    spine$ideogram_cartesian <- ticks$ideogram_cartesian <- text$ideogram_cartesian <- TRUE
  }
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
      geom = if (in_gap) GeomCircularGapText else GeomIdeogramAxisText,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = object$tick_length,
        label_gap = object$label_gap,
        size = object$size, colour = object$colour,
        family = object$family, na.rm = FALSE)
    )
  )
  if (in_gap) layers[[3]]$geom_params$gap_track <- object$track
  for (part in seq_along(layers)) {
    layers[[part]]$stat <- StatTrackAxis
    layers[[part]]$stat_params <- list(
      na.rm = FALSE, axis_spec = object,
      chromosomes = chromosomes, axis_part = part)
  }
  layers
}
