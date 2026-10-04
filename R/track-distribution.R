# Native distribution statistics and their chromosome projection.

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
