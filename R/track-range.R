# Training and validation of track value ranges.

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
