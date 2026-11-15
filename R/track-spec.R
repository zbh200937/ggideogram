# Declarative track specifications for the dimensionless renderer.

#' @noRd
track <- function(
    side = c("right", "left", "overlay", "outer", "inner"),
    width = 1,
    gap = NULL,
    value_scale = c("per_track", "global", "per_chr"),
    limits = NULL,
    transform = "identity",
    reverse = FALSE,
    offset = 0,
    .partial = FALSE) {
  side <- canonical_chr_side(match.arg(side))
  value_scale <- match.arg(value_scale)
  gap <- gap %||% if (side == "overlay") 0 else 0.2
  check_positive_layout(width, "width")
  check_nonnegative_layout(gap, "gap")
  if (!is.numeric(offset) || length(offset) != 1L || !is.finite(offset)) {
    stopf("`offset` must be one finite number.")
  }
  if (!.partial && side == "overlay" && gap != 0) {
    stopf("An overlay track must use `gap = 0` or the default `gap = NULL`.")
  }
  if (!.partial && side != "overlay" && offset != 0) {
    stopf("Only an overlay track can use a non-zero `offset`; set `offset = 0` for a beside track.")
  }
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    stopf("`reverse` must be `TRUE` or `FALSE`.")
  }
  if (!is.null(limits)) {
    if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) ||
        any(!is.finite(limits)) || limits[1] >= limits[2]) {
      stopf("`limits` must be NULL or two finite increasing numbers.")
    }
    limits <- as.numeric(limits)
  }
  transformer <- tryCatch(
    scales::as.transform(transform),
    error = function(error) {
      stopf("Invalid track `transform`: %s", conditionMessage(error))
    }
  )
  if (!is.null(limits)) {
    transformed <- transformer$transform(limits)
    if (anyNA(transformed) || any(!is.finite(transformed)) ||
        transformed[1] >= transformed[2]) {
      stopf("Track `limits` are invalid after the `%s` transformation.",
            transformer$name)
    }
  }

  structure(
    list(
      side = side,
      width = width,
      gap = gap,
      value_scale = value_scale,
      limits = limits,
      transform = transformer,
      transform_name = transformer$name,
      reverse = reverse,
      offset = offset
    ),
    class = "ideogram_track_spec"
  )
}

default_track_width <- function(side, layers) {
  if (side == "overlay") return(1)
  gene <- vapply(layers, function(layer) {
    if (inherits(layer, "ggideogram_object")) layer <- layer$component
    inherits(layer, "ggideogram_gene_component")
  }, logical(1))
  if (any(gene)) 6 else 3
}

validate_track_spec <- function(spec) {
  args <- setdiff(names(formals(track)), ".partial")
  geometry <- spec[args]
  if (isTRUE(spec$auto_width)) {
    layers <- if (is.null(spec$parts)) spec$layers else
      unlist(lapply(spec$parts, `[[`, "layers"), recursive = FALSE)
    geometry$width <- default_track_width(spec$side, layers)
  }
  if (isTRUE(spec$auto_gap)) geometry["gap"] <- list(NULL)
  geometry <- do.call(track, geometry)
  spec[names(geometry)] <- geometry
  spec
}

#' Combine named chromosome tracks
#'
#' Track names are stable identifiers used by every track layer. Tracks on the
#' same side are laid out in declaration order from the chromosome body
#' outwards; adding ggplot2 layers later never changes that order.
#'
#' @param ... Named [geom_track()] declarations. The list names supply each
#'   track identifier and its layers, data and mapping remain together.
#'
#' @return An `ideogram_track_layout` object.
#' @examples
#' tracks <- track_layout(
#'   markers = geom_track(side = "right", width = 0.6),
#'   density = geom_track(side = "right", width = 1.2, limits = c(0, 10)),
#'   counts = geom_track(side = "left", width = 1, limits = c(0, 100))
#' )
#' @export
track_layout <- function(...) {
  specs <- list(...)
  ids <- names(specs)
  if (length(specs) &&
      (is.null(ids) || anyNA(ids) || any(!nzchar(ids)))) {
    stopf("Every track in `track_layout()` must have a non-empty name.")
  }
  if (anyDuplicated(ids)) {
    stopf("Track identifiers must be unique; duplicated: %s.",
          paste0("`", unique(ids[duplicated(ids)]), "`", collapse = ", "))
  }
  bad <- !vapply(specs, inherits, logical(1), "ideogram_track_spec")
  if (any(bad)) {
    stopf("Track%s %s %s not created by `geom_track()`.",
          if (sum(bad) > 1L) "s" else "",
          paste0("`", ids[bad], "`", collapse = ", "),
          if (sum(bad) > 1L) "were" else "was")
  }
  structure(specs, class = c("ideogram_track_layout", "list"))
}

#' @export
print.ideogram_track_spec <- function(x, ...) {
  cat("<ideogram_track_spec>\n")
  cat("  side/value scale: ", x$side, "/", x$value_scale, "\n", sep = "")
  cat("  width/gap       : ", x$width, "/", x$gap, " body widths\n",
      sep = "")
  cat("  transform       : ", x$transform_name,
      if (x$reverse) " (reversed)" else "", "\n", sep = "")
  cat("  limits          : ",
      if (is.null(x$limits)) "automatic" else paste(x$limits, collapse = " .. "),
      "\n", sep = "")
  invisible(x)
}

#' @export
print.ideogram_track_layout <- function(x, ...) {
  cat("<ideogram_track_layout>\n")
  cat("  tracks: ", if (length(x)) paste(names(x), collapse = ", ") else "none",
      "\n", sep = "")
  invisible(x)
}

# Track and column boundaries in the layout's Cartesian coordinates.
track_column_region <- function(layout, chr, spec, column = 1L, columns = 1L) {
  g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
  offsets <- spec$low_offset + c(column - 1, column) / columns * (spec$high_offset - spec$low_offset)
  if (is_circular_layout(layout)) {
    positions <- seq(g$.start, g$.end, length.out = layout$curve_points)
    points <- project_positions_checked(layout, rep(chr, length(positions) * 2),
      rep(positions, 2))
    p <- offset_chr_points_signed(layout, rep(chr, length(points$x)), points,
      rep(offsets, each = length(positions)))
    xy <- as.data.frame(circular_xy(layout, p$x, p$y))
    polygon <- rbind(xy[seq_along(positions), ], xy[rev(seq_along(positions) + length(positions)), ])
    polygon <- rbind(polygon, polygon[1, ])
    return(list(xlim = range(xy$x), ylim = range(xy$y), offset = mean(offsets), polygon = polygon))
  }
  points <- project_positions_checked(layout, rep(chr, 4), rep(c(g$.start, g$.end), each = 2))
  corners <- offset_chr_points_signed(layout, rep(chr, 4), points, rep(offsets, 2))
  polygon <- data.frame(x = range(corners$x)[c(1, 2, 2, 1, 1)],
    y = range(corners$y)[c(1, 1, 2, 2, 1)])
  list(xlim = range(corners$x), ylim = range(corners$y), offset = offsets[1], polygon = polygon)
}

resolve_track_geometry <- function(tracks, chromosome_width, source_chr, visible_chr,
    compact = FALSE, min_gap = 0) {
  tracks <- normalise_chr_tracks(tracks)
  if (is.null(tracks)) tracks <- track_layout()
  marker_extent <- attr(tracks, "marker_extent") %||% c(left = 0, right = 0)
  if (!inherits(tracks, "ideogram_track_layout")) {
    stopf("`tracks` must be NULL or created by `track_layout()`.")
  }
  if (!length(tracks)) {
    table <- data.frame(
      id = character(), side = character(), width = numeric(),
      gap = numeric(), value_scale = character(), reverse = logical(),
      transform_name = character(), offset = numeric(), low_offset = numeric(),
      high_offset = numeric(), stringsAsFactors = FALSE
    )
    table$limits <- I(list())
    table$transform <- I(list())
    return(list(table = table, owners = list(),
      left = marker_extent[["left"]], right = marker_extent[["right"]]))
  }

  owners <- lapply(names(tracks), function(id) {
    spec <- tracks[[id]]
    unknown <- setdiff(spec$chr, source_chr)
    if (length(unknown)) stopf("Track `%s` selects unknown source chromosomes: %s.",
      id, format_chr_rows(unknown))
    selected <- spec$chr
    labels <- spec$.labels_source
    if (!is.null(labels)) {
      selected <- spec$.labels_chrs %||%
        as.character(mapped_fields(labels$data, labels$mapping, "chr", what = "mapping")$chr)
    }
    intersect(visible_chr, selected %||% visible_chr)
  })
  names(owners) <- names(tracks)
  rows <- vector("list", length(tracks))
  half <- chromosome_width / 2
  first <- half + marker_extent
  cursors <- lapply(as.list(first), rep, length.out = if (compact) length(visible_chr) else 1L)
  for (index in seq_along(tracks)) {
    spec <- tracks[[index]]
    gap <- spec$gap
    if (spec$side == "overlay") {
      if (abs(spec$offset) + spec$width / 2 > half + sqrt(.Machine$double.eps)) {
        stopf(paste0(
          "Overlay track `%s` extends beyond the chromosome body.\n",
          "  Adjust its `width` or `offset` to fit inside body width %s."),
          names(tracks)[index], format(chromosome_width))
      }
      low <- spec$offset - spec$width / 2
      high <- spec$offset + spec$width / 2
    } else {
      gap <- max(gap, min_gap)
      side <- spec$side
      keys <- owners[[index]]
      slots <- if (compact) match(keys, visible_chr) else if (length(keys)) 1L else integer()
      near <- max(c(first[[side]], cursors[[side]][slots])) + gap
      far <- near + spec$width
      sign <- if (side == "left") -1 else 1
      low <- sign * near
      high <- sign * far
      cursors[[side]][slots] <- far
    }
    rows[[index]] <- data.frame(
      id = names(tracks)[index], side = spec$side,
      width = spec$width, gap = gap,
      value_scale = spec$value_scale, reverse = spec$reverse,
      transform_name = spec$transform_name, offset = spec$offset,
      low_offset = low, high_offset = high,
      stringsAsFactors = FALSE
    )
    rows[[index]]$limits <- I(list(spec$limits))
    rows[[index]]$transform <- I(list(spec$transform))
  }
  table <- do.call(rbind, rows)
  row.names(table) <- NULL
  validate_global_track_specs(table)
  validate_overlay_tracks(table, owners)
  list(
    table = table, owners = owners,
    left = max(cursors$left) - half,
    right = max(cursors$right) - half
  )
}

validate_overlay_tracks <- function(table, owners) {
  rows <- which(table$side == "overlay" & lengths(owners) > 0L)
  if (length(rows) < 2L) return(invisible(NULL))
  for (j in seq.int(2L, length(rows))) {
    i <- rows[j]
    others <- rows[seq_len(j - 1L)]
    width <- pmin(table$high_offset[i], table$high_offset[others]) -
      pmax(table$low_offset[i], table$low_offset[others])
    for (other in others[width > sqrt(.Machine$double.eps)]) {
      shared <- intersect(owners[[i]], owners[[other]])
      if (length(shared)) stopf(paste0(
        "Overlay tracks `%s` and `%s` overlap on %s. ",
        "Adjust `width` or `offset`, or use one track for layers that should overlap."),
        table$id[other], table$id[i], format_chr_rows(shared))
    }
  }
}

track_chromosomes <- function(layout, id) {
  intersect(layout$chrom$.chr, layout$track_owners[[id]] %||% layout$chrom$.chr)
}

check_track_chromosomes <- function(layout, id, chr) {
  outside <- setdiff(chr, track_chromosomes(layout, id))
  if (length(outside)) stopf("Track `%s` is not assigned to chromosome%s %s.",
    id, if (length(outside) > 1L) "s" else "", format_chr_rows(outside))
}

validate_global_track_specs <- function(table) {
  global <- table[table$value_scale == "global", , drop = FALSE]
  if (nrow(global) < 2L) return(invisible(table))
  if (length(unique(global$transform_name)) != 1L) {
    stopf("Tracks sharing `value_scale = \"global\"` must use one transform.")
  }
  explicit <- Filter(Negate(is.null), global$limits)
  if (length(explicit) > 1L &&
      any(!vapply(explicit[-1], identical, logical(1), explicit[[1]]))) {
    stopf("Global tracks with explicit `limits` must use identical limits.")
  }
  invisible(table)
}

initialize_track_ranges <- function(track_table, chromosomes) {
  registry <- list()
  if (!nrow(track_table)) return(registry)
  for (index in seq_len(nrow(track_table))) {
    spec <- track_table[index, , drop = FALSE]
    keys <- track_range_keys(spec, chromosomes)
    limits <- spec$limits[[1]]
    for (key in keys) {
      existing <- registry[[key]]
      if (is.null(existing)) {
        registry[[key]] <- list(limits = limits, fixed = !is.null(limits))
      } else if (!is.null(limits)) {
        if (is.null(existing$limits)) {
          registry[[key]] <- list(limits = limits, fixed = TRUE)
        } else if (!identical(existing$limits, limits)) {
          stopf("Shared track range `%s` has inconsistent explicit limits.", key)
        }
      }
    }
  }
  registry
}

track_range_keys <- function(spec, chromosomes) {
  scope <- spec$value_scale[[1]]
  if (scope == "global") return("global")
  if (scope == "per_track") return(paste0("track:", spec$id[[1]]))
  paste0("track:", spec$id[[1]], ":chr:", chromosomes)
}

track_table_row <- function(layout, track_id) {
  if (!is.character(track_id) || length(track_id) != 1L ||
      is.na(track_id) || !nzchar(track_id)) {
    stopf("`track` must be one non-empty track identifier.")
  }
  index <- match(track_id, layout$tracks$id)
  if (is.na(index)) {
    available <- if (nrow(layout$tracks)) {
      paste0("`", layout$tracks$id, "`", collapse = ", ")
    } else {
      "none"
    }
    stopf("Unknown track `%s`; declared tracks: %s.", track_id, available)
  }
  layout$tracks[index, , drop = FALSE]
}

track_key_for_rows <- function(spec, chr) {
  scope <- spec$value_scale[[1]]
  if (scope == "global") return(rep("global", length(chr)))
  if (scope == "per_track") {
    return(rep(paste0("track:", spec$id[[1]]), length(chr)))
  }
  paste0("track:", spec$id[[1]], ":chr:", chr)
}

register_track_values <- function(layout, track_id, chr, value) {
  check_layout_v2(layout)
  spec <- track_table_row(layout, track_id)
  chr <- as.character(chr)
  match_layout_chr(layout, chr)
  check_track_chromosomes(layout, track_id, chr)
  if (!is.numeric(value) || length(value) != length(chr)) {
    stopf("Track `value` must be numeric and match the chromosome rows.")
  }
  finite <- !is.na(value)
  if (any(!is.finite(value[finite]))) {
    stopf("Track `value` must contain only finite numbers or missing values.")
  }
  if (!any(finite)) return(layout)

  transformed <- spec$transform[[1]]$transform(value[finite])
  if (anyNA(transformed) || any(!is.finite(transformed))) {
    stopf("Track `%s` has values outside the `%s` transform domain.",
          track_id, spec$transform_name)
  }
  keys <- track_key_for_rows(spec, chr)
  registry <- layout$track_ranges
  for (key in unique(keys[finite])) {
    local <- value[finite & keys == key]
    entry <- registry[[key]]
    if (is.null(entry)) {
      entry <- list(limits = NULL, fixed = FALSE)
    }
    if (isTRUE(entry$fixed)) {
      outside <- local < entry$limits[1] | local > entry$limits[2]
      if (any(outside)) {
        stopf(paste0(
          "Track `%s` has value%s outside explicit limits [%s, %s]: %s."),
          track_id, if (sum(outside) > 1L) "s" else "",
          format(entry$limits[1]), format(entry$limits[2]),
          paste(format(unique(local[outside]), trim = TRUE), collapse = ", "))
      }
    } else {
      entry$limits <- range(c(entry$limits, local), na.rm = TRUE)
    }
    registry[[key]] <- entry
  }
  layout$track_ranges <- registry
  layout
}

track_value_fraction <- function(layout, track_id, chr, value) {
  spec <- track_table_row(layout, track_id)
  check_track_chromosomes(layout, track_id, chr)
  keys <- track_key_for_rows(spec, as.character(chr))
  fraction <- numeric(length(value))
  for (key in unique(keys)) {
    index <- which(keys == key)
    entry <- layout$track_ranges[[key]]
    if (is.null(entry) || is.null(entry$limits)) {
      stopf("Track `%s` has no resolved value range for `%s`.", track_id, key)
    }
    transformed_limits <- spec$transform[[1]]$transform(entry$limits)
    transformed_value <- spec$transform[[1]]$transform(value[index])
    if (anyNA(transformed_value) || any(!is.finite(transformed_value))) {
      stopf("Track `%s` has values outside the `%s` transform domain.",
            track_id, spec$transform_name)
    }
    tolerance <- sqrt(.Machine$double.eps) *
      max(1, abs(transformed_limits), abs(transformed_value))
    outside <- transformed_value < transformed_limits[1] - tolerance |
      transformed_value > transformed_limits[2] + tolerance
    if (any(outside)) {
      stopf("Track `%s` geometry extends beyond its resolved value limits.",
            track_id)
    }
    span <- diff(transformed_limits)
    fraction[index] <- if (span == 0) {
      0.5
    } else {
      (transformed_value - transformed_limits[1]) / span
    }
  }
  if (isTRUE(spec$reverse)) fraction <- 1 - fraction
  fraction
}

project_track_values_raw <- function(layout, track_id, chr, position, value,
                                     check_position = TRUE) {
  if (!is.numeric(position) || anyNA(position) || any(!is.finite(position))) {
    stopf("Track `position` must contain finite numeric coordinates.")
  }
  index <- match_layout_chr(layout, chr)
  point <- if (check_position) {
    project_positions_checked(layout, chr, position)
  } else {
    point <- project_positions_raw(layout$chrom, chr, position)
    point$index <- index
    point
  }
  fraction <- track_value_fraction(layout, track_id, chr, value)
  spec <- track_table_row(layout, track_id)
  offset <- spec$low_offset + fraction * (spec$high_offset - spec$low_offset)
  g <- layout$chrom[index, , drop = FALSE]
  normal <- chromosome_right_normal(layout, g)
  list(
    x = point$x + normal$nx * offset,
    y = point$y + normal$ny * offset,
    index = index,
    fraction = fraction,
    offset = offset
  )
}

#' Project chromosome track observations
#'
#' This is the public low-level bridge for third-party ggplot2 extensions. The
#' track must be declared in the layout and its value range must already be
#' explicit or registered by a track component.
#'
#' @param layout An `ideogram_layout_v2` object or ggideogram plot.
#' @param data A data frame.
#' @param chr,position,value Column names in `data`.
#' @param track One declared track identifier.
#'
#' @return `data` with stable prefixed semantic and projected columns.
#' @export
project_chr_track <- function(layout, data, chr, position, value, track) {
  layout <- chr_projection_layout(layout)
  check_layout_v2(layout)
  check_projection_data(data)
  chr_value <- projection_column(data, chr, "chr")
  position_value <- projection_column(data, position, "position")
  value_value <- projection_column(data, value, "value")
  projected <- project_track_values_raw(
    layout, track, chr_value, position_value, value_value)
  data$.chr <- as.character(chr_value)
  data$.position <- as.numeric(position_value)
  data$.value <- as.numeric(value_value)
  data$.track <- track
  data$.x <- projected$x
  data$.y <- projected$y
  data$.chr_index <- projected$index
  data$.value_fraction <- projected$fraction
  data
}
