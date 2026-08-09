# Declarative track specifications for the dimensionless renderer.

#' Declare one chromosome track
#'
#' A track owns a stable transverse region beside (or over) every chromosome.
#' Its dimensions are normalized chromosome-body-width units, while point,
#' text and line sizes remain ordinary ggplot2 physical units.
#'
#' @param side Track placement: `"right"`, `"left"`, or `"overlay"`.
#' @param width Width of the value region in chromosome-body-width units.
#' @param gap Clear distance from the chromosome body or preceding track. If
#'   `NULL`, the documented default is `0.2` beside a chromosome and `0` for an
#'   overlay.
#' @param value_scale Scope of an automatically derived value range:
#'   `"global"` shares one range across all global tracks, `"per_track"`
#'   shares within this track, and `"per_chr"` resolves each chromosome
#'   separately.
#' @param limits Optional finite increasing raw-data limits. Explicit limits
#'   make a track independent of layer-addition order and are required by
#'   components whose standard ggplot2 Stat must run after global projection,
#'   including ribbon, area, boxplot and violin tracks.
#' @param transform A transformation name understood by
#'   [scales::as.transform()] or a scales transform object.
#' @param reverse Reverse the direction from low to high values.
#'
#' @return An `ideogram_track_spec` object, to be named inside
#'   [track_layout()].
#' @export
track <- function(
    side = c("right", "left", "overlay"),
    width = 1,
    gap = NULL,
    value_scale = c("per_track", "global", "per_chr"),
    limits = NULL,
    transform = "identity",
    reverse = FALSE) {
  side <- match.arg(side)
  value_scale <- match.arg(value_scale)
  gap <- gap %||% if (side == "overlay") 0 else 0.2
  check_positive_layout(width, "width")
  check_nonnegative_layout(gap, "gap")
  if (side == "overlay" && gap != 0) {
    stopf("An overlay track must use `gap = 0`.")
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
      reverse = reverse
    ),
    class = "ideogram_track_spec"
  )
}

#' Combine named chromosome tracks
#'
#' Track names are stable identifiers used by every track layer. Tracks on the
#' same side are laid out in declaration order from the chromosome body
#' outwards; adding ggplot2 layers later never changes that order.
#'
#' @param ... Named [track()] objects.
#'
#' @return An `ideogram_track_layout` object.
#' @examples
#' tracks <- track_layout(
#'   markers = track(side = "right", width = 0.6),
#'   density = track(side = "right", width = 1.2, limits = c(0, 10)),
#'   counts = track(side = "left", width = 1, limits = c(0, 100))
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
    stopf("Track%s %s %s not created by `track()`.",
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

resolve_track_geometry <- function(tracks, chromosome_width) {
  if (is.null(tracks)) tracks <- track_layout()
  if (!inherits(tracks, "ideogram_track_layout")) {
    stopf("`tracks` must be NULL or created by `track_layout()`.")
  }
  if (!length(tracks)) {
    table <- data.frame(
      id = character(), side = character(), width = numeric(),
      gap = numeric(), value_scale = character(), reverse = logical(),
      transform_name = character(), low_offset = numeric(),
      high_offset = numeric(), stringsAsFactors = FALSE
    )
    table$limits <- I(list())
    table$transform <- I(list())
    return(list(table = table, left = 0, right = 0))
  }

  rows <- vector("list", length(tracks))
  left_cursor <- chromosome_width / 2
  right_cursor <- chromosome_width / 2
  for (index in seq_along(tracks)) {
    spec <- tracks[[index]]
    if (spec$side == "overlay") {
      if (spec$width > chromosome_width) {
        stopf(paste0(
          "Overlay track `%s` is wider than the chromosome body.\n",
          "  Reduce its `width` to at most %s or place it beside the body."),
          names(tracks)[index], format(chromosome_width))
      }
      low <- -spec$width / 2
      high <- spec$width / 2
    } else if (spec$side == "right") {
      low <- right_cursor + spec$gap
      high <- low + spec$width
      right_cursor <- high
    } else {
      low <- -(left_cursor + spec$gap)
      high <- low - spec$width
      left_cursor <- abs(high)
    }
    rows[[index]] <- data.frame(
      id = names(tracks)[index], side = spec$side,
      width = spec$width, gap = spec$gap,
      value_scale = spec$value_scale, reverse = spec$reverse,
      transform_name = spec$transform_name,
      low_offset = low, high_offset = high,
      stringsAsFactors = FALSE
    )
    rows[[index]]$limits <- I(list(spec$limits))
    rows[[index]]$transform <- I(list(spec$transform))
  }
  table <- do.call(rbind, rows)
  row.names(table) <- NULL
  validate_global_track_specs(table)
  list(
    table = table,
    left = max(0, left_cursor - chromosome_width / 2),
    right = max(0, right_cursor - chromosome_width / 2)
  )
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
  if (!is.numeric(value) || length(value) != length(chr)) {
    stopf("Track `value` must be numeric and match the chromosome rows.")
  }
  finite <- !is.na(value)
  if (any(!is.finite(value[finite]))) {
    stopf("Track `value` must contain only finite numbers or missing values.")
  }
  if (!any(finite)) {
    stopf("Track `value` contains no finite observations.")
  }

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
#' @param layout An `ideogram_layout_v2` object.
#' @param data A data frame.
#' @param chr,position,value Column names in `data`.
#' @param track One declared track identifier.
#'
#' @return `data` with stable prefixed semantic and projected columns.
#' @export
project_chr_track <- function(layout, data, chr, position, value, track) {
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
