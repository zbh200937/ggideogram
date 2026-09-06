# Coordinate bridge for the dimensionless renderer.

coord_ideogram_next <- function(layout, padding = 0.5, clip = "off") {
  check_layout_v2(layout)
  check_nonnegative_layout(padding, "padding")
  if (!clip %in% c("on", "off")) {
    stopf("`clip` must be \"on\" or \"off\".")
  }
  xlim <- layout$bounds$x + c(-padding, padding)
  ylim <- layout$bounds$y + c(-padding, padding)

  ggplot2::ggproto(
    NULL, ggplot2::CoordCartesian,
    limits = list(x = xlim, y = ylim),
    expand = FALSE,
    clip = clip,
    ratio = 1,
    reverse = "none",
    layout = layout,
    padding = padding,
    setup_panel_params = function(self, scale_x, scale_y, params = list()) {
      panel <- ggplot2::ggproto_parent(ggplot2::CoordCartesian, self)$setup_panel_params(
        scale_x, scale_y, params)
      # Each build and panel owns a fresh registry shared by its inset layers.
      panel$ideogram_insets <- new.env(parent = emptyenv())
      panel
    },
    transform = function(self, data, panel_params) {
      semantic <- "ideogram_track" %in% names(data) ||
        all(c("ideogram_chr", "ideogram_position") %in% names(data))
      if (semantic && nrow(data)) {
        data <- transform_ideogram_semantics(self$layout, data)
      }
      ggplot2::ggproto_parent(ggplot2::CoordCartesian, self)$transform(
        data, panel_params)
    }
  )
}

transform_ideogram_semantics <- function(layout, data) {
  if ("ideogram_track" %in% names(data)) {
    return(transform_ideogram_track(layout, data))
  }
  if (all(c("ideogram_inset_start", "ideogram_inset_end") %in% names(data))) {
    project_positions_checked(layout, data$ideogram_chr, data$ideogram_inset_start)
    project_positions_checked(layout, data$ideogram_chr, data$ideogram_inset_end)
  }
  data <- apply_chr_repel_request(layout, data)
  marker <- all(c("ideogram_side", "ideogram_gap") %in% names(data))
  marker_track <- marker && "ideogram_marker_track" %in% names(data)
  marker_offset <- if (marker_track) {
    marker_track_offsets(layout, data$ideogram_marker_track, data)
  } else {
    NULL
  }
  interval <- "ideogram_interval" %in% names(data) &&
    any(data$ideogram_interval %in% TRUE)
  if (interval) {
    first <- project_positions_checked(
      layout, data$ideogram_chr, data$ideogram_position)
    last <- project_positions_checked(
      layout, data$ideogram_chr, data$ideogram_end_position)
    beside <- all(data$ideogram_interval_placement == "beside")
    if (beside) {
      if (marker_track) {
        first <- offset_chr_points_signed(
          layout, data$ideogram_chr, first, marker_offset)
        last <- offset_chr_points_signed(
          layout, data$ideogram_chr, last, marker_offset)
      } else {
        distance <- layout$chromosome_width / 2 + data$ideogram_gap
        first <- offset_chr_points(
          layout, data$ideogram_chr, first,
          side = data$ideogram_side, distance = distance)
        last <- offset_chr_points(
          layout, data$ideogram_chr, last,
          side = data$ideogram_side, distance = distance)
      }
    }
    data$x <- first$x
    data$y <- first$y
    data$xend <- last$x
    data$yend <- last$y
    return(data)
  }
  link <- marker && "ideogram_link" %in% names(data) &&
    any(data$ideogram_link %in% TRUE)

  if (link) {
    anchor <- project_positions_checked(
      layout, data$ideogram_chr, data$ideogram_anchor_position)
    endpoint <- project_positions_checked(
      layout, data$ideogram_chr, data$ideogram_position)
    if (marker_track) {
      anchor_offset <- sign(marker_offset) * layout$chromosome_width / 2
      anchor_offset[marker_offset == 0] <- 0
      anchor <- offset_chr_points_signed(
        layout, data$ideogram_chr, anchor, anchor_offset)
      endpoint <- offset_chr_points_signed(
        layout, data$ideogram_chr, endpoint, marker_offset)
    } else {
      anchor <- offset_chr_points(
        layout, data$ideogram_chr, anchor,
        side = data$ideogram_side,
        distance = rep(layout$chromosome_width / 2, nrow(data)))
      endpoint <- offset_chr_points(
        layout, data$ideogram_chr, endpoint,
        side = data$ideogram_side,
        distance = layout$chromosome_width / 2 + data$ideogram_gap)
    }
    data$x <- anchor$x
    data$y <- anchor$y
    data$xend <- endpoint$x
    data$yend <- endpoint$y
    return(data)
  }

  projected <- project_positions_checked(
    layout, data$ideogram_chr, data$ideogram_position)
  if (marker) {
    projected <- if (marker_track) {
      offset_chr_points_signed(
        layout, data$ideogram_chr, projected, marker_offset)
    } else {
      offset_chr_points(
        layout, data$ideogram_chr, projected,
        side = data$ideogram_side,
        distance = layout$chromosome_width / 2 + data$ideogram_gap)
    }
  }
  data$x <- projected$x
  data$y <- projected$y
  data
}

marker_track_offsets <- function(layout, track_id, data = NULL) {
  track_id <- as.character(track_id)
  if (length(unique(track_id)) != 1L) {
    stopf("One marker layer can target only one declared track.")
  }
  spec <- track_table_row(layout, track_id[1])
  use_inner_edge <- !is.null(data) &&
    "ideogram_inset_track_edge" %in% names(data) &&
    all(data$ideogram_inset_track_edge %in% TRUE)
  offset <- if (use_inner_edge && spec$side != "overlay") {
    spec$low_offset
  } else {
    mean(c(spec$low_offset, spec$high_offset))
  }
  rep(offset, length(track_id))
}

offset_chr_points_signed <- function(layout, chr, point, offset) {
  if (!is.numeric(offset) || length(offset) != length(chr) ||
      anyNA(offset) || any(!is.finite(offset))) {
    stopf("Marker track offsets must be finite numbers.")
  }
  index <- match_layout_chr(layout, chr)
  g <- layout$chrom[index, , drop = FALSE]
  normal <- chromosome_right_normal(layout, g)
  point$x <- point$x + normal$nx * offset
  point$y <- point$y + normal$ny * offset
  point
}

transform_ideogram_track <- function(layout, data) {
  if (!all(c("x", "y", "ideogram_chr", "ideogram_track") %in%
           names(data))) {
    stopf("A chromosome track layer lost its required semantic coordinates.")
  }
  chr <- as.character(data$ideogram_chr)
  track_id <- as.character(data$ideogram_track)
  if (length(unique(track_id)) != 1L) {
    stopf("One ggplot2 layer can target only one chromosome track.")
  }
  track_id <- track_id[1]

  raw <- data
  shift <- raw$ideogram_position_shift %||% rep(0, nrow(raw))
  main <- project_track_coordinate_pair(
    layout, track_id, chr, raw$x - shift, raw$y,
    check_position = TRUE)
  data$x <- main$x
  data$y <- main$y

  if (all(c("xend", "yend") %in% names(raw))) {
    endpoint <- project_track_coordinate_pair(
      layout, track_id, chr, raw$xend - shift, raw$yend,
      check_position = TRUE)
    data$xend <- endpoint$x
    data$yend <- endpoint$y
  }

  rectangle <- all(c("xmin", "xmax", "ymin", "ymax") %in% names(raw))
  if (rectangle) {
    corners <- list(
      project_track_coordinate_pair(
        layout, track_id, chr, raw$xmin - shift, raw$ymin,
        check_position = FALSE),
      project_track_coordinate_pair(
        layout, track_id, chr, raw$xmin - shift, raw$ymax,
        check_position = FALSE),
      project_track_coordinate_pair(
        layout, track_id, chr, raw$xmax - shift, raw$ymin,
        check_position = FALSE),
      project_track_coordinate_pair(
        layout, track_id, chr, raw$xmax - shift, raw$ymax,
        check_position = FALSE)
    )
    x <- do.call(cbind, lapply(corners, `[[`, "x"))
    y <- do.call(cbind, lapply(corners, `[[`, "y"))
    data$xmin <- apply(x, 1, min)
    data$xmax <- apply(x, 1, max)
    data$ymin <- apply(y, 1, min)
    data$ymax <- apply(y, 1, max)
  }
  data
}

project_track_coordinate_pair <- function(
    layout, track_id, chr, position, value, check_position) {
  missing <- is.na(position) | is.na(value)
  result <- list(x = rep(NA_real_, length(chr)),
                 y = rep(NA_real_, length(chr)))
  if (all(missing)) return(result)
  projected <- project_track_values_raw(
    layout, track_id, chr[!missing], position[!missing], value[!missing],
    check_position = check_position)
  result$x[!missing] <- projected$x
  result$y[!missing] <- projected$y
  result
}

offset_chr_points <- function(layout, chr, point, side, distance) {
  side <- as.character(side)
  if (length(side) == 1L) side <- rep(side, length(chr))
  if (length(distance) == 1L) distance <- rep(distance, length(chr))
  if (length(side) != length(chr) || anyNA(side) ||
      any(!side %in% c("left", "right"))) {
    stopf("Marker `side` must contain only \"left\" or \"right\".")
  }
  if (!is.numeric(distance) || length(distance) != length(chr) ||
      anyNA(distance) || any(!is.finite(distance)) || any(distance < 0)) {
    stopf("Marker transverse distances must be finite non-negative numbers.")
  }

  index <- match_layout_chr(layout, chr)
  g <- layout$chrom[index, , drop = FALSE]
  normal <- chromosome_right_normal(layout, g)
  sign <- ifelse(side == "right", 1, -1)
  point$x <- point$x + normal$nx * sign * distance
  point$y <- point$y + normal$ny * sign * distance
  point
}

apply_chr_repel_request <- function(layout, data) {
  if (!"ideogram_repel_distance" %in% names(data)) return(data)
  required <- c("ideogram_anchor_position", "ideogram_repel_units")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stopf("Chromosome repel data is missing %s.",
          paste0("`", missing, "`", collapse = ", "))
  }

  distance <- data$ideogram_repel_distance
  units <- as.character(data$ideogram_repel_units)
  if (!is.numeric(distance) || anyNA(distance) ||
      any(!is.finite(distance)) || any(distance <= 0)) {
    stopf("Chromosome repel distances must be positive finite numbers.")
  }
  if (anyNA(units) || any(!units %in% c("bp", "body_width"))) {
    stopf("Chromosome repel units must be \"bp\" or \"body_width\".")
  }

  chr <- as.character(data$ideogram_chr)
  side <- if ("ideogram_side" %in% names(data)) {
    as.character(data$ideogram_side)
  } else {
    rep("overlay", nrow(data))
  }
  key <- interaction(chr, side, drop = TRUE, lex.order = TRUE)
  display <- as.numeric(data$ideogram_anchor_position)

  for (index in split(seq_len(nrow(data)), key)) {
    chromosome <- chr[index[1]]
    chromosome_index <- match_layout_chr(layout, chromosome)
    g <- layout$chrom[chromosome_index, , drop = FALSE]
    local_distance <- unique(distance[index])
    local_units <- unique(units[index])
    if (length(local_distance) != 1L || length(local_units) != 1L) {
      stopf("One repel rule must be used within each chromosome side.")
    }
    distance_bp <- if (local_units == "bp") {
      local_distance
    } else {
      local_distance / g$.units_per_bp
    }
    display[index] <- repel_1d(
      display[index], distance_bp, lo = g$.start, hi = g$.end)
  }
  data$ideogram_position <- display
  data
}
