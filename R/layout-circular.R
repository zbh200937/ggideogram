# Circular layouts retain an unfolded arc-length/radius data domain.

is_circular_layout <- function(layout) identical(layout$orientation, "circular")

canonical_chr_side <- function(side) {
  side[side == "inner"] <- "left"
  side[side == "outer"] <- "right"
  side
}

circular_ideogram_layout <- function(data, k, tracks, track_geometry,
    chromosome_width, radius, start_angle, gap_angle, opening_angle, clockwise,
    scale_length, order_by, genome_order, homolog_order, reverse_chr,
    curve_points) {
  check_positive_layout(radius, "radius")
  if (!is.null(start_angle) &&
      (!is.numeric(start_angle) || length(start_angle) != 1L || !is.finite(start_angle))) {
    stopf("`start_angle` must be one finite angle in degrees.")
  }
  if (!is.logical(clockwise) || length(clockwise) != 1L || is.na(clockwise)) {
    stopf("`clockwise` must be TRUE or FALSE.")
  }
  inner <- radius - chromosome_width / 2 - track_geometry$left
  if (inner <= 0) {
    stopf("`radius` must accommodate the chromosome body and all inner tracks.")
  }
  order <- chromosome_layout_slots(k, NULL, order_by, genome_order, homolog_order)$index
  k <- k[order, , drop = FALSE]
  n <- nrow(k)
  gap <- circular_sector_gaps(gap_angle, k$.chr, data$karyotype$.chr, opening_angle)
  if (is.null(start_angle)) start_angle <- if (clockwise) gap[n] else 0
  available <- (360 - sum(gap)) * pi / 180
  length <- k$.end - k$.start
  sweep <- available * if (scale_length == "global") length / sum(length) else rep(1 / n, n)
  arc_start <- radius * c(0, utils::head(cumsum(sweep + gap * pi / 180), -1L))
  arc_end <- arc_start + radius * sweep
  reverse_chr <- validate_chr(reverse_chr, "reverse_chr")
  unknown <- setdiff(reverse_chr, k$.chr)
  if (length(unknown)) stopf("Unknown reverse chromosome: %s.", format_chr_rows(unknown))
  reverse <- k$.chr %in% reverse_chr
  chrom <- data.frame(
    .chr = k$.chr, .start = k$.start, .end = k$.end,
    .centromere_start = k$.centromere_start, .centromere_end = k$.centromere_end,
    .length = length, .display_length = radius * sweep,
    .units_per_bp = if (scale_length == "global") rep(radius * available / sum(length), n) else radius * sweep / length,
    .row = 0L, .col = seq_len(n) - 1L,
    .axis_start_x = ifelse(reverse, arc_end, arc_start), .axis_start_y = radius,
    .axis_end_x = ifelse(reverse, arc_start, arc_end), .axis_end_y = radius,
    .x_min = arc_start, .x_max = arc_end,
    .y_min = radius - chromosome_width / 2, .y_max = radius + chromosome_width / 2,
    .arc_start = arc_start, .arc_end = arc_end, .circular_radius = radius,
    .direction = ifelse(reverse, -1, 1),
    .name_x = (arc_start + arc_end) / 2, .name_y = radius,
    stringsAsFactors = FALSE
  )
  if (!is.null(data$view)) {
    chrom$.source_start <- k$.source_start
    chrom$.source_end <- k$.source_end
    chrom$.cut_start <- chrom$.start > chrom$.source_start
    chrom$.cut_end <- chrom$.end < chrom$.source_end
  }
  fields <- intersect(c(".chr_name", ".genome", ".assembly", ".homolog_group", ".display_name"), names(k))
  chrom[fields] <- k[fields]
  c1 <- project_positions_raw(chrom, chrom$.chr, chrom$.centromere_start, allow_na = TRUE)
  c2 <- project_positions_raw(chrom, chrom$.chr, chrom$.centromere_end, allow_na = TRUE)
  chrom$.centromere_start_x <- c1$x; chrom$.centromere_start_y <- c1$y
  chrom$.centromere_end_x <- c2$x; chrom$.centromere_end_y <- c2$y
  outer <- radius + chromosome_width / 2 + track_geometry$right
  layout <- structure(list(
    data = data, chrom = chrom, index = stats::setNames(seq_len(n), chrom$.chr),
    ncol = n, nrow = 1L, row_height = max(chrom$.display_length),
    chromosome_width = chromosome_width, chromosome_gap = gap,
    row_gap = 0, max_chr_length = max(chrom$.display_length),
    curve_points = curve_points, orientation = "circular", scale_length = scale_length,
    order_by = order_by, length_comparable = identical(scale_length, "global"),
    track_layout = tracks %||% track_layout(), tracks = track_geometry$table,
    track_ranges = initialize_track_ranges(track_geometry$table, chrom$.chr),
    track_extent = c(left = track_geometry$left, right = track_geometry$right),
    attachments = list(chr = list(), locus = list()),
    bounds = list(x = c(-outer, outer), y = c(-outer, outer)),
    units = "chromosome_width",
    circular = list(radius = radius, start_angle = start_angle, clockwise = clockwise,
      gap_angle = stats::setNames(gap, chrom$.chr), opening_angle = gap[n],
      inner_radius = inner,
      perimeter = 2 * pi * radius,
      semantic_separation = 2 * (max(k$.end) - min(k$.start) + 2 * pi * radius + 1))
  ), class = c("ideogram_layout_v2", "ideogram_layout"))
  left_arc <- if (clockwise) arc_end[n] else 0
  left_angle <- circular_theta(layout, left_arc) * 180 / pi
  layout$circular$label_gap <- list(angle = gap[n],
    left_arc = left_arc, right_arc = if (clockwise) 0 else arc_end[n],
    left_angle = left_angle, right_angle = left_angle + gap[n])
  layout
}

circular_sector_gaps <- function(gap, chr, source_chr, opening = NULL) {
  if (is.null(gap)) gap <- stats::setNames(c(rep(2, length(chr) - 1L), 20), chr)
  if (!is.numeric(gap) || !length(gap) || anyNA(gap) || any(!is.finite(gap)) || any(gap < 0)) {
    stopf("`gap_angle` must contain finite non-negative angles.")
  }
  if (length(gap) == 1L && is.null(names(gap))) {
    gap <- rep(gap, length(chr))
  } else {
    if (is.null(names(gap)) || anyNA(names(gap)) || anyDuplicated(names(gap)) ||
        !all(chr %in% names(gap)) || any(!names(gap) %in% source_chr)) {
      stopf("Name every `gap_angle` by its chromosome key, or supply one unnamed angle.")
    }
    gap <- unname(gap[chr])
  }
  if (!is.null(opening)) {
    check_nonnegative_layout(opening, "opening_angle")
    gap[length(gap)] <- opening
  }
  if (sum(gap) >= 360) stopf("The sum of chromosome gaps and `opening_angle` must be less than 360 degrees.")
  gap
}

circular_theta <- function(layout, arc) {
  arc <- scales::squish_infinite(arc, range = c(0, layout$circular$perimeter))
  layout$circular$start_angle * pi / 180 +
    (if (layout$circular$clockwise) 1 else -1) * arc / layout$circular$radius
}

circular_xy <- function(layout, arc, radius) {
  theta <- circular_theta(layout, arc)
  radius <- scales::squish_infinite(radius, range = c(0, max(abs(unlist(layout$bounds)))))
  list(x = radius * sin(theta), y = radius * cos(theta))
}

circular_transform_data <- function(layout, data) {
  if (!nrow(data) || ("ideogram_cartesian" %in% names(data) && all(data$ideogram_cartesian))) return(data)
  theta <- circular_theta(layout, data$x)
  xy <- circular_xy(layout, data$x, data$y)
  data$x <- xy$x; data$y <- xy$y
  if (all(c("xend", "yend") %in% names(data))) {
    xy <- circular_xy(layout, data$xend, data$yend)
    data$xend <- xy$x; data$yend <- xy$y
  }
  if (all(c("nx", "ny") %in% names(data))) {
    direction <- if (layout$circular$clockwise) 1 else -1
    nx <- data$nx; ny <- data$ny
    data$nx <- direction * nx * cos(theta) + ny * sin(theta)
    data$ny <- -direction * nx * sin(theta) + ny * cos(theta)
  }
  data
}

include_ideogram_coordinates <- function(layout, x, y, xend = NULL, yend = NULL,
    cartesian = FALSE) {
  x <- c(x, xend); y <- c(y, yend)
  if (is_circular_layout(layout)) {
    radius <- if (cartesian) sqrt(x^2 + y^2) else abs(y)
    extent <- max(abs(unlist(layout$bounds)), radius, na.rm = TRUE)
    layout$bounds <- list(x = c(-extent, extent), y = c(-extent, extent))
    return(layout)
  }
  layout$bounds$x <- range(layout$bounds$x, x, na.rm = TRUE)
  layout$bounds$y <- range(layout$bounds$y, y, na.rm = TRUE)
  layout
}

circular_coordinate_distance <- function(layout, x, y) {
  dx <- abs(diff(x))
  semantic <- pmax(abs(utils::head(x, -1L)), abs(utils::tail(x, -1L))) >
    layout$circular$semantic_separation / 2
  semantic[is.na(semantic)] <- FALSE
  dx[semantic] <- dx[semantic] * max(layout$chrom$.units_per_bp)
  diameter <- 2 * max(abs(unlist(layout$bounds)))
  # Native coordinate munching subdivides long edges; the bound keeps a raw
  # track's value units from asking for excessive subdivisions.
  pmin(1, sqrt((dx / diameter)^2 + (diff(y) / diameter)^2))
}

StatCircularName <- ggplot2::ggproto("StatCircularName", ggplot2::StatIdentity,
  extra_params = c("na.rm", "name_gap", "name_size"),
  compute_layer = function(self, data, params, layout) {
    object <- layout$coord$layout
    offset <- object$chromosome_width / 2 + max(object$track_extent[["right"]],
      (object$marker_extent %||% c(right = 0))[["right"]])
    axis <- object$circular_name_axis
    if (!is.null(axis)) offset <- offset + axis$gap
    data$y <- object$circular$radius + offset
    angle <- -circular_theta(object, data$x) * 180 / pi
    angle <- (angle + 180) %% 360 - 180
    flip <- angle < -90 | angle > 90
    data$angle <- ifelse(flip, (angle + 360) %% 360 - 180, angle)
    data$hjust <- 0.5
    data$vjust <- 0.5
    data$nx <- 0; data$ny <- 1
    data
  })

circular_chromosome_name_layer <- function(layout, size, gap, family, colour,
                                           position = "auto") {
  x <- switch(position, start = layout$chrom$.arc_start,
    end = layout$chrom$.arc_end, layout$chrom$.name_x)
  data <- data.frame(x = x, y = layout$circular$radius,
    label = layout$chrom$.display_name %||% layout$chrom$.chr)
  ggplot2::layer(data = data,
    mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
    stat = StatCircularName, geom = GeomCircularName, position = "identity",
    params = list(name_gap = gap, name_size = size, size = size, family = family,
      colour = colour, na.rm = FALSE), inherit.aes = FALSE)
}

circular_text_dimensions <- function(label, size, family = "", fontface = 1) {
  grobs <- lapply(as.character(label), function(value) grid::textGrob(value,
    gp = grid::gpar(fontsize = size * ggplot2::.pt, fontfamily = family, fontface = fontface)))
  list(width = max(do.call(grid::unit.c, lapply(grobs, grid::grobWidth))),
       height = max(do.call(grid::unit.c, lapply(grobs, grid::grobHeight))))
}

GeomCircularName <- ggplot2::ggproto("GeomCircularName", ggplot2::GeomText,
  draw_panel = function(data, panel_params, coord, name_gap = 0.5,
      name_size = ideogram_text_size("chromosome"), na.rm = FALSE) {
    text <- ggplot2::GeomText$draw_panel(data, panel_params, coord, na.rm = na.rm)
    xy <- coord$transform(data, panel_params)
    dimensions <- chromosome_text_dimensions(data)
    offset <- max(dimensions$height) / 2 + grid::unit(name_gap * name_size, "mm")
    axis <- coord$layout$circular_name_axis
    if (!is.null(axis)) {
      offset <- offset + grid::unit(axis$tick_length, "mm")
      if (isTRUE(axis$labels) && length(axis$label)) {
        dimensions <- circular_text_dimensions(axis$label, axis$size, axis$family)
        theta <- circular_theta(coord$layout, axis$label_x %||% data$x)
        extent <- max(abs(sin(theta)) * dimensions$width + abs(cos(theta)) * dimensions$height)
        offset <- offset + grid::unit(axis$label_gap * axis$size, "mm") + extent
      }
    }
    text$x <- text$x + xy$nx * offset
    text$y <- text$y + xy$ny * offset
    text
  })

reserve_circular_text_space <- function(plot) {
  layout <- ideogram_plot_layout(plot)
  if (!is_circular_layout(layout)) return(plot)
  needed <- grid::unit(rep(0, 4), "mm")
  axis_offset <- grid::unit(0, "mm")
  axis <- layout$circular_name_axis
  if (!is.null(axis)) {
    axis_offset <- grid::unit(axis$tick_length, "mm")
    if (isTRUE(axis$labels) && length(axis$label)) {
      dimensions <- circular_text_dimensions(axis$label, axis$size, axis$family)
      theta <- circular_theta(layout, axis$label_x %||% layout$chrom$.name_x)
      axis_offset <- axis_offset + max(abs(sin(theta)) * dimensions$width +
        abs(cos(theta)) * dimensions$height) +
        grid::unit(axis$label_gap * axis$size, "mm")
    }
  }
  for (layer in plot$layers) {
    if (!inherits(layer$geom, "GeomCircularName")) next
    data <- layer$data
    data$size <- layer$aes_params$size
    data$family <- layer$aes_params$family %||% ""
    data$fontface <- 1; data$lineheight <- 1.2
    dimensions <- chromosome_text_dimensions(data)
    theta <- circular_theta(layout, data$x)
    nx <- sin(theta); ny <- cos(theta)
    offset <- axis_offset + max(dimensions$height) / 2 +
      grid::unit(layer$geom_params$name_gap * layer$aes_params$size, "mm")
    half_x <- abs(cos(theta)) * dimensions$width / 2 + abs(sin(theta)) * dimensions$height / 2
    half_y <- abs(sin(theta)) * dimensions$width / 2 + abs(cos(theta)) * dimensions$height / 2
    overhang <- list(pmax(ny, 0) * offset + half_y,
      pmax(nx, 0) * offset + half_x, pmax(-ny, 0) * offset + half_y,
      pmax(-nx, 0) * offset + half_x)
    for (i in seq_len(4)) needed[i] <- max(needed[i], max(overhang[[i]]))
  }
  margin <- plot$theme$plot.margin %||% theme_ideogram()$plot.margin
  plot <- plot + ggplot2::theme(plot.margin = grid::unit.pmax(margin, needed))
  spacing <- plot$theme$legend.box.spacing
  # Only preset spacing grows automatically; explicit theme units keep priority.
  if (isTRUE(attr(spacing, "ideogram_auto_spacing"))) {
    spacing <- grid::unit.pmax(spacing, needed[3])
    attr(spacing, "ideogram_auto_spacing") <- TRUE
    plot <- plot + ggplot2::theme(legend.box.spacing = spacing)
  }
  reserve_ideogram_heading_space(plot, needed[1])
}
