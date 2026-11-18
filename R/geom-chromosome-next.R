# Dimensionless chromosome geometry for the next-generation renderer.

validate_curve_points <- function(points) {
  if (!is.numeric(points) || length(points) != 1L || is.na(points) ||
      !is.finite(points) || points < 4 || points != floor(points)) {
    stopf("`curve_points` must be one whole number of at least 4.")
  }
  as.integer(points)
}

chromosome_axis_basis <- function(g) {
  dx <- g$.axis_end_x - g$.axis_start_x
  dy <- g$.axis_end_y - g$.axis_start_y
  length <- sqrt(dx^2 + dy^2)
  list(tx = dx / length, ty = dy / length,
       nx = -dy / length, ny = dx / length,
       length = length)
}

chromosome_profile_breaks <- function(g, width, curve_points,
                                      from = 0, to = g$.display_length,
                                      cap_ratio = 1) {
  radius <- width / 2 * cap_ratio
  circular <- !is.null(g$.circular_radius)
  if (!circular && g$.display_length < width) {
    stopf(paste0(
      "Chromosome %s is shorter than `chromosome_width` in display units.\n",
      "  Increase `max_chr_length`, reduce `chromosome_width`, or use ",
      "`scale_length = \"per_chr\"`."), format_chr_rows(g$.chr))
  }
  cap_start <- seq(0, radius, length.out = curve_points)
  cap_end <- seq(g$.display_length - radius, g$.display_length,
                 length.out = curve_points)
  if (isTRUE(g$.cut_start)) cap_start <- numeric(0)
  if (isTRUE(g$.cut_end)) cap_end <- numeric(0)
  arc <- if (circular) seq(from, to, length.out = curve_points) else numeric(0)
  if (circular) cap_start <- cap_end <- numeric(0)
  waist <- numeric(0)
  if (!is.na(g$.centromere_start) &&
      g$.centromere_end > g$.centromere_start) {
    c1 <- (g$.centromere_start - g$.start) * g$.units_per_bp
    c2 <- (g$.centromere_end - g$.start) * g$.units_per_bp
    waist <- c(c1, mean(c(c1, c2)), c2)
  }
  values <- sort(unique(c(from, to, cap_start, cap_end, waist, arc)))
  values[values >= from & values <= to]
}

chromosome_halfwidth_next <- function(g, s, width, cap_ratio = 1) {
  radius <- width / 2
  cap <- radius * cap_ratio
  total <- g$.display_length
  halfwidth <- rep(radius, length(s))

  circular <- !is.null(g$.circular_radius)
  top <- !circular & !isTRUE(g$.cut_start) & s < cap
  halfwidth[top] <- radius * sqrt(pmax(0, 1 - (1 - s[top] / cap)^2))
  bottom <- !circular & !isTRUE(g$.cut_end) & s > total - cap
  halfwidth[bottom] <- sqrt(pmax(
    0, 1 - ((s[bottom] - (total - cap)) / cap)^2)) * radius

  if (!is.na(g$.centromere_start) &&
      g$.centromere_end > g$.centromere_start) {
    c1 <- (g$.centromere_start - g$.start) * g$.units_per_bp
    c2 <- (g$.centromere_end - g$.start) * g$.units_per_bp
    at_waist <- s >= c1 & s <= c2
    halfwidth[at_waist] <- pmin(halfwidth[at_waist],
      radius * abs(2 * (s[at_waist] - c1) / (c2 - c1) - 1))
  }
  halfwidth
}

chromosome_section_polygon <- function(g, width, curve_points,
                                       start = g$.start, end = g$.end,
                                       cap_ratio = 1) {
  basis <- chromosome_axis_basis(g)
  from <- (start - g$.start) * g$.units_per_bp
  to <- (end - g$.start) * g$.units_per_bp
  s <- chromosome_profile_breaks(g, width, curve_points, from, to, cap_ratio)
  halfwidth <- chromosome_halfwidth_next(g, s, width, cap_ratio)
  centre_x <- g$.axis_start_x + basis$tx * s
  centre_y <- g$.axis_start_y + basis$ty * s
  right <- data.frame(
    x = centre_x + basis$nx * halfwidth,
    y = centre_y + basis$ny * halfwidth
  )
  left <- data.frame(
    x = centre_x - basis$nx * halfwidth,
    y = centre_y - basis$ny * halfwidth
  )
  polygon <- rbind(right, left[rev(seq_len(nrow(left))), , drop = FALSE])
  rbind(polygon, polygon[1, , drop = FALSE])
}

chromosome_polygon_data <- function(layout, curve_points = 32) {
  check_layout_v2(layout)
  curve_points <- validate_curve_points(curve_points)
  parts <- lapply(seq_len(nrow(layout$chrom)), function(index) {
    g <- layout$chrom[index, , drop = FALSE]
    polygon <- chromosome_section_polygon(
      g, layout$chromosome_width, curve_points)
    polygon$.chr <- g$.chr
    polygon$.group <- index
    polygon
  })
  do.call(rbind, parts)
}

cytoband_polygon_data <- function(layout, curve_points = 32,
                                  scheme = "circos", palette = NULL,
                                  bleach = 0) {
  bands <- layout$data$cytoband
  if (is.null(bands) || !nrow(bands)) return(NULL)
  bands <- bands[bands$.chr %in% layout$chrom$.chr, , drop = FALSE]
  if (!nrow(bands)) return(NULL)
  bands$.start <- bands$.start - 1
  if (!is.null(layout$data$view)) {
    window <- layout$data$view
    bands$.start <- pmax(bands$.start, window$start)
    bands$.end <- pmin(bands$.end, window$end)
    bands <- bands[bands$.start < bands$.end, , drop = FALSE]
    if (!nrow(bands)) return(NULL)
  }
  curve_points <- validate_curve_points(curve_points)
  stain_column <- intersect(c("Stain", "gieStain"), names(bands))
  if (!length(stain_column)) {
    stopf("`cytoband` needs a `Stain` or `gieStain` column for drawing.")
  }
  fill <- cytoband_colours(
    bands[[stain_column[1]]], scheme = scheme,
    palette = palette, bleach = bleach)

  parts <- lapply(seq_len(nrow(bands)), function(index) {
    chromosome_index <- layout$index[[bands$.chr[index]]]
    g <- layout$chrom[chromosome_index, , drop = FALSE]
    polygon <- chromosome_section_polygon(
      g, layout$chromosome_width, curve_points,
      start = bands$.start[index], end = bands$.end[index])
    polygon$.chr <- bands$.chr[index]
    polygon$.band <- index
    polygon$.start <- bands$.start[index]
    polygon$.end <- bands$.end[index]
    polygon$.fill <- fill[index]
    polygon
  })
  do.call(rbind, parts)
}

chromosome_body_layers <- function(layout, fill, colour, linewidth,
                                   curve_points, cytoband_scheme,
                                   cytoband_palette, cytoband_bleach,
                                   alpha = NA, linetype = 1) {
  layers <- list(chromosome_fill_layer(layout, fill, curve_points, alpha))
  bands <- chromosome_cytoband_layer(
    layout, curve_points, cytoband_scheme,
    cytoband_palette, cytoband_bleach, alpha)
  if (!is.null(bands)) layers <- c(layers, list(bands))
  c(layers, list(chromosome_outline_layer(
    layout, colour, linewidth, curve_points, alpha, linetype)))
}

chromosome_fill_layer <- function(layout, fill, curve_points, alpha = NA) {
  body <- chromosome_polygon_data(layout, curve_points)
  ggplot2::layer(geom = GeomChromosomePolygon, stat = "identity", position = "identity",
    data = body,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
                           group = .data$.group, ideogram_chr = .data$.chr),
    params = list(fill = fill, colour = NA, alpha = alpha, curve_points = curve_points), inherit.aes = FALSE
  )
}

chromosome_cytoband_layer <- function(
    layout, curve_points, scheme, palette, bleach, alpha = NA) {
  bands <- cytoband_polygon_data(
    layout, curve_points, scheme, palette, bleach)
  if (is.null(bands)) return(NULL)
  ggplot2::layer(geom = GeomChromosomePolygon, stat = "identity", position = "identity",
    data = bands,
    mapping = ggplot2::aes(
      x = .data$x, y = .data$y, group = .data$.band,
      ideogram_chr = .data$.chr, ideogram_start = .data$.start,
      ideogram_end = .data$.end,
      fill = I(.data$.fill), colour = I(.data$.fill)
    ),
    params = list(linewidth = 0, alpha = alpha, curve_points = curve_points), inherit.aes = FALSE, show.legend = FALSE
  )
}

chromosome_outline_layer <- function(
    layout, colour, linewidth, curve_points, alpha = NA, linetype = 1) {
  body <- chromosome_polygon_data(layout, curve_points)
  ggplot2::layer(geom = GeomChromosomePath, stat = "identity", position = "identity",
    data = body,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
                           group = .data$.group, ideogram_chr = .data$.chr),
    params = list(colour = colour, linewidth = linewidth, alpha = alpha, linetype = linetype,
      linejoin = "round", lineend = "round", curve_points = curve_points), inherit.aes = FALSE
  )
}

chromosome_name_layer <- function(layout, size, gap, family, colour,
                                  position = "auto", align = "body") {
  check_nonnegative_layout(gap, "name_gap")
  position <- match.arg(position, c("auto", "start", "end", "middle"))
  if (is_circular_layout(layout)) return(circular_chromosome_name_layer(
    layout, size, gap, family, colour, position))
  data <- linear_chromosome_name_data(layout, position, align)
  layer <- ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
      label = .data$label, name_chr = .data$name_chr),
    stat = StatChromosomeName, geom = GeomChromosomeName,
    position = "identity",
    params = list(name_position = position, name_align = align, name_gap = gap,
      size = size, family = family, colour = colour, na.rm = FALSE),
    inherit.aes = FALSE
  )
  layer$ideogram_name_gap <- gap
  layer
}

linear_chromosome_name_data <- function(layout, position = "auto", align = "body") {
  g <- layout$chrom
  if (position == "auto") position <- "start"
  vertical <- layout$orientation == "vertical"
  left <- layout$chromosome_width / 2 + annotation_extent(layout, g$.chr, "left")
  right <- layout$chromosome_width / 2 + annotation_extent(layout, g$.chr, "right")
  data <- data.frame(name_chr = g$.chr,
    x = (g$.axis_start_x + g$.axis_end_x) / 2,
    y = (g$.axis_start_y + g$.axis_end_y) / 2,
    label = g$.display_name %||% g$.chr, nx = 0, ny = 0)
  if (vertical) {
    if (align == "tracks") data$x <- data$x + (right - left) / 2
    if (position == "middle") { data$x <- g$.axis_start_x + right; data$nx <- 1 }
    else {
      data$y <- if (position == "start") pmax(g$.axis_start_y, g$.axis_end_y) else
        pmin(g$.axis_start_y, g$.axis_end_y)
      data$ny <- if (position == "start") 1 else -1
    }
  } else {
    if (align == "tracks") data$y <- data$y + (left - right) / 2
    if (position == "middle") { data$y <- g$.axis_start_y - right; data$ny <- -1 }
    else {
      data$x <- if (position == "start") pmin(g$.axis_start_x, g$.axis_end_x) else
        pmax(g$.axis_start_x, g$.axis_end_x)
      data$nx <- if (position == "start") -1 else 1
    }
  }
  data
}

StatChromosomeName <- ggplot2::ggproto("StatChromosomeName", ggplot2::StatIdentity,
  required_aes = c("name_chr", "label"),
  extra_params = c("na.rm", "name_position", "name_align", "name_gap"),
  compute_layer = function(self, data, params, layout) {
    anchors <- linear_chromosome_name_data(layout$coord$layout, params$name_position,
      params$name_align)
    i <- match(data$name_chr, anchors$name_chr)
    data[c("x", "y", "nx", "ny")] <- anchors[i, c("x", "y", "nx", "ny")]
    data$hjust <- data$vjust <- 0.5
    data
  })

GeomChromosomeName <- ggplot2::ggproto("GeomChromosomeName", ggplot2::GeomText,
  draw_panel = function(data, panel_params, coord, name_gap = 0.5,
      name_axis_clearance = NULL, name_margin_clearance = NULL,
      name_end_gap = NULL, name_title_alignment = NULL, na.rm = FALSE) {
    text <- ggplot2::GeomText$draw_panel(data, panel_params, coord, na.rm = na.rm)
    xy <- coord$transform(data, panel_params)
    dimensions <- chromosome_text_dimensions(data)
    gap <- grid::unit(name_gap * data$size, "mm")
    offset <- abs(xy$nx) * dimensions$width / 2 + abs(xy$ny) * dimensions$height / 2 + gap
    if (!is.null(name_axis_clearance)) offset <- offset + name_axis_clearance
    if (any(name_end_gap > 0)) {
      end <- coord$transform(data.frame(x = data$x + data$nx * name_end_gap,
        y = data$y + data$ny * name_end_gap), panel_params)
      text$x <- text$x + grid::unit(end$x - xy$x, "native")
      text$y <- text$y + grid::unit(end$y - xy$y, "native")
    }
    text$x <- text$x + xy$nx * offset
    text$y <- text$y + xy$ny * offset
    linear_name_alignment_grob(text, data, xy, dimensions, coord, panel_params, gap,
      name_title_alignment)
  })

GeomIdeogramTick <- ggplot2::ggproto(
  "GeomIdeogramTick", ggplot2::Geom,
  required_aes = c("x", "y", "nx", "ny"),
  default_aes = ggplot2::aes(colour = ideogram_colour("axis"), linewidth = ideogram_linewidth("axis"),
                             alpha = NA),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord, tick_length = 1.5,
                        na.rm = FALSE) {
    coordinates <- coord$transform(data, panel_params)
    grid::segmentsGrob(
      x0 = grid::unit(coordinates$x, "native"),
      y0 = grid::unit(coordinates$y, "native"),
      x1 = grid::unit(coordinates$x, "native") +
        grid::unit(coordinates$nx * tick_length, "mm"),
      y1 = grid::unit(coordinates$y, "native") +
        grid::unit(coordinates$ny * tick_length, "mm"),
      gp = grid::gpar(
        col = scales::alpha(coordinates$colour, coordinates$alpha),
        lwd = coordinates$linewidth * ggplot2::.pt,
        lineend = "butt"
      )
    )
  }
)

GeomIdeogramAxisText <- ggplot2::ggproto(
  "GeomIdeogramAxisText", ggplot2::GeomText,
  required_aes = c("x", "y", "label", "nx", "ny"),
  default_aes = utils::modifyList(ggplot2::GeomText$default_aes,
    ggplot2::aes(axis_chr = NA_character_)),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
      data, panel_params, coord, tick_length = 1.5, label_gap = 0.25,
      check_overlap = FALSE, na.rm = FALSE, text_clearance = NULL,
      track_axis_width = NULL, title_hjust = NULL, title_vjust = NULL) {
    coordinates <- coord$transform(data, panel_params)
    direction_length <- sqrt(coordinates$nx^2 + coordinates$ny^2)
    coordinates$nx <- coordinates$nx / direction_length
    coordinates$ny <- coordinates$ny / direction_length

    # `size` follows ggplot2's physical millimetre text convention.  Measuring
    # `label_gap` in em therefore keeps the requested clear space independent
    # of panel scaling, while `tick_length` starts that clearance at the tick
    # tip rather than at the axis spine.
    offset <- grid::unit(tick_length + label_gap * coordinates$size, "mm") +
      (text_clearance %||% grid::unit(0, "mm"))
    rotated <- coordinates$angle %% 360 != 0
    if (any(rotated)) {
      dimensions <- chromosome_text_dimensions(coordinates)
      radians <- coordinates$angle * pi / 180
      # Center a rotated label, keeping its nearest edge outside the tick tip.
      nx <- coordinates$nx; ny <- coordinates$ny
      offset <- offset + as.numeric(rotated) * (abs(nx * cos(radians) + ny * sin(radians)) *
        dimensions$width + abs(ny * cos(radians) - nx * sin(radians)) * dimensions$height) / 2
    }
    x <- grid::unit(coordinates$x, "native") +
      coordinates$nx * offset
    y <- grid::unit(coordinates$y, "native") +
      coordinates$ny * offset
    hjust <- ifelse(coordinates$nx > 0, 0,
                    ifelse(coordinates$nx < 0, 1, 0.5))
    vjust <- ifelse(coordinates$ny > 0, 0,
                    ifelse(coordinates$ny < 0, 1, 0.5))
    if (!is.null(title_hjust)) hjust <- rep_len(title_hjust, nrow(coordinates))
    if (!is.null(title_vjust)) vjust <- rep_len(title_vjust, nrow(coordinates))
    hjust[rotated] <- vjust[rotated] <- 0.5

    text <- grid::textGrob(
      coordinates$label,
      x = x, y = y,
      hjust = hjust, vjust = vjust,
      check.overlap = check_overlap,
      rot = coordinates$angle,
      gp = grid::gpar(
        col = scales::alpha(coordinates$colour, coordinates$alpha),
        fontsize = coordinates$size * ggplot2::.pt,
        fontfamily = coordinates$family,
        fontface = coordinates$fontface,
        lineheight = coordinates$lineheight
      )
    )
    if (!is.null(track_axis_width)) {
      return(adaptive_track_axis_grob(text, coordinates, coord, panel_params,
        track_axis_width, label_gap))
    }
    layout <- coord$layout
    if (is.null(layout) || is_circular_layout(layout) || all(is.na(data$axis_chr))) return(text)
    axis_text_clearance_grob(text, data, coord, panel_params, check_overlap)
  }
)

chromosome_axis_layers <- function(
    layout, chr, side, breaks, n, units, labels,
    gap, tick_length, label_gap, size, family,
    colour, linewidth, angle = 0, check_overlap = NULL) {
  if (!is.numeric(angle) || length(angle) != 1L || !is.finite(angle)) {
    stopf("Axis label angle must be one finite number in degrees.")
  }
  if (!is.null(check_overlap) && (!is.logical(check_overlap) ||
      length(check_overlap) != 1L || is.na(check_overlap))) {
    stopf("Axis overlap checking must be TRUE, FALSE or NULL.")
  }
  if (identical(chr, FALSE)) return(list())
  selected <- if (isTRUE(chr)) layout$chrom$.chr else as.character(chr)
  unknown <- setdiff(selected, layout$chrom$.chr)
  if (length(unknown)) {
    stopf("Unknown axis chromosome%s: %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  shared <- can_share_microsynteny_axis(layout) && side == "auto" && isTRUE(chr) && is.null(breaks)
  if (shared) selected <- utils::tail(layout$chrom$.chr, 1)
  spine_parts <- list()
  tick_parts <- list()
  text_parts <- list()

  for (index in seq_along(selected)) {
    g <- layout$chrom[layout$index[[selected[index]]], , drop = FALSE]
    axis_side <- if (side == "auto") microsynteny_side(layout, g$.chr, "axis") else side
    side_sign <- if (axis_side == "left") -1 else 1
    normal <- chromosome_right_normal(layout, g)
    if (shared) {
      g$.start <- 0
      g$.end <- max(layout$chrom$.end - layout$chrom$.start)
      g$.axis_start_x <- min(layout$chrom$.x_min)
      g$.axis_end_x <- max(layout$chrom$.x_max)
      g$.direction <- 1
    }
    span <- g$.end - g$.start
    local_breaks <- if (is.function(breaks)) {
      breaks(c(g$.start, g$.end))
    } else if (is.null(breaks)) {
      if (is_circular_layout(layout)) g$.start + axis_breaks(span, n) else
        axis_breaks(span, n, start = g$.start)
    } else {
      breaks
    }
    if (!is.null(breaks) && (!is.numeric(local_breaks) || anyNA(local_breaks) || any(!is.finite(local_breaks)))) {
      stopf("Axis `breaks` must contain finite numeric genomic coordinates.")
    }
    local_breaks <- sort(unique(local_breaks[
      local_breaks >= g$.start & local_breaks <= g$.end]))
    if (is_circular_layout(layout) && is.null(breaks)) {
      local_breaks <- local_breaks[local_breaks < g$.end]
    }
    if (!length(local_breaks)) next
    projected <- project_positions_raw(if (shared) g else layout$chrom, g$.chr, local_breaks)
    ext <- if (isTRUE(layout$microsynteny)) {
      annotation_extent(layout, g$.chr, axis_side)
    } else {
      marker_extent <- (layout$marker_extent %||% c(left = 0, right = 0))[[axis_side]]
      max(annotation_extent(layout, g$.chr, axis_side), marker_extent)
    }
    offset <- layout$chromosome_width / 2 + ext + gap
    if (is_circular_layout(layout) && axis_side == "left" && offset >= layout$circular$radius) {
      stopf("The inner bp axis reaches the circle center. Increase `radius` or reduce `axis_gap`.")
    }
    nx <- normal$nx * side_sign
    ny <- normal$ny * side_sign
    tick_parts[[index]] <- data.frame(
      x = projected$x + nx * offset,
      y = projected$y + ny * offset,
      nx = nx, ny = ny,
      stringsAsFactors = FALSE
    )
    spine_parts[[index]] <- data.frame(
      x = g$.axis_start_x + nx * offset,
      y = g$.axis_start_y + ny * offset,
      xend = g$.axis_end_x + nx * offset,
      yend = g$.axis_end_y + ny * offset,
      stringsAsFactors = FALSE
    )
    step <- if (length(local_breaks) > 1L) min(diff(local_breaks)) else span
    unit <- axis_unit(step, units)
    text_parts[[index]] <- data.frame(
      x = projected$x + nx * offset,
      y = projected$y + ny * offset,
      nx = nx, ny = ny,
      axis_chr = g$.chr,
      label = axis_labels(local_breaks, unit$div, unit$suffix),
      stringsAsFactors = FALSE
    )
  }

  spine <- do.call(rbind, spine_parts)
  ticks <- do.call(rbind, tick_parts)
  text <- do.call(rbind, text_parts)
  if (is.null(ticks) || !nrow(ticks)) return(list())

  layers <- list(ggplot2::geom_segment(
    data = spine,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
                           xend = .data$xend, yend = .data$yend),
    colour = colour, linewidth = linewidth, inherit.aes = FALSE
  ))
  if (tick_length > 0) {
    layers <- c(layers, list(ggplot2::layer(
      data = ticks,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                             nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity", geom = GeomIdeogramTick,
      inherit.aes = FALSE,
      params = list(tick_length = tick_length, colour = colour,
                    linewidth = linewidth, na.rm = FALSE)
    )))
  }
  if (labels) {
    layers <- c(layers, list(ggplot2::layer(
      data = text,
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                             label = .data$label,
                             axis_chr = .data$axis_chr,
                             nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity",
      geom = GeomIdeogramAxisText,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = tick_length, label_gap = label_gap,
        check_overlap = check_overlap %||% is.null(breaks),
        size = size, angle = angle, family = family, colour = colour, na.rm = FALSE)
    )))
  }
  spec <- list(chr = chr, side = side, breaks = breaks, n = n, units = units,
    labels = labels, gap = gap, tick_length = tick_length, label_gap = label_gap,
    size = size, family = family, colour = colour, linewidth = linewidth,
    angle = angle, check_overlap = check_overlap)
  for (part in seq_along(layers)) {
    layers[[part]]$ideogram_axis_spec <- spec
    layers[[part]]$ideogram_axis_part <- part
  }
  layers
}
