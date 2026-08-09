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
                                      from = 0, to = g$.display_length) {
  radius <- width / 2
  if (g$.display_length < width) {
    stopf(paste0(
      "Chromosome %s is shorter than `chromosome_width` in display units.\n",
      "  Increase `max_chr_length`, reduce `chromosome_width`, or use ",
      "`scale_length = \"per_chr\"`."), format_chr_rows(g$.chr))
  }
  cap_start <- seq(0, radius, length.out = curve_points)
  cap_end <- seq(g$.display_length - radius, g$.display_length,
                 length.out = curve_points)
  waist <- numeric(0)
  if (!is.na(g$.centromere_start) &&
      g$.centromere_end > g$.centromere_start) {
    c1 <- (g$.centromere_start - g$.start) * g$.units_per_bp
    c2 <- (g$.centromere_end - g$.start) * g$.units_per_bp
    waist <- c(c1, mean(c(c1, c2)), c2)
  }
  values <- sort(unique(c(from, to, cap_start, cap_end, waist)))
  values[values >= from & values <= to]
}

chromosome_halfwidth_next <- function(g, s, width) {
  radius <- width / 2
  total <- g$.display_length
  halfwidth <- rep(radius, length(s))

  top <- s < radius
  halfwidth[top] <- sqrt(pmax(0, radius^2 - (radius - s[top])^2))
  bottom <- s > total - radius
  halfwidth[bottom] <- sqrt(pmax(
    0, radius^2 - (s[bottom] - (total - radius))^2))

  if (!is.na(g$.centromere_start) &&
      g$.centromere_end > g$.centromere_start) {
    c1 <- (g$.centromere_start - g$.start) * g$.units_per_bp
    c2 <- (g$.centromere_end - g$.start) * g$.units_per_bp
    at_waist <- s >= c1 & s <= c2
    halfwidth[at_waist] <- radius *
      abs(2 * (s[at_waist] - c1) / (c2 - c1) - 1)
  }
  halfwidth
}

chromosome_section_polygon <- function(g, width, curve_points,
                                       start = g$.start, end = g$.end) {
  basis <- chromosome_axis_basis(g)
  from <- (start - g$.start) * g$.units_per_bp
  to <- (end - g$.start) * g$.units_per_bp
  s <- chromosome_profile_breaks(g, width, curve_points, from, to)
  halfwidth <- chromosome_halfwidth_next(g, s, width)
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
    polygon$.fill <- fill[index]
    polygon
  })
  do.call(rbind, parts)
}

chromosome_body_layers <- function(layout, fill, colour, linewidth,
                                   curve_points, cytoband_scheme,
                                   cytoband_palette, cytoband_bleach) {
  layers <- list(chromosome_fill_layer(layout, fill, curve_points))
  bands <- chromosome_cytoband_layer(
    layout, curve_points, cytoband_scheme,
    cytoband_palette, cytoband_bleach)
  if (!is.null(bands)) layers <- c(layers, list(bands))
  c(layers, list(chromosome_outline_layer(
    layout, colour, linewidth, curve_points)))
}

chromosome_fill_layer <- function(layout, fill, curve_points) {
  body <- chromosome_polygon_data(layout, curve_points)
  ggplot2::geom_polygon(
    data = body,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
                           group = .data$.group),
    fill = fill, colour = NA, inherit.aes = FALSE
  )
}

chromosome_cytoband_layer <- function(
    layout, curve_points, scheme, palette, bleach) {
  bands <- cytoband_polygon_data(
    layout, curve_points, scheme, palette, bleach)
  if (is.null(bands)) return(NULL)
  ggplot2::geom_polygon(
    data = bands,
    mapping = ggplot2::aes(
      x = .data$x, y = .data$y, group = .data$.band,
      fill = I(.data$.fill), colour = I(.data$.fill)
    ),
    linewidth = 0, inherit.aes = FALSE, show.legend = FALSE
  )
}

chromosome_outline_layer <- function(
    layout, colour, linewidth, curve_points) {
  body <- chromosome_polygon_data(layout, curve_points)
  ggplot2::geom_path(
    data = body,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
                           group = .data$.group),
    colour = colour, linewidth = linewidth,
    linejoin = "round", lineend = "round", inherit.aes = FALSE
  )
}

chromosome_name_layer <- function(layout, size, gap, family, colour) {
  check_nonnegative_layout(gap, "name_gap")
  data <- data.frame(
    x = layout$chrom$.axis_end_x,
    y = layout$chrom$.axis_end_y,
    label = layout$chrom$.chr,
    stringsAsFactors = FALSE
  )
  vertical <- identical(layout$orientation, "vertical")
  ggplot2::geom_text(
    data = data,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
                           label = .data$label),
    hjust = if (vertical) 0.5 else -gap,
    vjust = if (vertical) 1 + gap else 0.5,
    size = size, family = family, colour = colour,
    inherit.aes = FALSE
  )
}

GeomIdeogramTick <- ggplot2::ggproto(
  "GeomIdeogramTick", ggplot2::Geom,
  required_aes = c("x", "y", "nx", "ny"),
  default_aes = ggplot2::aes(colour = "#666666", linewidth = 0.3,
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
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
      data, panel_params, coord, tick_length = 1.5, label_gap = 0.25,
      na.rm = FALSE) {
    coordinates <- coord$transform(data, panel_params)
    direction_length <- sqrt(coordinates$nx^2 + coordinates$ny^2)
    coordinates$nx <- coordinates$nx / direction_length
    coordinates$ny <- coordinates$ny / direction_length

    # `size` follows ggplot2's physical millimetre text convention.  Measuring
    # `label_gap` in em therefore keeps the requested clear space independent
    # of panel scaling, while `tick_length` starts that clearance at the tick
    # tip rather than at the axis spine.
    offset <- tick_length + label_gap * coordinates$size
    x <- grid::unit(coordinates$x, "native") +
      grid::unit(coordinates$nx * offset, "mm")
    y <- grid::unit(coordinates$y, "native") +
      grid::unit(coordinates$ny * offset, "mm")
    hjust <- ifelse(coordinates$nx > 0, 0,
                    ifelse(coordinates$nx < 0, 1, 0.5))
    vjust <- ifelse(coordinates$ny > 0, 0,
                    ifelse(coordinates$ny < 0, 1, 0.5))

    grid::textGrob(
      coordinates$label,
      x = x, y = y,
      hjust = hjust, vjust = vjust,
      rot = coordinates$angle,
      gp = grid::gpar(
        col = scales::alpha(coordinates$colour, coordinates$alpha),
        fontsize = coordinates$size * ggplot2::.pt,
        fontfamily = coordinates$family,
        fontface = coordinates$fontface,
        lineheight = coordinates$lineheight
      )
    )
  }
)

chromosome_axis_layers <- function(
    layout, chr, side, breaks, n, units, labels,
    gap, tick_length, label_gap, size, family,
    colour, linewidth) {
  if (identical(chr, FALSE)) return(list())
  selected <- if (isTRUE(chr)) layout$chrom$.chr else as.character(chr)
  unknown <- setdiff(selected, layout$chrom$.chr)
  if (length(unknown)) {
    stopf("Unknown axis chromosome%s: %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  side_sign <- if (side == "left") -1 else 1
  spine_parts <- list()
  tick_parts <- list()
  text_parts <- list()

  for (index in seq_along(selected)) {
    g <- layout$chrom[layout$index[[selected[index]]], , drop = FALSE]
    normal <- chromosome_right_normal(layout, g)
    span <- g$.end - g$.start
    local_breaks <- if (is.function(breaks)) {
      breaks(c(g$.start, g$.end))
    } else if (is.null(breaks)) {
      g$.start + axis_breaks(span, n)
    } else {
      breaks
    }
    local_breaks <- local_breaks[
      local_breaks >= g$.start & local_breaks <= g$.end]
    if (!length(local_breaks)) next
    projected <- project_positions_raw(layout$chrom, g$.chr, local_breaks)
    offset <- layout$chromosome_width / 2 + gap
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
                             nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity",
      geom = GeomIdeogramAxisText,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = tick_length, label_gap = label_gap,
        size = size, family = family, colour = colour, na.rm = FALSE)
    )))
  }
  layers
}
