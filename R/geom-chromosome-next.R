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

chromosome_halfwidth_next <- function(g, s, width) {
  radius <- width / 2
  total <- g$.display_length
  halfwidth <- rep(radius, length(s))

  circular <- !is.null(g$.circular_radius)
  top <- !circular & !isTRUE(g$.cut_start) & s < radius
  halfwidth[top] <- sqrt(pmax(0, radius^2 - (radius - s[top])^2))
  bottom <- !circular & !isTRUE(g$.cut_end) & s > total - radius
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
  bands <- bands[bands$.chr %in% layout$chrom$.chr, , drop = FALSE]
  if (!nrow(bands)) return(NULL)
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

chromosome_name_layer <- function(layout, size, gap, family, colour,
                                  position = "auto") {
  check_nonnegative_layout(gap, "name_gap")
  position <- match.arg(position, c("auto", "start", "end", "middle"))
  if (is_circular_layout(layout)) return(circular_chromosome_name_layer(
    layout, size, gap, family, colour, position))
  data <- linear_chromosome_name_data(layout, position)
  layer <- ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(x = .data$x, y = .data$y,
      label = .data$label, name_chr = .data$name_chr),
    stat = StatChromosomeName, geom = GeomChromosomeName,
    position = "identity",
    params = list(name_position = position, name_gap = gap,
      size = size, family = family, colour = colour, na.rm = FALSE),
    inherit.aes = FALSE
  )
  layer$ideogram_name_gap <- gap
  layer
}

linear_chromosome_name_data <- function(layout, position = "auto") {
  g <- layout$chrom
  if (position == "auto") position <- "start"
  vertical <- layout$orientation == "vertical"
  left <- layout$chromosome_width / 2 + layout$track_extent[["left"]]
  right <- layout$chromosome_width / 2 + layout$track_extent[["right"]]
  data <- data.frame(name_chr = g$.chr,
    x = (g$.axis_start_x + g$.axis_end_x) / 2,
    y = (g$.axis_start_y + g$.axis_end_y) / 2,
    label = g$.display_name %||% g$.chr, nx = 0, ny = 0)
  if (vertical) {
    data$x <- data$x + (right - left) / 2
    if (position == "middle") { data$x <- g$.axis_start_x + right; data$nx <- 1 }
    else {
      data$y <- if (position == "start") pmax(g$.axis_start_y, g$.axis_end_y) else
        pmin(g$.axis_start_y, g$.axis_end_y)
      data$ny <- if (position == "start") 1 else -1
    }
  } else {
    data$y <- data$y + (left - right) / 2
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
  extra_params = c("na.rm", "name_position", "name_gap"),
  compute_layer = function(self, data, params, layout) {
    anchors <- linear_chromosome_name_data(layout$coord$layout, params$name_position)
    i <- match(data$name_chr, anchors$name_chr)
    data[c("x", "y", "nx", "ny")] <- anchors[i, c("x", "y", "nx", "ny")]
    data$hjust <- data$vjust <- 0.5
    data
  })

GeomChromosomeName <- ggplot2::ggproto("GeomChromosomeName", ggplot2::GeomText,
  draw_panel = function(data, panel_params, coord, name_gap = 0.5,
      name_axis_clearance = NULL, na.rm = FALSE) {
    text <- ggplot2::GeomText$draw_panel(data, panel_params, coord, na.rm = na.rm)
    xy <- coord$transform(data, panel_params)
    dimensions <- chromosome_text_dimensions(data)
    offset <- abs(xy$nx) * dimensions$width / 2 + abs(xy$ny) * dimensions$height / 2 +
      grid::unit(name_gap * data$size, "mm")
    if (!is.null(name_axis_clearance)) offset <- offset + name_axis_clearance
    text$x <- text$x + xy$nx * offset
    text$y <- text$y + xy$ny * offset
    text
  })

linear_name_axis_clearance <- function(plot, name_layer) {
  layout <- ideogram_plot_layout(plot)
  names <- name_layer$data
  clearance <- grid::unit(rep(0, nrow(names)), "mm")
  previous <- name_layer$geom_params$name_axis_clearance
  hidden <- !names$name_chr %in% layout$chrom$.chr
  if (any(hidden) && !is.null(previous) && length(previous) == nrow(names))
    clearance[hidden] <- previous[hidden]
  axes <- Filter(function(layer) inherits(layer$geom, "GeomIdeogramAxisText") &&
    (inherits(layer$stat, "StatTrackAxis") || !is.null(layer$ideogram_track_title)), plot$layers)
  for (layer in axes) {
    title <- layer$ideogram_track_title
    object <- if (is.null(title)) layer$stat_params$axis_spec else list(
      size = layer$aes_params$size, family = layer$aes_params$family %||% "",
      fontface = layer$aes_params$fontface %||% 1,
      tick_length = layer$geom_params$tick_length, label_gap = layer$geom_params$label_gap)
    selected <- intersect(if (is.null(title)) layer$stat_params$chromosomes else title$chr,
      layout$chrom$.chr)
    if (!length(selected)) next
    for (i in seq_len(nrow(names))) {
      chr <- as.character(names$name_chr[i])
      if (!chr %in% selected) next
      data <- if (!is.null(title)) layer$data else {
        parts <- build_track_axis_layers(layout, object, chr)
        if (!length(parts)) next
        utils::tail(parts, 1)[[1]]$data
      }
      same_end <- data$nx * names$nx[i] + data$ny * names$ny[i] > 0.5
      if (!any(same_end)) next
      dimensions <- circular_text_dimensions(data$label[same_end], object$size,
        object$family, object$fontface %||% 1)
      extent <- abs(names$nx[i]) * dimensions$width + abs(names$ny[i]) * dimensions$height
      clearance[i] <- max(clearance[i], extent +
        grid::unit(object$tick_length + object$label_gap * object$size, "mm"))
    }
  }
  clearance
}

chromosome_text_dimensions <- function(data, parse = FALSE) {
  grobs <- lapply(seq_len(nrow(data)), function(i) {
    label <- if (parse) parse(text = as.character(data$label[i])) else data$label[i]
    grid::textGrob(label, gp = grid::gpar(fontsize = data$size[i] * ggplot2::.pt,
      fontfamily = data$family[i], fontface = data$fontface[i],
      lineheight = data$lineheight[i]))
  })
  list(width = do.call(grid::unit.c, lapply(grobs, grid::grobWidth)),
    height = do.call(grid::unit.c, lapply(grobs, grid::grobHeight)))
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
    if (is_circular_layout(layout) && is.null(breaks)) {
      local_breaks <- local_breaks[local_breaks < g$.end]
    }
    if (!length(local_breaks)) next
    projected <- project_positions_raw(layout$chrom, g$.chr, local_breaks)
    marker_extent <- (layout$marker_extent %||% c(left = 0, right = 0))[[side]]
    offset <- layout$chromosome_width / 2 +
      max(layout$track_extent[[side]], marker_extent) + gap
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
  spec <- list(chr = chr, side = side, breaks = breaks, n = n, units = units,
    labels = labels, gap = gap, tick_length = tick_length, label_gap = label_gap,
    size = size, family = family, colour = colour, linewidth = linewidth)
  for (part in seq_along(layers)) {
    layers[[part]]$ideogram_axis_spec <- spec
    layers[[part]]$ideogram_axis_part <- part
  }
  layers
}

# Include the axis spine in data bounds and reserve physical label overhang in
# the standard plot margin. No device dimensions enter the genomic layout.
reserve_chromosome_axis_space <- function(plot, layers) {
  layout <- ideogram_plot_layout(plot)
  for (layer in layers) {
    data <- layer$data
    layout <- include_ideogram_coordinates(layout, data$x, data$y, data$xend, data$yend)
    if (is_circular_layout(layout) && identical(layer$ideogram_axis_spec$side, "right")) {
      layout$circular_name_axis <- layer$ideogram_axis_spec
      if ("label" %in% names(data)) {
        layout$circular_name_axis$label <- data$label
        layout$circular_name_axis$label_x <- data$x
      }
    }
  }
  plot <- update_plot_ideogram_layout(plot, layout)
  needed <- grid::unit(rep(0, 4), "mm") # top, right, bottom, left
  reserve <- function(current, extent, at) {
    if (any(at)) max(current, extent[at]) else current
  }
  for (layer in layers) {
    data <- layer$data
    if (!all(c("nx", "ny") %in% names(data))) next
    if (is_circular_layout(layout)) data <- circular_transform_data(layout, data)
    tick <- layer$geom_params$tick_length %||% 0
    width <- height <- grid::unit(rep(0, nrow(data)), "mm")
    clearance <- grid::unit(tick, "mm")
    if ("label" %in% names(data)) {
      size <- layer$aes_params$size
      family <- layer$aes_params$family %||% ""
      fontface <- layer$aes_params$fontface %||% 1
      grobs <- lapply(as.character(data$label), function(label) {
        grid::textGrob(label, gp = grid::gpar(
          fontsize = size * ggplot2::.pt, fontfamily = family, fontface = fontface))
      })
      width <- do.call(grid::unit.c, lapply(grobs, grid::grobWidth))
      height <- do.call(grid::unit.c, lapply(grobs, grid::grobHeight))
      clearance <- grid::unit(tick + layer$geom_params$label_gap * size, "mm")
    }
    left <- data$x <= layout$bounds$x[1]
    right <- data$x >= layout$bounds$x[2]
    bottom <- data$y <= layout$bounds$y[1]
    top <- data$y >= layout$bounds$y[2]
    needed[4] <- reserve(needed[4], width + clearance, left & data$nx < 0)
    needed[4] <- reserve(needed[4], width / 2, left & data$nx == 0)
    needed[2] <- reserve(needed[2], width + clearance, right & data$nx > 0)
    needed[2] <- reserve(needed[2], width / 2, right & data$nx == 0)
    needed[3] <- reserve(needed[3], height + clearance, bottom & data$ny < 0)
    needed[3] <- reserve(needed[3], height / 2, bottom & data$ny == 0)
    needed[1] <- reserve(needed[1], height + clearance, top & data$ny > 0)
    needed[1] <- reserve(needed[1], height / 2, top & data$ny == 0)
  }
  if (!is_circular_layout(layout)) {
    names <- Filter(function(layer) !is.null(layer$ideogram_name_gap), plot$layers)
    for (layer in names) {
      axis_clearance <- linear_name_axis_clearance(plot, layer)
      layer$geom_params$name_axis_clearance <- axis_clearance
      dimensions <- circular_text_dimensions(layer$data$label, layer$aes_params$size,
        layer$aes_params$family %||% "")
      clearance <- grid::unit(layer$ideogram_name_gap * layer$aes_params$size, "mm") +
        max(axis_clearance)
      normals <- layer$data[c("nx", "ny")]
      if (any(normals$ny > 0)) needed[1] <- max(needed[1], dimensions$height + clearance)
      if (any(normals$nx > 0)) needed[2] <- max(needed[2], dimensions$width + clearance)
      if (any(normals$ny < 0)) needed[3] <- max(needed[3], dimensions$height + clearance)
      if (any(normals$nx < 0)) needed[4] <- max(needed[4], dimensions$width + clearance)
    }
  }
  margin <- plot$theme$plot.margin %||% theme_ideogram()$plot.margin
  # Keep text measurements lazy so the eventual device supplies font metrics.
  plot <- plot + ggplot2::theme(plot.margin = grid::unit.pmax(margin, needed))
  reserve_circular_text_space(reserve_ideogram_heading_space(plot, needed[1]))
}

reserve_ideogram_heading_space <- function(plot, needed) {
  clearance <- grid::convertHeight(needed, "mm", valueOnly = TRUE)
  heading_margin <- function(element) {
    margin <- element$margin %||% ggplot2::margin(0, 0, 0, 0)
    values <- grid::convertUnit(margin, "mm", valueOnly = TRUE)
    values[3] <- max(values[3], clearance)
    ggplot2::margin(values[1], values[2], values[3], values[4], unit = "mm")
  }
  plot + ggplot2::theme(
    plot.title = ggplot2::element_text(margin = heading_margin(plot$theme$plot.title)),
    plot.subtitle = ggplot2::element_text(margin = heading_margin(plot$theme$plot.subtitle)))
}
