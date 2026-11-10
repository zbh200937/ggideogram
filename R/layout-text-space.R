# Physical text extents and automatic spacing for chromosome components.

track_title_axis_clearance <- function(plot, title) {
  layout <- ideogram_plot_layout(plot)
  owner <- title$ideogram_track_title
  clearance <- grid::unit(0, "mm")
  gap <- 0
  for (layer in plot$layers) {
    if (!inherits(layer$stat, "StatTrackAxis") || layer$stat_params$axis_part != 3L) next
    object <- layer$stat_params$axis_spec
    if (!identical(object$track, owner$track) ||
        !owner$chr %in% layer$stat_params$chromosomes) next
    parts <- build_track_axis_layers(layout, object, owner$chr)
    if (!length(parts)) next
    ticks <- parts[[2]]$data
    if (!any(ticks$nx * title$data$nx + ticks$ny * title$data$ny > 0.5)) next
    gap <- max(gap, track_axis_gap(object$gap))
    data <- parts[[3]]$data
    extent <- grid::unit(object$tick_length, "mm")
    if (nrow(data)) {
      dimensions <- circular_text_dimensions(data$label, object$size,
        object$family %||% layout$base_spec$base_family %||% "")
      extent <- extent + grid::unit(object$label_gap * object$size, "mm") +
        abs(title$data$nx) * dimensions$width + abs(title$data$ny) * dimensions$height
    }
    clearance <- max(clearance, extent)
  }
  list(clearance = clearance, gap = gap)
}

align_linear_track_titles <- function(plot, layout) {
  indices <- which(vapply(plot$layers, function(layer)
    !is.null(layer$ideogram_track_title), logical(1)))
  if (!length(indices)) return(plot)
  titles <- plot$layers[indices]
  owners <- vapply(titles, function(layer) layer$ideogram_track_title$chr, character(1))
  if (any(!owners %in% layout$chrom$.chr)) return(plot)
  anchors <- lapply(titles, function(layer)
    layer$ideogram_track_title$anchor %||% layer$data[c("x", "y")])
  axis <- if (layout$orientation == "horizontal") "x" else "y"
  starts <- vapply(anchors, function(anchor) anchor[[axis]], numeric(1))
  for (group in split(seq_along(titles), match(starts, unique(starts)))) {
    spaces <- lapply(titles[group], function(layer) track_title_axis_clearance(plot, layer))
    gap <- max(vapply(spaces, function(space) space$gap, numeric(1)))
    label_offsets <- grid::unit(vapply(titles[group], function(layer) {
      layer$geom_params$tick_length + layer$geom_params$label_gap * layer$aes_params$size
    }, numeric(1)), "mm")
    axis_clearances <- do.call(grid::unit.c, lapply(spaces, function(space) space$clearance))
    clearance <- max(axis_clearances + label_offsets)
    for (j in seq_along(group)) {
      i <- group[j]
      layer <- copy_ideogram_layer(titles[[i]])
      anchor <- anchors[[i]]
      layer$ideogram_track_title$anchor <- anchor
      layer$data$x <- anchor$x + layer$data$nx * gap
      layer$data$y <- anchor$y + layer$data$ny * gap
      layer$geom_params$text_clearance <- clearance - label_offsets[j]
      plot$layers[[indices[i]]] <- layer
      layout <- include_ideogram_coordinates(layout, layer$data$x, layer$data$y)
    }
  }
  update_plot_ideogram_layout(plot, layout)
}

linear_name_axis_clearance <- function(plot, name_layer) {
  layout <- ideogram_plot_layout(plot)
  names <- name_layer$data
  clearance <- grid::unit(rep(0, nrow(names)), "mm")
  end_gap <- numeric(nrow(names))
  previous <- name_layer$geom_params$name_axis_clearance
  hidden <- !names$name_chr %in% layout$chrom$.chr
  if (any(hidden) && !is.null(previous) && length(previous) == nrow(names))
    clearance[hidden] <- previous[hidden]
  previous_end <- name_layer$geom_params$name_end_gap
  if (any(hidden) && length(previous_end) == nrow(names)) end_gap[hidden] <- previous_end[hidden]
  gene_labels <- Filter(function(layer) isTRUE(layer$ideogram_gene_label), plot$layers)
  for (layer in gene_labels) {
    for (i in seq_len(nrow(names))) {
      same_chr <- layer$data$.model_chr == names$name_chr[i]
      if (!any(same_chr)) next
      dimensions <- circular_text_dimensions(layer$data$.model_label[same_chr],
        layer$aes_params$size, layer$aes_params$family %||% "",
        layer$aes_params$fontface %||% 1)
      if (names$nx[i] < 0) {
        clearance[i] <- max(clearance[i], dimensions$width * layer$aes_params$hjust)
      } else if (names$ny[i] > 0) {
        clearance[i] <- max(clearance[i], dimensions$height * (1 - layer$aes_params$vjust))
      }
    }
  }
  axes <- Filter(function(layer) inherits(layer$geom, "GeomIdeogramAxisText") &&
    inherits(layer$stat, "StatTrackAxis"), plot$layers)
  for (layer in axes) {
    object <- layer$stat_params$axis_spec
    selected <- names$name_chr %in% layer$stat_params$chromosomes & !hidden
    for (i in which(selected)) {
      parts <- build_track_axis_layers(layout, object, names$name_chr[i])
      if (!length(parts)) next
      data <- parts[[3]]$data
      same_end <- data$nx * names$nx[i] + data$ny * names$ny[i] > 0.5
      if (!any(same_end)) next
      data <- data[same_end, , drop = FALSE]
      end_gap[i] <- max(end_gap[i], (data$x - names$x[i]) * names$nx[i] +
        (data$y - names$y[i]) * names$ny[i])
      dimensions <- circular_text_dimensions(data$label, object$size,
        object$family, object$fontface %||% 1)
      extent <- abs(names$nx[i]) * dimensions$width + abs(names$ny[i]) * dimensions$height
      clearance[i] <- max(clearance[i], extent +
        grid::unit(object$tick_length + object$label_gap * object$size, "mm"))
    }
  }

  alignment <- stats::setNames(vector("list", nrow(names)), names$name_chr)
  previous_alignment <- name_layer$geom_params$name_title_alignment
  if (any(hidden) && length(previous_alignment)) {
    alignment[hidden] <- previous_alignment[names$name_chr[hidden]]
  }
  titles <- Filter(function(layer) !is.null(layer$ideogram_track_title), plot$layers)
  if (length(titles)) {
    vertical <- layout$orientation == "vertical"
    along <- if (vertical) "y" else "x"
    anchors <- do.call(rbind, lapply(titles, function(layer) layer$data[c("x", "y")]))
    starts <- vapply(titles, function(layer)
      layer$ideogram_track_title$anchor[[along]], numeric(1))
    title_near <- do.call(grid::unit.c, lapply(titles, function(layer) {
      params <- layer$geom_params
      grid::unit(params$tick_length + params$label_gap * layer$aes_params$size, "mm") +
        (params$text_clearance %||% grid::unit(0, "mm"))
    }))
    dimensions <- lapply(titles, function(layer) {
      aes <- layer$aes_params
      circular_text_dimensions(layer$data$label, aes$size, aes$family, aes$fontface)
    })
    title_width <- do.call(grid::unit.c, lapply(dimensions, "[[", "width"))
    title_height <- do.call(grid::unit.c, lapply(dimensions, "[[", "height"))
    for (i in which(!hidden)) {
      outward <- if (vertical) names$ny[i] else -names$nx[i]
      group <- which(abs(starts - names[[along]][i]) < 1e-8)
      if (outward <= 0 || !length(group)) next
      alignment[[i]] <- list(anchors = anchors[group, , drop = FALSE],
        near = title_near[group], half = (if (vertical) title_width else title_height)[group] / 2,
        min_gap = end_gap[i], min_clearance = clearance[i])
      extent <- if (vertical) title_height else title_width
      end_gap[i] <- max(end_gap[i], abs(anchors[[along]][group] - starts[group]))
      clearance[i] <- max(clearance[i], title_near[group] + extent[group])
    }
  }
  list(clearance = clearance, end_gap = end_gap, alignment = alignment)
}

linear_name_alignment_grob <- function(text, data, xy, dimensions, coord, panel, gap, alignment) {
  alignment <- alignment[as.character(data$name_chr)]
  if (!any(lengths(alignment) > 0)) return(text)
  vertical <- coord$layout$orientation == "vertical"
  along <- if (vertical) "y" else "x"
  across <- if (vertical) "x" else "y"
  sign <- if (vertical) 1 else -1
  name_along <- if (vertical) dimensions$height else dimensions$width
  name_across <- if (vertical) dimensions$width else dimensions$height
  candidates <- lapply(seq_len(nrow(data)), function(i) {
    group <- alignment[[i]]
    if (is.null(group)) return(NULL)
    titles <- coord$transform(group$anchors, panel)
    end <- data[i, c("x", "y"), drop = FALSE]
    end[[along]] <- end[[along]] + sign * group$min_gap
    end <- coord$transform(end, panel)
    list(near = grid::unit(titles[[along]], "native") + sign * group$near,
      minimum = grid::unit(end[[along]], "native") + sign * (group$min_clearance + gap[i]),
      delta = grid::unit(titles[[across]] - xy[[across]][i], "native"),
      required = group$half + name_across[i] / 2,
      half = name_along[i] / 2)
  })
  grid::gTree(text = text, candidates = candidates, vertical = vertical,
    cl = "ideogram_name_alignment")
}

#' @method makeContent ideogram_name_alignment
#' @export
makeContent.ideogram_name_alignment <- function(x) {
  text <- x$text
  along <- if (x$vertical) "y" else "x"
  sign <- if (x$vertical) 1 else -1
  along_mm <- if (x$vertical) grid::convertHeight else grid::convertWidth
  across_mm <- if (x$vertical) grid::convertWidth else grid::convertHeight
  for (i in seq_along(x$candidates)) {
    z <- x$candidates[[i]]
    if (is.null(z)) next
    separated <- abs(across_mm(z$delta, "mm", TRUE)) >= across_mm(z$required, "mm", TRUE)
    if (!all(separated)) next
    edge <- if (x$vertical) max(z$near) else min(z$near)
    if (sign * along_mm(edge - z$minimum, "mm", TRUE) < -1e-8) next
    text[[along]][i] <- edge + sign * z$half
  }
  grid::setChildren(x, grid::gList(text))
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

rotated_text_dimensions <- function(dimensions, angle) {
  radians <- angle * pi / 180
  cosine <- abs(cos(radians))
  sine <- abs(sin(radians))
  list(width = cosine * dimensions$width + sine * dimensions$height,
    height = sine * dimensions$width + cosine * dimensions$height)
}

axis_text_clearance_grob <- function(text, data, coord, panel, drop_labels) {
  layout <- coord$layout
  chrom <- layout$chrom
  vertical <- layout$orientation == "vertical"
  direction <- if (vertical) data$nx else -data$ny
  chr <- match(data$axis_chr, chrom$.chr)
  neighbor <- match(paste(chrom$.row[chr], chrom$.col[chr] + sign(direction)),
    paste(chrom$.row, chrom$.col))
  rows <- which(!is.na(neighbor))
  if (!length(rows)) return(text)
  other <- chrom[neighbor[rows], , drop = FALSE]
  right <- direction[rows] > 0
  side <- ifelse(right, "left", "right")
  space <- annotation_extent(layout, other$.chr, side)
  boundary <- data[rows, c("x", "y"), drop = FALSE]
  if (vertical) {
    boundary$x <- ifelse(right, other$.x_min - space, other$.x_max + space)
    boundary$y <- other$.y_min
    boundary$yend <- other$.y_max
  } else {
    boundary$y <- ifelse(right, other$.y_max + space, other$.y_min - space)
    boundary$x <- other$.x_min
    boundary$xend <- other$.x_max
  }
  boundary <- coord$transform(boundary, panel)
  grid::gTree(text = text, data = data[rows, , drop = FALSE], rows = rows,
    direction = if (vertical) direction[rows] else -direction[rows],
    limit = if (vertical) boundary$x else boundary$y,
    from = if (vertical) boundary$y else boundary$x,
    to = if (vertical) boundary$yend else boundary$xend,
    vertical = vertical, drop_labels = drop_labels, cl = "ideogram_axis_clearance")
}

#' @method makeContent ideogram_axis_clearance
#' @export
makeContent.ideogram_axis_clearance <- function(x) {
  text <- x$text
  dimensions <- rotated_text_dimensions(chromosome_text_dimensions(x$data), text$rot[x$rows])
  limit <- grid::unit(x$limit, "native")
  if (x$vertical) {
    room <- grid::convertWidth(limit - text$x[x$rows], "mm", TRUE)
    extent <- grid::convertWidth(dimensions$width, "mm", TRUE) *
      ifelse(text$hjust[x$rows] == 0.5, 0.5, 1)
    along <- grid::convertY(text$y[x$rows], "mm", TRUE)
    half <- grid::convertHeight(dimensions$height, "mm", TRUE) / 2
    convert <- grid::convertY
  } else {
    room <- grid::convertHeight(limit - text$y[x$rows], "mm", TRUE)
    extent <- grid::convertHeight(dimensions$height, "mm", TRUE) *
      ifelse(text$vjust[x$rows] == 0.5, 0.5, 1)
    along <- grid::convertX(text$x[x$rows], "mm", TRUE)
    half <- grid::convertWidth(dimensions$width, "mm", TRUE) / 2
    convert <- grid::convertX
  }
  from <- convert(grid::unit(x$from, "native"), "mm", TRUE)
  to <- convert(grid::unit(x$to, "native"), "mm", TRUE)
  crowded <- along + half > from & along - half < to & extent > x$direction * room + 1e-8
  if (any(crowded)) {
    if (x$drop_labels) {
      text$label[x$rows[crowded]] <- ""
    } else {
      warning(sprintf(paste0("BP axis labels for %s extend into the space reserved for an adjacent chromosome. ",
        "Adjust the axis label angle, increase `chromosome_gap` or use a larger output size."),
        format_chr_rows(unique(x$data$axis_chr[crowded]))), call. = FALSE)
    }
  }
  grid::setChildren(x, grid::gList(text))
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
  if (!is_circular_layout(layout)) {
    plot <- align_linear_track_titles(plot, layout)
    layout <- ideogram_plot_layout(plot)
  }
  layers <- Filter(function(layer) !is.null(layer$ideogram_axis_spec) ||
    inherits(layer$stat, "StatTrackAxis") || !is.null(layer$ideogram_track_title), plot$layers)
  needed <- grid::unit(rep(0, 4), "mm") # top, right, bottom, left
  reserve <- function(current, extent, at) {
    if (any(at)) max(current, extent[at]) else current
  }
  for (layer in layers) {
    data <- layer$data
    if (!nrow(data) || !all(c("nx", "ny") %in% names(data))) next
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
      dimensions <- rotated_text_dimensions(list(width = width, height = height),
        layer$aes_params$angle %||% 0)
      width <- dimensions$width; height <- dimensions$height
      clearance <- grid::unit(tick + layer$geom_params$label_gap * size, "mm") +
        (layer$geom_params$text_clearance %||% grid::unit(0, "mm"))
    }
    left <- data$x <= layout$bounds$x[1]
    right <- data$x >= layout$bounds$x[2]
    bottom <- data$y <= layout$bounds$y[1]
    top <- data$y >= layout$bounds$y[2]
    needed[4] <- reserve(needed[4], width + clearance, left & data$nx < 0)
    needed[4] <- reserve(needed[4], width / 2, data$nx == 0)
    needed[2] <- reserve(needed[2], width + clearance, right & data$nx > 0)
    needed[2] <- reserve(needed[2], width / 2, data$nx == 0)
    needed[3] <- reserve(needed[3], height + clearance, bottom & data$ny < 0)
    needed[3] <- reserve(needed[3], height / 2, data$ny == 0)
    needed[1] <- reserve(needed[1], height + clearance, top & data$ny > 0)
    needed[1] <- reserve(needed[1], height / 2, data$ny == 0)
  }
  if (!is_circular_layout(layout)) {
    gene_labels <- Filter(function(layer) isTRUE(layer$ideogram_gene_label), plot$layers)
    for (layer in gene_labels) {
      dimensions <- circular_text_dimensions(layer$data$.model_label,
        layer$aes_params$size, layer$aes_params$family %||% "",
        layer$aes_params$fontface %||% 1)
      if (layout$orientation == "horizontal") {
        needed[4] <- max(needed[4], dimensions$width * layer$aes_params$hjust)
      } else {
        needed[1] <- max(needed[1], dimensions$height * (1 - layer$aes_params$vjust))
      }
    }
    names <- which(vapply(plot$layers, function(layer)
      !is.null(layer$ideogram_name_gap), logical(1)))
    for (i in names) {
      layer <- copy_ideogram_layer(plot$layers[[i]])
      name_space <- linear_name_axis_clearance(plot, layer)
      layer$geom_params$name_axis_clearance <- name_space$clearance
      layer$geom_params$name_end_gap <- name_space$end_gap
      layer$geom_params$name_title_alignment <- name_space$alignment
      plot$layers[[i]] <- layer
      dimensions <- circular_text_dimensions(layer$data$label, layer$aes_params$size,
        layer$aes_params$family %||% "")
      clearance <- grid::unit(layer$ideogram_name_gap * layer$aes_params$size, "mm") +
        max(name_space$clearance)
      normals <- layer$data[c("nx", "ny")]
      if (any(normals$nx == 0)) {
        needed[c(2, 4)] <- grid::unit.pmax(needed[c(2, 4)], max(dimensions$width) / 2)
      }
      if (any(normals$ny == 0)) {
        needed[c(1, 3)] <- grid::unit.pmax(needed[c(1, 3)], max(dimensions$height) / 2)
      }
      if (any(normals$ny > 0)) needed[1] <- max(needed[1], dimensions$height + clearance)
      if (any(normals$nx > 0)) needed[2] <- max(needed[2], dimensions$width + clearance)
      if (any(normals$ny < 0)) needed[3] <- max(needed[3], dimensions$height + clearance)
      if (any(normals$nx < 0)) needed[4] <- max(needed[4], dimensions$width + clearance)
    }
  }
  needed <- grid::unit.pmax(needed, circular_text_space(plot))
  margin <- plot$theme$plot.margin %||% theme_ideogram()$plot.margin
  # Keep text measurements lazy so the eventual device supplies font metrics.
  plot <- plot + ggplot2::theme(plot.margin = auto_text_space(margin, needed))
  spacing <- plot$theme$legend.box.spacing
  if (isTRUE(attr(spacing, "ideogram_auto_spacing"))) {
    spacing <- auto_text_space(spacing, needed[3])
    attr(spacing, "ideogram_auto_spacing") <- TRUE
    plot <- plot + ggplot2::theme(legend.box.spacing = spacing)
  }
  reserve_ideogram_heading_space(plot, needed[1])
}

auto_text_space <- function(current, needed) {
  # Keep the user's base separately from the allowance for current contents.
  base <- attr(current, "ideogram_space_base") %||% current
  result <- grid::unit.pmax(base, needed)
  attr(result, "ideogram_space_base") <- base
  result
}

reserve_ideogram_heading_space <- function(plot, needed) {
  clearance <- grid::convertHeight(needed, "mm", valueOnly = TRUE)
  heading_margin <- function(element) {
    margin <- element$margin %||% ggplot2::margin(0, 0, 0, 0)
    margin <- attr(margin, "ideogram_space_base") %||% margin
    values <- grid::convertUnit(margin, "mm", valueOnly = TRUE)
    values[3] <- max(values[3], clearance)
    result <- ggplot2::margin(values[1], values[2], values[3], values[4], unit = "mm")
    attr(result, "ideogram_space_base") <- margin
    result
  }
  plot + ggplot2::theme(
    plot.title = ggplot2::element_text(margin = heading_margin(plot$theme$plot.title)),
    plot.subtitle = ggplot2::element_text(margin = heading_margin(plot$theme$plot.subtitle)))
}

circular_text_space <- function(plot) {
  layout <- ideogram_plot_layout(plot)
  needed <- grid::unit(rep(0, 4), "mm")
  if (!is_circular_layout(layout)) return(needed)
  axis_offset <- grid::unit(0, "mm")
  axis <- layout$circular_name_axis
  if (!is.null(axis)) {
    axis_offset <- grid::unit(axis$tick_length, "mm")
    if (isTRUE(axis$labels) && length(axis$label)) {
      dimensions <- rotated_text_dimensions(circular_text_dimensions(axis$label, axis$size, axis$family),
        axis$angle %||% 0)
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
  needed
}
