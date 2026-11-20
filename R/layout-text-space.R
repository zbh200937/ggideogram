# Physical text extents and automatic spacing for chromosome components.

track_title_axis_clearance <- function(plot, title) {
  layout <- ideogram_plot_layout(plot)
  owner <- title$ideogram_track_title
  clearance <- grid::unit(0, "mm")
  gap <- 0
  has_start_labels <- FALSE
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
      has_start_labels <- TRUE
      size <- parts[[3]]$aes_params$size
      dimensions <- circular_text_dimensions(data$label, size,
        object$family %||% layout$base_spec$base_family %||% "")
      extent <- extent + grid::unit(object$label_gap * size, "mm") +
        (if (abs(title$data$nx) > 0) dimensions$width else dimensions$height)
    }
    clearance <- max(clearance, extent)
  }
  for (layer in plot$layers) {
    if (!isTRUE(layer$ideogram_gene_label_start) ||
        !identical(layer$ideogram_gene_track, owner$track)) next
    same_chr <- layer$data$.model_chr == owner$chr
    if (!any(same_chr)) next
    outward <- if (layout$orientation == "vertical") title$data$ny else -title$data$nx
    if (outward <= 0) next
    dimensions <- circular_text_dimensions(layer$data$.model_label[same_chr],
      layer$aes_params$size, layer$aes_params$family %||% "",
      layer$aes_params$fontface %||% 1)
    gene_ext <- if (abs(title$data$nx) > 0) {
      max(dimensions$width) * layer$aes_params$hjust
    } else {
      max(dimensions$height) * (1 - layer$aes_params$vjust)
    }
    clearance <- max(clearance, gene_ext)
  }
  list(clearance = clearance, gap = gap, has_start_labels = has_start_labels)
}

align_linear_track_titles <- function(plot, layout) {
  indices <- which(vapply(plot$layers, function(layer)
    !is.null(layer$ideogram_track_title), logical(1)))
  if (!length(indices)) return(plot)
  titles <- plot$layers[indices]
  owners <- vapply(titles, function(layer) layer$ideogram_track_title$chr, character(1))
  if (any(!owners %in% layout$chrom$.chr)) return(plot)
  for (i in seq_along(titles)) {
    layer <- copy_ideogram_layer(titles[[i]])
    meta <- layer$ideogram_track_title
    anchor <- meta$anchor %||% layer$data[c("x", "y")]
    meta$anchor <- anchor
    if (identical(layout$orientation, "horizontal") && identical(meta$position, "top")) {
      if (isTRUE(meta$auto_angle)) layer$aes_params$angle <- 0
      layer$data$x <- anchor$x
      layer$data$y <- anchor$y
      layer$geom_params$text_clearance <- grid::unit(0, "mm")
    } else {
      space <- track_title_axis_clearance(plot, layer)
      if (isTRUE(meta$auto_angle)) {
        spec_row <- track_table_row(layout, meta$track)
        disc_levels <- if (!is.null(spec_row)) {
          rng_key <- track_key_for_rows(spec_row, meta$chr)[1]
          layout$track_ranges[[rng_key]]$discrete_levels
        } else NULL
        use_vert <- identical(layout$orientation, "horizontal") &&
          identical(meta$position %||% "side", "side") &&
          isTRUE(space$has_start_labels) &&
          length(disc_levels) >= 3L &&
          !is.null(spec_row) && isTRUE(spec_row$width >= 5)
        layer$aes_params$angle <- if (use_vert) 90 else 0
      }
      layer$data$x <- anchor$x + layer$data$nx * space$gap
      layer$data$y <- anchor$y + layer$data$ny * space$gap
      layer$geom_params$text_clearance <- space$clearance
    }
    layer$ideogram_track_title <- meta
    plot$layers[[indices[i]]] <- layer
    layout <- include_ideogram_coordinates(layout, layer$data$x, layer$data$y)
  }
  update_plot_ideogram_layout(plot, layout)
}

per_label_text_dimensions <- function(label, size, family = "", fontface = 1) {
  lbl_chr <- as.character(label)
  if (!length(lbl_chr)) {
    z <- grid::unit(numeric(0), "mm")
    return(list(width = z, height = z))
  }
  ulab <- unique(lbl_chr)
  gp <- grid::gpar(fontsize = size * ggplot2::.pt,
    fontfamily = family %||% "", fontface = fontface %||% 1)
  grobs <- lapply(ulab, function(value) grid::textGrob(value, gp = gp))
  u_w <- do.call(grid::unit.c, lapply(grobs, grid::grobWidth))
  u_h <- do.call(grid::unit.c, lapply(grobs, grid::grobHeight))
  m_idx <- match(lbl_chr, ulab)
  list(width = u_w[m_idx], height = u_h[m_idx])
}

linear_name_axis_clearance <- function(plot, name_layer) {
  layout <- ideogram_plot_layout(plot)
  names <- name_layer$data
  vertical <- layout$orientation == "vertical"
  along <- if (vertical) "y" else "x"
  clearance_init <- grid::unit(rep(0, nrow(names)), "mm")
  margin_init <- grid::unit(rep(0, nrow(names)), "mm")
  bp_list <- replicate(nrow(names), grid::unit(0, "mm"), simplify = FALSE)
  end_gap <- numeric(nrow(names))
  bp_end_gap <- numeric(nrow(names))
  previous <- name_layer$geom_params$name_axis_clearance
  previous_margin <- name_layer$geom_params$name_margin_clearance
  hidden <- !names$name_chr %in% layout$chrom$.chr
  clearance_list <- lapply(seq_len(nrow(names)), function(i) {
    if (hidden[i] && !is.null(previous) && length(previous) == nrow(names)) previous[i] else clearance_init[i]
  })
  margin_list <- lapply(seq_len(nrow(names)), function(i) {
    if (hidden[i] && !is.null(previous_margin) && length(previous_margin) == nrow(names)) previous_margin[i] else margin_init[i]
  })
  previous_end <- name_layer$geom_params$name_end_gap
  if (any(hidden) && length(previous_end) == nrow(names)) end_gap[hidden] <- previous_end[hidden]

  obs_list <- vector("list", nrow(names))
  add_obs <- function(idx, anch, reach_u, half_u) {
    n_a <- nrow(anch)
    if (!n_a) return()
    if (length(reach_u) == 1L && n_a > 1L) reach_u <- rep(reach_u, n_a)
    if (length(half_u) == 1L && n_a > 1L) half_u <- rep(half_u, n_a)
    if (is.null(obs_list[[idx]])) {
      obs_list[[idx]] <<- list(anchors = anch[c("x", "y")], reach = reach_u, half = half_u)
    } else {
      obs_list[[idx]]$anchors <<- rbind(obs_list[[idx]]$anchors, anch[c("x", "y")])
      obs_list[[idx]]$reach <<- grid::unit.c(obs_list[[idx]]$reach, reach_u)
      obs_list[[idx]]$half <<- grid::unit.c(obs_list[[idx]]$half, half_u)
    }
  }

  bp_axes <- Filter(function(layer) !is.null(layer$ideogram_axis_spec) &&
    inherits(layer$geom, "GeomIdeogramAxisText"), plot$layers)
  for (layer in bp_axes) {
    d_ax <- layer$data
    if (!nrow(d_ax) || !"axis_chr" %in% names(d_ax)) next
    spec_ax <- layer$ideogram_axis_spec
    for (i in which(!hidden)) {
      hit_ax <- !is.na(d_ax$axis_chr) & d_ax$axis_chr == names$name_chr[i] &
        (d_ax$nx * names$nx[i] + d_ax$ny * names$ny[i] > 0.5)
      if (!any(hit_ax)) next
      sub_ax <- d_ax[hit_ax, , drop = FALSE]
      gap_val <- max((sub_ax$x - names$x[i]) * names$nx[i] +
        (sub_ax$y - names$y[i]) * names$ny[i])
      bp_end_gap[i] <- max(bp_end_gap[i], gap_val)
      end_gap[i] <- max(end_gap[i], gap_val)
      dims_ax <- rotated_text_dimensions(
        circular_text_dimensions(sub_ax$label, spec_ax$size, spec_ax$family %||% ""),
        spec_ax$angle %||% 0)
      ext_ax <- (if (abs(names$nx[i]) > 0) dims_ax$width else dims_ax$height) +
        grid::unit(spec_ax$tick_length + spec_ax$label_gap * spec_ax$size, "mm")
      bp_list[[i]] <- max(grid::unit.c(bp_list[[i]], ext_ax))
      clearance_list[[i]] <- max(grid::unit.c(clearance_list[[i]], ext_ax))
      margin_list[[i]] <- max(grid::unit.c(margin_list[[i]], ext_ax))
    }
  }

  gene_labels <- Filter(function(layer) isTRUE(layer$ideogram_gene_label_start), plot$layers)
  for (layer in gene_labels) {
    for (i in which(!hidden)) {
      if (!(names$nx[i] < 0 || names$ny[i] > 0)) next
      same_chr <- layer$data$.model_chr == names$name_chr[i]
      if (!any(same_chr)) next
      dimensions <- per_label_text_dimensions(layer$data$.model_label[same_chr],
        layer$aes_params$size, layer$aes_params$family %||% "",
        layer$aes_params$fontface %||% 1)
      g_reach <- if (names$nx[i] < 0) {
        dimensions$width * layer$aes_params$hjust
      } else {
        dimensions$height * (1 - layer$aes_params$vjust)
      }
      g_half <- (if (names$nx[i] < 0) dimensions$height else dimensions$width) / 2
      g_anch <- data.frame(x = layer$data$.lx[same_chr], y = layer$data$.ly[same_chr])
      add_obs(i, g_anch, g_reach, g_half)
      clearance_list[[i]] <- max(grid::unit.c(clearance_list[[i]], g_reach))
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
      ticks <- parts[[2]]$data
      same_tick <- ticks$nx * names$nx[i] + ticks$ny * names$ny[i] > 0.5
      if (!any(same_tick)) next
      ticks <- ticks[same_tick, , drop = FALSE]
      end_gap[i] <- max(end_gap[i], (ticks$x - names$x[i]) * names$nx[i] +
        (ticks$y - names$y[i]) * names$ny[i])
      tick_u <- grid::unit(rep(object$tick_length, nrow(ticks)), "mm")
      zero_u <- grid::unit(rep(0, nrow(ticks)), "mm")
      add_obs(i, ticks, tick_u, zero_u)
      clearance_list[[i]] <- max(grid::unit.c(clearance_list[[i]], grid::unit(object$tick_length, "mm")))
      data <- parts[[3]]$data
      same_end <- data$nx * names$nx[i] + data$ny * names$ny[i] > 0.5
      if (any(same_end)) {
        data <- data[same_end, , drop = FALSE]
        size <- parts[[3]]$aes_params$size
        dimensions <- per_label_text_dimensions(data$label, size,
          object$family, object$fontface %||% 1)
        extent <- (if (abs(names$nx[i]) > 0) dimensions$width else dimensions$height) +
          grid::unit(object$tick_length + object$label_gap * size, "mm")
        lbl_half <- (if (vertical) dimensions$width else dimensions$height) / 2
        add_obs(i, data, extent, lbl_half)
        clearance_list[[i]] <- max(grid::unit.c(clearance_list[[i]], extent))
      }
    }
  }

  titles <- Filter(function(layer) !is.null(layer$ideogram_track_title), plot$layers)
  if (length(titles)) {
    side_titles <- Filter(function(layer) {
      !(identical(layout$orientation, "horizontal") &&
        identical(layer$ideogram_track_title$position, "top"))
    }, titles)
    if (length(side_titles)) {
      anchors <- do.call(rbind, lapply(side_titles, function(layer) layer$data[c("x", "y")]))
      starts <- vapply(side_titles, function(layer)
        layer$ideogram_track_title$anchor[[along]], numeric(1))
      title_near <- do.call(grid::unit.c, lapply(side_titles, function(layer) {
        params <- layer$geom_params
        grid::unit(params$tick_length + params$label_gap * layer$aes_params$size, "mm") +
          (params$text_clearance %||% grid::unit(0, "mm"))
      }))
      dimensions <- lapply(side_titles, function(layer) {
        aes <- layer$aes_params
        rotated_text_dimensions(
          per_label_text_dimensions(layer$data$label, aes$size, aes$family, aes$fontface),
          aes$angle %||% 0)
      })
      title_width <- do.call(grid::unit.c, lapply(dimensions, "[[", "width"))
      title_height <- do.call(grid::unit.c, lapply(dimensions, "[[", "height"))
      t_reach <- title_near + (if (vertical) title_height else title_width)
      t_half <- (if (vertical) title_width else title_height) / 2
      for (i in which(!hidden)) {
        outward <- if (vertical) names$ny[i] else -names$nx[i]
        group <- which(abs(starts - names[[along]][i]) < 1e-8)
        if (outward <= 0 || !length(group)) next
        add_obs(i, anchors[group, , drop = FALSE], t_reach[group], t_half[group])
        end_gap[i] <- max(end_gap[i], abs(anchors[[along]][group] - starts[group]))
        clearance_list[[i]] <- max(grid::unit.c(clearance_list[[i]], t_reach[group]))
      }
    }
  }

  alignment <- stats::setNames(vector("list", nrow(names)), names$name_chr)
  previous_alignment <- name_layer$geom_params$name_title_alignment
  if (any(hidden) && length(previous_alignment)) {
    alignment[hidden] <- previous_alignment[names$name_chr[hidden]]
  }
  for (i in which(!hidden)) {
    dir_sign <- if (vertical) names$ny[i] else names$nx[i]
    if (!is.null(obs_list[[i]]) && dir_sign != 0) {
      alignment[[i]] <- list(
        sign = dir_sign,
        min_gap = bp_end_gap[i],
        min_clearance = bp_list[[i]],
        obstacles = obs_list[[i]]
      )
      obs_ac <- if (vertical) obs_list[[i]]$anchors$x - names$x[i] else
        obs_list[[i]]$anchors$y - names$y[i]
      tol_ac <- if (vertical) layout$chromosome_width * 3 else layout$chromosome_width * 0.75
      near_body <- which(abs(obs_ac) <= tol_ac)
      if (length(near_body)) {
        margin_list[[i]] <- max(grid::unit.c(margin_list[[i]], obs_list[[i]]$reach[near_body]))
      }
    }
  }
  clearance <- do.call(grid::unit.c, clearance_list)
  margin_clearance <- do.call(grid::unit.c, margin_list)
  list(clearance = clearance, margin_clearance = margin_clearance,
    end_gap = end_gap, alignment = alignment)
}

linear_name_alignment_grob <- function(text, data, xy, dimensions, coord, panel, gap, alignment) {
  alignment <- alignment[as.character(data$name_chr)]
  if (!any(lengths(alignment) > 0)) return(text)
  vertical <- coord$layout$orientation == "vertical"
  along <- if (vertical) "y" else "x"
  across <- if (vertical) "x" else "y"
  name_along <- if (vertical) dimensions$height else dimensions$width
  name_across <- if (vertical) dimensions$width else dimensions$height
  candidates <- lapply(seq_len(nrow(data)), function(i) {
    group <- alignment[[i]]
    if (is.null(group)) return(NULL)
    s <- group$sign %||% (if (vertical) 1 else -1)
    end <- data[i, c("x", "y"), drop = FALSE]
    end[[along]] <- end[[along]] + s * (group$min_gap %||% 0)
    end <- coord$transform(end, panel)
    cand <- list(
      sign = s,
      half = name_along[i] / 2,
      minimum = grid::unit(end[[along]], "native") +
        s * ((group$min_clearance %||% grid::unit(0, "mm")) + gap[i])
    )
    if (!is.null(group$obstacles)) {
      obs_xy <- coord$transform(group$obstacles$anchors, panel)
      cand$obstacles <- list(
        delta = grid::unit(obs_xy[[across]] - xy[[across]][i], "native"),
        required = group$obstacles$half + name_across[i] / 2,
        edge = grid::unit(obs_xy[[along]], "native") + s * (group$obstacles$reach + gap[i])
      )
    }
    cand
  })
  grid::gTree(text = text, candidates = candidates, vertical = vertical,
    cl = "ideogram_name_alignment")
}

#' @method makeContent ideogram_name_alignment
#' @export
makeContent.ideogram_name_alignment <- function(x) {
  text <- x$text
  along <- if (x$vertical) "y" else "x"
  pos_mm <- if (x$vertical) grid::convertY else grid::convertX
  along_mm <- if (x$vertical) grid::convertHeight else grid::convertWidth
  across_mm <- if (x$vertical) grid::convertWidth else grid::convertHeight
  new_coords <- lapply(seq_len(length(text[[along]])), function(i) text[[along]][i])
  for (i in seq_along(x$candidates)) {
    z <- x$candidates[[i]]
    if (is.null(z)) next
    s <- z$sign
    min_mm <- pos_mm(z$minimum, "mm", TRUE)
    if (!is.null(z$obstacles)) {
      hit <- abs(across_mm(z$obstacles$delta, "mm", TRUE)) <
        across_mm(z$obstacles$required, "mm", TRUE)
      if (any(hit)) {
        edges_mm <- pos_mm(z$obstacles$edge[hit], "mm", TRUE)
        obs_mm <- if (s > 0) max(edges_mm) else min(edges_mm)
        if (s * (obs_mm - min_mm) > 0) {
          min_mm <- obs_mm
        }
      }
    }
    half_mm <- along_mm(z$half, "mm", TRUE)
    new_coords[[i]] <- grid::unit(min_mm + s * half_mm, "mm")
  }
  text[[along]] <- do.call(grid::unit.c, new_coords)
  grid::setChildren(x, grid::gList(text))
}

chromosome_text_dimensions <- function(data, parse = FALSE) {
  n <- nrow(data)
  if (!n) {
    z <- grid::unit(numeric(0), "mm")
    return(list(width = z, height = z))
  }
  lbl <- as.character(data$label)
  key <- paste(is.na(lbl), lbl, data$size, data$family, data$fontface, data$lineheight, sep = "\r")
  ukey <- unique(key)
  first_i <- match(ukey, key)
  grobs <- lapply(first_i, function(i) {
    label <- if (parse) parse(text = lbl[i]) else data$label[i]
    grid::textGrob(label, gp = grid::gpar(fontsize = data$size[i] * ggplot2::.pt,
      fontfamily = data$family[i], fontface = data$fontface[i],
      lineheight = data$lineheight[i]))
  })
  u_w <- do.call(grid::unit.c, lapply(grobs, grid::grobWidth))
  u_h <- do.call(grid::unit.c, lapply(grobs, grid::grobHeight))
  if (length(ukey) == n) return(list(width = u_w, height = u_h))
  m <- match(key, ukey)
  list(width = u_w[m], height = u_h[m])
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
    is_cart <- "ideogram_cartesian" %in% names(data) && any(data$ideogram_cartesian %in% TRUE)
    layout <- include_ideogram_coordinates(
      layout, data$x, data$y, data$xend, data$yend, cartesian = is_cart)
    if (is_circular_layout(layout) && identical(layer$ideogram_axis_spec$side, "right")) {
      spec <- layer$ideogram_axis_spec
      prev <- layout$circular_name_axis
      if (is.null(prev)) {
        layout$circular_name_axis <- spec
        if ("label" %in% names(data)) {
          layout$circular_name_axis$label <- as.character(data$label)
          layout$circular_name_axis$label_x <- data$x
        }
      } else {
        prev$gap <- max(prev$gap, spec$gap)
        prev$tick_length <- max(prev$tick_length, spec$tick_length)
        prev$label_gap <- max(prev$label_gap, spec$label_gap)
        prev$size <- max(prev$size, spec$size)
        prev$labels <- isTRUE(prev$labels) || isTRUE(spec$labels)
        if ("label" %in% names(data)) {
          prev$label <- c(prev$label %||% character(), as.character(data$label))
          prev$label_x <- c(prev$label_x %||% numeric(), data$x)
        }
        layout$circular_name_axis <- prev
      }
    }
  }
  plot <- update_plot_ideogram_layout(plot, layout)
  if (!is_circular_layout(layout)) {
    plot <- align_linear_track_titles(plot, layout)
    layout <- ideogram_plot_layout(plot)
  }
  layers <- Filter(function(layer) !is.null(layer$ideogram_axis_spec) ||
    inherits(layer$stat, "StatTrackAxis") || !is.null(layer$ideogram_track_title) ||
    isTRUE(layer$ideogram_serial_row_label), plot$layers)
  needed <- grid::unit(rep(0, 4), "mm") # top, right, bottom, left
  reserve <- function(current, extent, at) {
    if (any(at)) max(current, extent[at]) else current
  }
  tol <- if (is_circular_layout(layout)) {
    max(layout$chromosome_width, 0.15 * max(abs(unlist(layout$bounds))))
  } else layout$chromosome_width
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
      lbl_chr <- as.character(data$label)
      ulab <- unique(lbl_chr)
      gp <- grid::gpar(fontsize = size * ggplot2::.pt, fontfamily = family, fontface = fontface)
      grobs <- lapply(ulab, function(label) grid::textGrob(label, gp = gp))
      u_w <- do.call(grid::unit.c, lapply(grobs, grid::grobWidth))
      u_h <- do.call(grid::unit.c, lapply(grobs, grid::grobHeight))
      m_idx <- match(lbl_chr, ulab)
      width <- u_w[m_idx]
      height <- u_h[m_idx]
      dimensions <- rotated_text_dimensions(list(width = width, height = height),
        layer$aes_params$angle %||% 0)
      width <- dimensions$width; height <- dimensions$height
      clearance <- grid::unit(tick + layer$geom_params$label_gap * size, "mm") +
        (layer$geom_params$text_clearance %||% grid::unit(0, "mm"))
    }
    left <- data$x <= layout$bounds$x[1] + tol
    right <- data$x >= layout$bounds$x[2] - tol
    bottom <- data$y <= layout$bounds$y[1] + tol
    top <- data$y >= layout$bounds$y[2] - tol
    zero_mm <- grid::unit(0, "mm")
    span_x <- max(1e-6, diff(layout$bounds$x))
    span_y <- max(1e-6, diff(layout$bounds$y))
    d_left_mm <- grid::unit(pmax(0, (data$x - layout$bounds$x[1]) / span_x) * 65, "mm")
    d_right_mm <- grid::unit(pmax(0, (layout$bounds$x[2] - data$x) / span_x) * 65, "mm")
    d_bot_mm <- grid::unit(pmax(0, (data$y - layout$bounds$y[1]) / span_y) * 55, "mm")
    d_top_mm <- grid::unit(pmax(0, (layout$bounds$y[2] - data$y) / span_y) * 55, "mm")
    h_align <- layer$geom_params$title_hjust
    rotated_layer <- ((layer$aes_params$angle %||% 0) %% 360) != 0
    left_overhang <- grid::unit.pmax(zero_mm,
      (if (is.null(h_align) || rotated_layer) width / 2 else width * h_align) - d_left_mm)
    right_overhang <- grid::unit.pmax(zero_mm,
      (if (is.null(h_align) || rotated_layer) width / 2 else width * (1 - h_align)) - d_right_mm)
    bot_overhang <- grid::unit.pmax(zero_mm, height / 2 - d_bot_mm)
    top_overhang <- grid::unit.pmax(zero_mm, height / 2 - d_top_mm)
    needed[4] <- reserve(needed[4], width + clearance, left & data$nx < 0)
    needed[4] <- reserve(needed[4], left_overhang, data$nx == 0)
    needed[2] <- reserve(needed[2], width + clearance, right & data$nx > 0)
    needed[2] <- reserve(needed[2], right_overhang, data$nx == 0)
    needed[3] <- reserve(needed[3], height + clearance, bottom & data$ny < 0)
    needed[3] <- reserve(needed[3], bot_overhang, data$ny == 0)
    needed[1] <- reserve(needed[1], height + clearance, top & data$ny > 0)
    needed[1] <- reserve(needed[1], top_overhang, data$ny == 0)
  }
  if (!is_circular_layout(layout)) {
    gene_labels <- Filter(function(layer) isTRUE(layer$ideogram_gene_label_start), plot$layers)
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
      layer$geom_params$name_margin_clearance <- name_space$margin_clearance
      layer$geom_params$name_end_gap <- name_space$end_gap
      layer$geom_params$name_title_alignment <- name_space$alignment
      plot$layers[[i]] <- layer
      dimensions <- circular_text_dimensions(layer$data$label, layer$aes_params$size,
        layer$aes_params$family %||% "")
      clearance <- grid::unit(layer$ideogram_name_gap * layer$aes_params$size, "mm") +
        max(name_space$margin_clearance)
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
  if (isTRUE(attr(margin, "ideogram_auto_margin"))) {
    margin <- auto_text_space(margin, needed)
    attr(margin, "ideogram_auto_margin") <- TRUE
    plot <- plot + ggplot2::theme(plot.margin = margin)
  }
  spacing <- plot$theme$legend.box.spacing
  if (isTRUE(attr(spacing, "ideogram_auto_spacing"))) {
    spacing <- auto_text_space(spacing, needed[3])
    attr(spacing, "ideogram_auto_spacing") <- TRUE
    plot <- plot + ggplot2::theme(legend.box.spacing = spacing)
  }
  reserve_ideogram_heading_space(plot, needed)
}

auto_text_space <- function(current, needed) {
  # Keep the user's base separately from the allowance for current contents.
  base <- attr(current, "ideogram_space_base") %||% current
  result <- grid::unit.pmax(base, needed)
  attr(result, "ideogram_space_base") <- base
  result
}

as_ideogram_margin <- function(x) {
  proto <- ggplot2::margin(0, 0, 0, 0, unit = "mm")
  if ("S7_object" %in% class(proto)) {
    attr(x, "S7_class") <- attr(proto, "S7_class")
  }
  class(x) <- unique(c(
    setdiff(class(proto), c("simpleUnit", "unit", "unit_v2", "S7_object")),
    class(x),
    intersect("S7_object", class(proto))
  ))
  x
}

auto_element_margin <- function(element, side, needed) {
  if (inherits(element, "element_blank")) return(NULL)
  margin <- element$margin %||% ggplot2::margin(0, 0, 0, 0, unit = "mm")
  base <- attr(margin, "ideogram_space_base") %||% margin
  zero <- grid::unit(0, "mm")
  pad <- list(zero, zero, zero, zero)
  pad[[side]] <- needed
  result <- as_ideogram_margin(grid::unit.pmax(base, do.call(grid::unit.c, pad)))
  attr(result, "ideogram_space_base") <- base
  ggplot2::element_text(margin = result)
}

reserve_ideogram_heading_space <- function(plot, needed) {
  has_subtitle <- !is.null(plot$labels$subtitle) &&
    !inherits(plot$theme$plot.subtitle, "element_blank")
  spec <- list(
    plot.title = list(side = 3L, pad = if (has_subtitle) grid::unit(0, "mm") else needed[1]),
    plot.subtitle = list(side = 3L, pad = needed[1]),
    axis.title.x.top = list(side = 3L, pad = needed[1]),
    axis.title.x = list(side = 1L, pad = needed[3]),
    plot.caption = list(side = 1L, pad = needed[3]),
    axis.title.y = list(side = 2L, pad = needed[4]),
    axis.title.y.right = list(side = 4L, pad = needed[2])
  )
  updates <- list()
  for (field in names(spec)) {
    entry <- spec[[field]]
    updated <- auto_element_margin(plot$theme[[field]], entry$side, entry$pad)
    if (!is.null(updated)) updates[[field]] <- updated
  }
  if (!length(updates)) return(plot)
  plot + do.call(ggplot2::theme, updates)
}

circular_text_space <- function(plot) {
  layout <- ideogram_plot_layout(plot)
  needed <- grid::unit(rep(0, 4), "mm")
  if (!is_circular_layout(layout)) return(needed)
  axis_offset <- grid::unit(0, "mm")
  axis <- layout$circular_name_axis
  if (!is.null(axis)) {
    axis_offset <- grid::unit(axis$tick_length, "mm")
    theta_axis <- circular_theta(layout, axis$label_x %||% layout$chrom$.name_x)
    if (isTRUE(axis$labels) && length(axis$label)) {
      raw_dims <- circular_text_dimensions(axis$label, axis$size, axis$family)
      alpha <- (axis$angle %||% 0) * pi / 180
      axis_offset <- axis_offset + max(abs(sin(theta_axis + alpha)) * raw_dims$width +
        abs(cos(theta_axis + alpha)) * raw_dims$height) +
        grid::unit(axis$label_gap * axis$size, "mm")
    }
    nx_a <- sin(theta_axis); ny_a <- cos(theta_axis)
    axis_overhang <- list(
      pmax(ny_a, 0) * axis_offset,
      pmax(nx_a, 0) * axis_offset,
      pmax(-ny_a, 0) * axis_offset,
      pmax(-nx_a, 0) * axis_offset
    )
    for (i in seq_len(4)) needed[i] <- max(needed[i], max(axis_overhang[[i]]))
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
  max_r <- max(abs(unlist(layout$bounds)))
  for (layer in Filter(function(l) isTRUE(l$ideogram_gene_label), plot$layers)) {
    data <- layer$data
    if (!nrow(data) || max_r - max(abs(data$.ly)) > layout$chromosome_width * 2) next
    dimensions <- rotated_text_dimensions(
      circular_text_dimensions(data$.model_label, layer$aes_params$size,
        layer$aes_params$family %||% "", layer$aes_params$fontface %||% 1),
      layer$aes_params$angle %||% 0)
    theta <- circular_theta(layout, data$.lx)
    dir <- if (layout$circular$clockwise) 1 else -1
    tx <- -dir * cos(theta); ty <- dir * sin(theta)
    reach <- pmax(abs(layer$aes_params$hjust), abs(1 - layer$aes_params$hjust))
    ext_x <- abs(tx) * dimensions$width * reach + dimensions$height / 2
    ext_y <- abs(ty) * dimensions$width * reach + dimensions$height / 2
    xy <- circular_xy(layout, data$.lx, data$.ly)
    at_top <- xy$y >= max_r - layout$chromosome_width * 2 | ty > 0.5
    at_right <- xy$x >= max_r - layout$chromosome_width * 2 | tx > 0.5
    at_bottom <- xy$y <= -max_r + layout$chromosome_width * 2 | ty < -0.5
    at_left <- xy$x <= -max_r + layout$chromosome_width * 2 | tx < -0.5
    if (any(at_top)) needed[1] <- max(needed[1], max(ext_y[at_top]))
    if (any(at_right)) needed[2] <- max(needed[2], max(ext_x[at_right]))
    if (any(at_bottom)) needed[3] <- max(needed[3], max(ext_y[at_bottom]))
    if (any(at_left)) needed[4] <- max(needed[4], max(ext_x[at_left]))
  }
  needed
}
