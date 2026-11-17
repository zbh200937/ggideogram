#' @noRd
geom_chr_connection <- function(mapping = NULL, data = NULL,
    side1 = 'auto', side2 = 'auto', gap = NULL, bend = 0, curvature = NULL,
    curve_points = 32, stat = 'identity',
    position = 'identity', ..., na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
    nodes = NULL) {
  chr_pair_component(mapping, data, 'point', side1, side2, gap, bend, curvature,
    curve_points, '+', stat,
    position, list(...), na.rm, show.legend, inherit.aes, nodes, !missing(bend))
}

#' @noRd
geom_chr_synteny <- function(mapping = NULL, data = NULL,
    side1 = 'auto', side2 = 'auto', gap = NULL, bend = 0, orientation = '+',
    curvature = NULL, curve_points = 32, stat = 'identity',
    position = 'identity', ..., na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
    nodes = NULL) {
  chr_pair_component(mapping, data, 'interval', side1, side2, gap, bend, curvature,
    curve_points, orientation,
    stat, position, list(...), na.rm, show.legend, inherit.aes, nodes, !missing(bend))
}

chr_pair_component <- function(mapping, data, type, side1, side2, gap, bend,
    curvature, curve_points, orientation, stat, position, params, na.rm, show.legend,
    inherit.aes, nodes, bend_explicit) {
  if (!identical(stat, 'identity')) stopf('Chromosome pairs require stat = "identity".')
  choices <- c('auto', 'left', 'right', 'inner', 'outer', 'center')
  side1 <- canonical_chr_side(match.arg(side1, choices))
  side2 <- canonical_chr_side(match.arg(side2, choices))
  if (!is.null(gap) && (!is.numeric(gap) || !length(gap) %in% c(1, 2) || anyNA(gap) ||
      any(!is.finite(gap)) || any(gap < 0))) stopf('`gap` must contain one or two non-negative offsets.')
  if (!is.numeric(bend) || length(bend) != 1L || !is.finite(bend)) stopf('`bend` must be one finite layout distance.')
  if (!is.null(curvature) && (!is.numeric(curvature) || length(curvature) != 1L ||
      !is.finite(curvature) || curvature < 0 || curvature > 1)) stopf('`curvature` must be NULL or a number from 0 to 1.')
  curve_points <- validate_curve_points(curve_points)
  structure(list(mapping = mapping, data = data, nodes = nodes, type = type, side1 = side1,
    side2 = side2, gap = if (is.null(gap)) NULL else rep(gap, length.out = 2), bend = bend, orientation = orientation,
    bend_explicit = bend_explicit, curvature = curvature, curve_points = curve_points,
    position = position, params = params, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes), class = 'ggideogram_pair_component')
}

prepare_chr_pairs <- function(object, layout) {
  d <- object$data
  if (!is.data.frame(d)) stopf('Chromosome pairs require data-frame `data`.')
  required <- if (object$type == 'point') c('chr1', 'position1', 'chr2', 'position2') else
    c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2')
  f <- mapped_fields(d, object$mapping, required, what = 'mapping')
  source <- layout$data$karyotype
  visible <- rep(TRUE, nrow(d))
  scope <- list()
  for (j in 1:2) {
    chr <- validate_chr(f[[paste0('chr', j)]], paste0('chr', j))
    i <- match(chr, source$.chr)
    if (anyNA(i)) stopf('Chromosome pairs contain an unknown source chromosome.')
    lo <- validate_coordinate(f[[paste0(if (object$type == 'point') 'position' else 'start', j)]], 'pair start')
    hi <- if (object$type == 'point') lo else validate_coordinate(f[[paste0('end', j)]], 'pair end')
    closed <- object$type == 'interval'
    if (any(lo > hi | lo < source$.start[i] + as.numeric(closed) | hi > source$.end[i]) ||
        (closed && any(lo != floor(lo) | hi != floor(hi)))) stopf('Pair coordinates are outside source bounds or invalid.')
    d[[paste0('.pair_chr', j)]] <- chr
    d[[paste0('.pair_start', j)]] <- lo
    d[[paste0('.pair_end', j)]] <- hi
    d[[paste0('.pair_low', j)]] <- lo - as.numeric(closed)
    d[[paste0('.pair_high', j)]] <- hi
    v <- match(chr, layout$chrom$.chr)
    visible <- visible & !is.na(v)
    scope[[j]] <- list(low = layout$chrom$.start[v], high = layout$chrom$.end[v])
    if (!closed) visible <- visible & !is.na(v) & lo >= scope[[j]]$low & lo <= scope[[j]]$high
  }
  d$.pair_orientation <- if ('orientation' %in% names(object$mapping))
    recycle_semantic(rlang::eval_tidy(object$mapping$orientation, d), nrow(d), 'orientation') else
    rep(object$orientation, length.out = nrow(d))
  if (anyNA(d$.pair_orientation) || any(!d$.pair_orientation %in% c('+', '-'))) {
    stopf('Pair orientation must be explicitly + or -.')
  }
  if (object$type == 'interval') {
    from <- list(d$.pair_low1, ifelse(d$.pair_orientation == '+', d$.pair_low2, d$.pair_high2))
    to <- list(d$.pair_high1, ifelse(d$.pair_orientation == '+', d$.pair_high2, d$.pair_low2))
    lower <- rep(0, nrow(d)); upper <- rep(1, nrow(d))
    for (j in 1:2) {
      a <- (scope[[j]]$low - from[[j]]) / (to[[j]] - from[[j]])
      b <- (scope[[j]]$high - from[[j]]) / (to[[j]] - from[[j]])
      lower <- pmax(lower, pmin(a, b)); upper <- pmin(upper, pmax(a, b))
    }
    visible <- visible & !is.na(lower) & !is.na(upper) & lower < upper
    for (j in 1:2) {
      raw_f <- from[[j]] + lower * (to[[j]] - from[[j]])
      raw_t <- from[[j]] + upper * (to[[j]] - from[[j]])
      d[[paste0('.pair_from', j)]] <- pmin(pmax(raw_f, scope[[j]]$low), scope[[j]]$high)
      d[[paste0('.pair_to', j)]] <- pmin(pmax(raw_t, scope[[j]]$low), scope[[j]]$high)
    }
  }
  d$.pair_id <- seq_len(nrow(d))
  d[visible %in% TRUE, , drop = FALSE]
}

pair_endpoint_offset <- function(layout, chr, point, other, side, gap) {
  automatic_gap <- is.null(gap)
  gap <- gap %||% 0
  if (is_circular_layout(layout) && side == 'auto') {
    return(offset_chr_points_signed(layout, chr, point,
      rep(-layout$chromosome_width / 2 - layout$track_extent[['left']] - gap, length(chr))))
  }
  g <- layout$chrom[match_layout_chr(layout, chr), , drop = FALSE]
  normal <- chromosome_right_normal(layout, g)
  sign <- switch(side, left = -1, right = 1, center = 0,
    auto = ifelse((other$x - point$x) * normal$nx + (other$y - point$y) * normal$ny >= 0, 1, -1))
  if (automatic_gap && isTRUE(layout$microsynteny)) {
    gap <- microsynteny_link_gap(layout, chr, rep(sign, length.out = length(chr)))
  }
  offset_chr_points_signed(layout, chr, point,
    rep(sign, length.out = length(chr)) * (layout$chromosome_width / 2 + gap))
}

#' @method ggplot_add ggideogram_pair_component
#' @export
ggplot_add.ggideogram_pair_component <- function(object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  object <- resolve_chr_pair_nodes(object, layout)
  d <- prepare_chr_pairs(object, layout)
  if (!nrow(d)) return(plot)
  first <- project_positions_checked(layout, d$.pair_chr1,
    if (object$type == 'point') d$.pair_start1 else (d$.pair_from1 + d$.pair_to1) / 2)
  second <- project_positions_checked(layout, d$.pair_chr2,
    if (object$type == 'point') d$.pair_start2 else (d$.pair_from2 + d$.pair_to2) / 2)
  mapping <- track_mapping_without(object$mapping,
    c('chr1', 'position1', 'start1', 'end1', 'chr2', 'position2', 'start2', 'end2', 'orientation'))
  if (is_circular_layout(layout)) {
    return(add_circular_chr_pairs(object, plot, layout, d, first, second, mapping))
  }
  d$.pair_curvature <- linear_pair_curvature(object, d)
  if (object$type == 'point') {
    a <- pair_endpoint_offset(layout, d$.pair_chr1, first, second, object$side1, object$gap[1])
    b <- pair_endpoint_offset(layout, d$.pair_chr2, second, first, object$side2, object$gap[2])
    d$.pair_x <- a$x; d$.pair_y <- a$y; d$.pair_xend <- b$x; d$.pair_yend <- b$y
    if (object$bend != 0) {
      normal <- chromosome_right_normal(layout, layout$chrom[match_layout_chr(layout, d$.pair_chr1), , drop = FALSE])
      mx <- (a$x + b$x) / 2 + object$bend * normal$nx
      my <- (a$y + b$y) / 2 + object$bend * normal$ny
      d <- d[rep(seq_len(nrow(d)), each = 2), , drop = FALSE]
      d$.pair_x <- as.vector(rbind(a$x, mx)); d$.pair_y <- as.vector(rbind(a$y, my))
      d$.pair_xend <- as.vector(rbind(mx, b$x)); d$.pair_yend <- as.vector(rbind(my, b$y))
    } else if (any(d$.pair_curvature != 0)) {
      normal <- linear_pair_curve_normal(layout, d$.pair_chr1, first, a, b)
      auto_corridor <- isTRUE(layout$microsynteny) && is.null(object$gap) &&
        object$side1 == 'auto' && object$side2 == 'auto'
      adjacent <- if (auto_corridor) {
        abs(layout$index[d$.pair_chr1] - layout$index[d$.pair_chr2]) == 1L
      } else rep(FALSE, nrow(d))
      curved <- d$.pair_curvature != 0
      k_pts <- object$curve_points
      seg_counts <- ifelse(curved, k_pts - 1L, 1L)
      idx <- rep(seq_len(nrow(d)), seg_counts)
      out <- d[idx, , drop = FALSE]
      is_straight <- !curved[idx]
      out$.pair_x <- numeric(length(idx))
      out$.pair_y <- numeric(length(idx))
      out$.pair_xend <- numeric(length(idx))
      out$.pair_yend <- numeric(length(idx))
      if (any(!curved)) {
        out$.pair_x[is_straight] <- a$x[!curved]
        out$.pair_y[is_straight] <- a$y[!curved]
        out$.pair_xend[is_straight] <- b$x[!curved]
        out$.pair_yend[is_straight] <- b$y[!curved]
      }
      ci <- which(curved)
      c_mat <- sample_linear_curves(a$x[ci], a$y[ci], b$x[ci], b$y[ci],
        d$.pair_curvature[ci], k_pts, normal[ci, , drop = FALSE],
        adjacent[ci], a$y[ci], b$y[ci])
      out$.pair_x[!is_straight] <- as.vector(c_mat$x[-k_pts, , drop = FALSE])
      out$.pair_y[!is_straight] <- as.vector(c_mat$y[-k_pts, , drop = FALSE])
      out$.pair_xend[!is_straight] <- as.vector(c_mat$x[-1L, , drop = FALSE])
      out$.pair_yend[!is_straight] <- as.vector(c_mat$y[-1L, , drop = FALSE])
      d <- out
    }
    mapping <- combine_track_mapping(mapping, ggplot2::aes(x = .data$.pair_x, y = .data$.pair_y,
      xend = .data$.pair_xend, yend = .data$.pair_yend))
    geom <- ggplot2::geom_segment
  } else {
    a <- pair_endpoint_offset(layout, d$.pair_chr1, first, second, object$side1, object$gap[1])
    b <- pair_endpoint_offset(layout, d$.pair_chr2, second, first, object$side2, object$gap[2])
    ends <- list()
    for (j in 1:2) for (edge in c('from', 'to')) {
      point <- project_positions_checked(layout, d[[paste0('.pair_chr', j)]], d[[paste0('.pair_', edge, j)]])
      center <- if (j == 1) first else second
      shifted <- if (j == 1) a else b
      point$x <- point$x + shifted$x - center$x
      point$y <- point$y + shifted$y - center$y
      ends[[paste0(edge, j)]] <- point
    }
    d <- linear_pair_polygons(d, ends, object, layout, first, a, b)
    mapping <- combine_track_mapping(mapping, ggplot2::aes(x = .data$.pair_x,
      y = .data$.pair_y, group = .data$.pair_id))
    geom <- ggplot2::geom_polygon
  }
  layout$bounds$x <- range(layout$bounds$x, d$.pair_x, d$.pair_xend)
  layout$bounds$y <- range(layout$bounds$y, d$.pair_y, d$.pair_yend)
  plot <- update_plot_ideogram_layout(plot, layout)
  add_chr_pair_native_layer(plot, object, mapping, d, geom, 'identity')
}

linear_pair_polygons <- function(data, ends, object, layout, center, a, b) {
  ends <- ends[c('from1', 'from2', 'to2', 'to1')]
  x <- do.call(rbind, lapply(ends, `[[`, 'x'))
  y <- do.call(rbind, lapply(ends, `[[`, 'y'))
  if (object$bend != 0) {
    normal <- chromosome_right_normal(layout,
      layout$chrom[match_layout_chr(layout, data$.pair_chr1), , drop = FALSE])
    bend_edges <- function(m, delta) rbind(m[1, ], (m[1, ] + m[2, ]) / 2 + delta,
      m[2, ], m[3, ], (m[3, ] + m[4, ]) / 2 + delta, m[4, ])
    x <- bend_edges(x, object$bend * normal$nx)
    y <- bend_edges(y, object$bend * normal$ny)
  }

  # Straight blocks need only four vertices. Sample the curved blocks alone,
  # then expand the source table once, in its original drawing order.
  curved <- data$.pair_curvature != 0
  counts <- ifelse(curved, 2L * object$curve_points, nrow(x))
  index <- rep(seq_len(nrow(data)), counts)
  result <- data[index, , drop = FALSE]
  straight <- !curved[index]
  result$.pair_x <- numeric(length(index))
  result$.pair_y <- numeric(length(index))
  result$.pair_x[straight] <- as.vector(x[, !curved, drop = FALSE])
  result$.pair_y[straight] <- as.vector(y[, !curved, drop = FALSE])
  if (any(curved)) {
    normal <- linear_pair_curve_normal(layout, data$.pair_chr1, center, a, b)
    auto_corridor <- isTRUE(layout$microsynteny) && is.null(object$gap) &&
      object$side1 == 'auto' && object$side2 == 'auto'
    adjacent <- if (auto_corridor) {
      abs(layout$index[data$.pair_chr1] - layout$index[data$.pair_chr2]) == 1L
    } else rep(FALSE, nrow(data))
    ci <- which(curved)
    k_pts <- object$curve_points
    from_mat <- sample_linear_curves(x[1, ci], y[1, ci], x[2, ci], y[2, ci],
      data$.pair_curvature[ci], k_pts, normal[ci, , drop = FALSE],
      adjacent[ci], y[1, ci], y[2, ci])
    to_mat <- sample_linear_curves(x[3, ci], y[3, ci], x[4, ci], y[4, ci],
      data$.pair_curvature[ci], k_pts, normal[ci, , drop = FALSE],
      adjacent[ci], y[1, ci], y[2, ci])
    result$.pair_x[!straight] <- as.vector(rbind(from_mat$x, to_mat$x))
    result$.pair_y[!straight] <- as.vector(rbind(from_mat$y, to_mat$y))
  }
  result
}

add_chr_pair_native_layer <- function(plot, object, mapping, data, geom, stat) {
  if (object$type == "interval") {
    defaults <- list(fill = ideogram_colour("synteny"), alpha = 0.35)
    supplied <- union(names(object$params), names(mapping))
    if (isTRUE(object$inherit.aes)) supplied <- union(supplied, names(plot$mapping))
    object$params <- c(defaults[!names(defaults) %in% supplied], object$params)
  }
  add <- function(rows, params) {
    if (nrow(rows)) plot <<- plot + do.call(geom, c(list(mapping = mapping, data = rows,
      stat = stat, position = object$position, na.rm = object$na.rm,
      show.legend = object$show.legend, inherit.aes = object$inherit.aes), params))
  }
  arrow <- object$params$arrow
  if (object$type != 'point' || is.null(arrow) || !anyDuplicated(data$.pair_id)) {
    add(data, object$params)
    return(plot)
  }
  first <- !duplicated(data$.pair_id)
  last <- !duplicated(data$.pair_id, fromLast = TRUE)
  head_first <- first & arrow$ends %in% c(1, 3)
  head_last <- last & arrow$ends %in% c(2, 3)
  params <- object$params
  params$arrow <- NULL
  add(data[!head_first & !head_last, , drop = FALSE], params)
  arrow_at <- function(ends) grid::arrow(angle = arrow$angle, length = arrow$length,
    ends = ends, type = c('open', 'closed')[arrow$type])
  for (ends in c('first', 'last', 'both')) {
    keep <- switch(ends, first = head_first & !head_last,
      last = head_last & !head_first, both = head_first & head_last)
    params$arrow <- arrow_at(ends)
    add(data[keep, , drop = FALSE], params)
  }
  plot
}

linear_pair_curvature <- function(object, d) {
  if (object$bend != 0) return(rep(0, nrow(d)))
  if (!is.null(object$curvature)) return(rep(object$curvature, nrow(d)))
  if (object$bend_explicit) return(rep(0, nrow(d)))
  ifelse(d$.pair_chr1 == d$.pair_chr2, 0.25, 0)
}

linear_pair_curve_normal <- function(layout, chr, center, a, b) {
  dx <- b$x - a$x; dy <- b$y - a$y
  distance <- sqrt(dx^2 + dy^2)
  normal <- cbind(-dy, dx) / ifelse(distance == 0, 1, distance)
  right <- chromosome_right_normal(layout, layout$chrom[match_layout_chr(layout, chr), , drop = FALSE])
  outward <- cbind(a$x - center$x, a$y - center$y)
  on_center <- rowSums(outward^2) == 0
  outward[on_center, ] <- cbind(right$nx, right$ny)[on_center, , drop = FALSE]
  zero_dist <- distance == 0
  if (any(zero_dist)) {
    out_len <- sqrt(rowSums(outward[zero_dist, , drop = FALSE]^2))
    normal[zero_dist, ] <- outward[zero_dist, , drop = FALSE] / ifelse(out_len == 0, 1, out_len)
  }
  flip <- rowSums(normal * outward) < 0
  normal[flip, ] <- -normal[flip, , drop = FALSE]
  normal
}

sample_linear_curves <- function(ax, ay, bx, by, curvature, points, normal,
                                 has_corridor, cor_y1, cor_y2) {
  t <- seq(0, 1, length.out = points)
  w0 <- (1 - t)^2; w1 <- 2 * t * (1 - t); w2 <- t^2
  dist <- sqrt((bx - ax)^2 + (by - ay)^2)
  cx <- (ax + bx) / 2 + curvature * dist * normal[, 1]
  cy <- (ay + by) / 2 + curvature * dist * normal[, 2]
  if (any(has_corridor)) {
    lo <- pmin(cor_y1[has_corridor], cor_y2[has_corridor])
    hi <- pmax(cor_y1[has_corridor], cor_y2[has_corridor])
    cy[has_corridor] <- pmin(hi, pmax(lo, cy[has_corridor]))
  }
  list(
    x = outer(w0, ax) + outer(w1, cx) + outer(w2, bx),
    y = outer(w0, ay) + outer(w1, cy) + outer(w2, by)
  )
}

StatCartesianPair <- ggplot2::ggproto('StatCartesianPair', ggplot2::StatIdentity,
  compute_layer = function(self, data, params, layout) {
    data$ideogram_cartesian <- TRUE
    data
  })

sample_circular_curves <- function(ax, ay, bx, by, curvature, points, bend, nx, ny) {
  k <- if (curvature == 0 && bend == 0) 2L else points
  t <- seq(0, 1, length.out = k)
  w0 <- (1 - t)^2; w1 <- 2 * t * (1 - t); w2 <- t^2
  cx <- (ax + bx) / 2 * (1 - curvature) + bend * nx
  cy <- (ay + by) / 2 * (1 - curvature) + bend * ny
  list(
    x = outer(w0, ax) + outer(w1, cx) + outer(w2, bx),
    y = outer(w0, ay) + outer(w1, cy) + outer(w2, by)
  )
}

sample_circular_arcs <- function(layout, arc_from, arc_to, radius_vec, points) {
  s <- seq(0, 1, length.out = points)
  arc_mat <- outer(1 - s, arc_from) + outer(s, arc_to)
  rad_mat <- matrix(rep(radius_vec, each = points), nrow = points, ncol = length(radius_vec))
  pts <- circular_xy(layout, as.vector(arc_mat), as.vector(rad_mat))
  list(
    x = matrix(pts$x, nrow = points, ncol = length(radius_vec)),
    y = matrix(pts$y, nrow = points, ncol = length(radius_vec))
  )
}

add_circular_chr_pairs <- function(object, plot, layout, d, first, second, mapping) {
  curvature <- if (!is.null(object$curvature)) object$curvature else if (isTRUE(object$bend_explicit)) 0 else 1
  endpoints <- list()
  for (j in 1:2) {
    center <- if (j == 1) first else second
    other <- if (j == 1) second else first
    shifted <- pair_endpoint_offset(layout, d[[paste0('.pair_chr', j)]], center, other,
      object[[paste0('side', j)]], object$gap[j])
    edges <- if (object$type == 'point') 'start' else c('from', 'to')
    for (edge in edges) {
      point <- project_positions_checked(layout, d[[paste0('.pair_chr', j)]],
        d[[paste0('.pair_', edge, j)]])
      point$y <- point$y + shifted$y - center$y
      endpoints[[paste0(edge, j)]] <- point
    }
  }
  xy <- lapply(endpoints, function(p) circular_xy(layout, p$x, p$y))
  theta <- circular_theta(layout, first$x)
  nx <- sin(theta); ny <- cos(theta)
  n_pairs <- nrow(d)
  if (object$type == 'point') {
    c_mat <- sample_circular_curves(xy$start1$x, xy$start1$y, xy$start2$x, xy$start2$y,
      curvature, object$curve_points, object$bend, nx, ny)
    k_pts <- nrow(c_mat$x)
    d <- d[rep(seq_len(n_pairs), each = k_pts - 1L), , drop = FALSE]
    d$.pair_x <- as.vector(c_mat$x[-k_pts, , drop = FALSE])
    d$.pair_y <- as.vector(c_mat$y[-k_pts, , drop = FALSE])
    d$.pair_xend <- as.vector(c_mat$x[-1L, , drop = FALSE])
    d$.pair_yend <- as.vector(c_mat$y[-1L, , drop = FALSE])
  } else {
    c1 <- sample_circular_curves(xy$from1$x, xy$from1$y, xy$from2$x, xy$from2$y,
      curvature, object$curve_points, object$bend, nx, ny)
    a1 <- sample_circular_arcs(layout, endpoints$from2$x, endpoints$to2$x,
      endpoints$from2$y, object$curve_points)
    c2 <- sample_circular_curves(xy$to2$x, xy$to2$y, xy$to1$x, xy$to1$y,
      curvature, object$curve_points, object$bend, nx, ny)
    a2 <- sample_circular_arcs(layout, endpoints$to1$x, endpoints$from1$x,
      endpoints$to1$y, object$curve_points)
    poly_x <- rbind(c1$x, a1$x, c2$x, a2$x)
    poly_y <- rbind(c1$y, a1$y, c2$y, a2$y)
    d <- d[rep(seq_len(n_pairs), each = nrow(poly_x)), , drop = FALSE]
    d$.pair_x <- as.vector(poly_x)
    d$.pair_y <- as.vector(poly_y)
  }
  if (object$type == 'point') {
    mapping <- combine_track_mapping(mapping, ggplot2::aes(x = .data$.pair_x, y = .data$.pair_y,
      xend = .data$.pair_xend, yend = .data$.pair_yend))
    geom <- ggplot2::geom_segment
  } else {
    mapping <- combine_track_mapping(mapping, ggplot2::aes(x = .data$.pair_x,
      y = .data$.pair_y, group = .data$.pair_id))
    geom <- ggplot2::geom_polygon
  }
  layout <- include_ideogram_coordinates(layout, d$.pair_x, d$.pair_y,
    d$.pair_xend, d$.pair_yend, cartesian = TRUE)
  plot <- update_plot_ideogram_layout(plot, layout)
  add_chr_pair_native_layer(plot, object, mapping, d, geom, StatCartesianPair)
}
