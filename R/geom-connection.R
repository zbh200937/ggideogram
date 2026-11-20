#' @noRd
geom_chr_connection <- function(mapping = NULL, data = NULL,
    side1 = 'auto', side2 = 'auto', gap = NULL, bend = 0, curvature = NULL,
    curve_points = NULL, curve_style = c('auto', 'straight', 'quadratic', 'bezier'),
    unknown_chr = c('error', 'drop'),
    stat = 'identity',
    position = 'identity', ..., na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
    nodes = NULL) {
  chr_pair_component(mapping, data, 'point', side1, side2, gap, bend, curvature,
    curve_points, !missing(curve_points) && !is.null(curve_points),
    match.arg(curve_style), match.arg(unknown_chr), '+', stat,
    position, list(...), na.rm, show.legend, inherit.aes, nodes, !missing(bend))
}

#' @noRd
geom_chr_synteny <- function(mapping = NULL, data = NULL,
    side1 = 'auto', side2 = 'auto', gap = NULL, bend = 0, orientation = '+',
    curvature = NULL, curve_points = NULL,
    curve_style = c('auto', 'straight', 'quadratic', 'bezier'),
    unknown_chr = c('error', 'drop'),
    stat = 'identity',
    position = 'identity', ..., na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
    nodes = NULL) {
  chr_pair_component(mapping, data, 'interval', side1, side2, gap, bend, curvature,
    curve_points, !missing(curve_points) && !is.null(curve_points),
    match.arg(curve_style), match.arg(unknown_chr), orientation,
    stat, position, list(...), na.rm, show.legend, inherit.aes, nodes, !missing(bend))
}

chr_pair_component <- function(mapping, data, type, side1, side2, gap, bend,
    curvature, curve_points, curve_points_explicit, curve_style, unknown_chr,
    orientation, stat, position, params, na.rm, show.legend,
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
  curve_points <- if (is.null(curve_points)) 32L else validate_curve_points(curve_points)
  structure(list(mapping = mapping, data = data, nodes = nodes, type = type, side1 = side1,
    side2 = side2, gap = if (is.null(gap)) NULL else rep(gap, length.out = 2), bend = bend,
    orientation = orientation, bend_explicit = bend_explicit, curvature = curvature,
    curve_points = curve_points, curve_points_explicit = isTRUE(curve_points_explicit),
    curve_style = curve_style, unknown_chr = unknown_chr,
    position = position, params = params, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes), class = 'ggideogram_pair_component')
}

prepare_chr_pairs <- function(object, layout) {
  d <- object$data
  if (!is.data.frame(d)) stopf('Chromosome pairs require data-frame `data`.')
  if (!nrow(d)) return(d)
  required <- if (object$type == 'point') c('chr1', 'position1', 'chr2', 'position2') else
    c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2')
  f <- mapped_fields(d, object$mapping, required, what = 'mapping')
  source <- layout$data$karyotype
  drop_unknown <- identical(object$unknown_chr %||% 'error', 'drop')
  chr1 <- validate_chr(f$chr1, 'chr1')
  chr2 <- validate_chr(f$chr2, 'chr2')
  i1 <- match(chr1, source$.chr)
  i2 <- match(chr2, source$.chr)
  if (drop_unknown) {
    known_mask <- !is.na(i1) & !is.na(i2)
    if (!any(known_mask)) return(d[0, , drop = FALSE])
    if (!all(known_mask)) {
      d <- d[known_mask, , drop = FALSE]
      f <- lapply(f, `[`, known_mask)
      chr1 <- chr1[known_mask]; chr2 <- chr2[known_mask]
      i1 <- i1[known_mask]; i2 <- i2[known_mask]
    }
  } else if (anyNA(i1) || anyNA(i2)) {
    stopf('Chromosome pairs contain an unknown source chromosome.')
  }

  visible <- rep(TRUE, nrow(d))
  scope <- list()
  chrs <- list(chr1, chr2)
  src_idx <- list(i1, i2)
  closed <- object$type == 'interval'
  for (j in 1:2) {
    chr <- chrs[[j]]
    i <- src_idx[[j]]
    lo <- validate_coordinate(f[[paste0(if (object$type == 'point') 'position' else 'start', j)]], 'pair start')
    hi <- if (object$type == 'point') lo else validate_coordinate(f[[paste0('end', j)]], 'pair end')
    if (any(lo > hi | lo < source$.start[i] + as.numeric(closed) | hi > source$.end[i]) ||
        (closed && any(lo != floor(lo) | hi != floor(hi)))) {
      stopf('Pair coordinates are outside source bounds or invalid.')
    }
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

# Extract column names referenced by an aesthetic mapping to avoid duplicating
# unused source metadata columns 32-128x when expanding dense link tables.
mapping_referenced_columns <- function(mapping, col_names) {
  if (!length(mapping)) return(character(0))
  expr_vars <- function(expr) {
    if (rlang::is_quosure(expr)) expr <- rlang::quo_get_expr(expr)
    if (is.null(expr)) return(character(0))
    v <- all.vars(expr)
    # Also detect .data[["col"]] and .data$col patterns
    if (is.call(expr)) {
      if ((identical(expr[[1L]], as.name("$")) || identical(expr[[1L]], as.name("[["))) &&
          identical(expr[[2L]], as.name(".data"))) {
        rhs <- expr[[3L]]
        if (is.character(rhs)) v <- c(v, rhs) else if (is.symbol(rhs)) v <- c(v, as.character(rhs))
      } else {
        for (i in seq_along(expr)) v <- c(v, expr_vars(expr[[i]]))
      }
    }
    v
  }
  intersect(unique(unlist(lapply(mapping, expr_vars), use.names = FALSE)), col_names)
}

prune_pair_data_for_expansion <- function(d, mapping, n_pairs) {
  if (n_pairs <= 200L) return(d)
  keep <- unique(c('.pair_id', mapping_referenced_columns(mapping, names(d))))
  d[, keep, drop = FALSE]
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
  proj <- (other$x - point$x) * normal$nx + (other$y - point$y) * normal$ny
  auto_sign <- ifelse(proj > 1e-9, 1, ifelse(proj < -1e-9, -1, 1))
  if (is_serial_layout(layout) && side == 'auto') {
    # Same-row serial arcs arch outward away from the opposite row(s):
    # top row arches upward (left side, -1), bottom row arches downward (right side, +1).
    same_row <- abs(proj) <= 1e-9
    if (any(same_row)) {
      max_r <- max(layout$chrom$.row)
      auto_sign[same_row] <- ifelse(max_r > 0L & g$.row[same_row] == max_r, 1, -1)
    }
  }
  sign <- switch(side, left = -1, right = 1, center = 0, auto = auto_sign)
  if (automatic_gap && isTRUE(layout$microsynteny)) {
    gap <- microsynteny_link_gap(layout, chr, rep(sign, length.out = length(chr)))
  }
  offset_chr_points_signed(layout, chr, point,
    rep(sign, length.out = length(chr)) * (layout$chromosome_width / 2 + gap))
}

resolve_linear_pair_style <- function(object, layout, d, first, second) {
  style <- object$curve_style %||% 'auto'
  n <- nrow(d)
  if (object$bend != 0 || isTRUE(object$bend_explicit) && is.null(object$curvature) && style == 'auto') {
    return(list(mode = rep('straight', n), curvature = rep(0, n)))
  }
  same_chr <- d$.pair_chr1 == d$.pair_chr2
  vertical <- identical(layout$orientation, 'vertical')
  same_line <- if (vertical) abs(first$x - second$x) <= 1e-8 else abs(first$y - second$y) <= 1e-8

  if (style == 'straight') {
    return(list(mode = rep('straight', n), curvature = rep(0, n)))
  }
  if (style == 'bezier') {
    tens <- object$curvature %||% 0.5
    if (tens == 0) return(list(mode = rep('straight', n), curvature = rep(0, n)))
    mode <- ifelse(same_line, 'quadratic', 'bezier')
    curv <- ifelse(same_line, min(tens, 0.35), tens)
    return(list(mode = mode, curvature = curv))
  }
  if (style == 'quadratic') {
    curv <- rep(object$curvature %||% 0.25, n)
    mode <- ifelse(curv == 0, 'straight', 'quadratic')
    return(list(mode = mode, curvature = curv))
  }
  # style == 'auto'
  if (is_serial_layout(layout)) {
    if (!is.null(object$curvature) && object$curvature == 0) {
      return(list(mode = rep('straight', n), curvature = rep(0, n)))
    }
    tens <- object$curvature %||% 0.5
    mode <- ifelse(same_line, 'quadratic', 'bezier')
    curv <- ifelse(same_line, min(tens, 0.28), tens)
    return(list(mode = mode, curvature = curv))
  }
  curv <- if (!is.null(object$curvature)) rep(object$curvature, n) else ifelse(same_chr, 0.25, 0)
  mode <- ifelse(curv == 0, 'straight', 'quadratic')
  list(mode = mode, curvature = curv)
}

resolve_adaptive_curve_points <- function(object, n_pairs) {
  if (isTRUE(object$curve_points_explicit)) return(object$curve_points)
  if (n_pairs > 2000L) return(16L)
  if (n_pairs > 200L) return(24L)
  object$curve_points
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
  style_info <- resolve_linear_pair_style(object, layout, d, first, second)
  d$.pair_curvature <- style_info$curvature
  pair_mode <- style_info$mode
  n_pairs <- nrow(d)
  source_d <- d

  if (object$type == 'point') {
    a <- pair_endpoint_offset(layout, d$.pair_chr1, first, second, object$side1, object$gap[1])
    b <- pair_endpoint_offset(layout, d$.pair_chr2, second, first, object$side2, object$gap[2])
    d$.pair_x <- a$x; d$.pair_y <- a$y; d$.pair_xend <- b$x; d$.pair_yend <- b$y
    if (object$bend != 0) {
      normal <- chromosome_right_normal(layout, layout$chrom[match_layout_chr(layout, d$.pair_chr1), , drop = FALSE])
      mx <- (a$x + b$x) / 2 + object$bend * normal$nx
      my <- (a$y + b$y) / 2 + object$bend * normal$ny
      base_d <- prune_pair_data_for_expansion(d, mapping, n_pairs)
      d <- base_d[rep(seq_len(n_pairs), each = 2), , drop = FALSE]
      d$.pair_x <- as.vector(rbind(a$x, mx)); d$.pair_y <- as.vector(rbind(a$y, my))
      d$.pair_xend <- as.vector(rbind(mx, b$x)); d$.pair_yend <- as.vector(rbind(my, b$y))
    } else if (any(pair_mode != 'straight')) {
      k_pts <- resolve_adaptive_curve_points(object, n_pairs)
      curved <- pair_mode != 'straight'
      seg_counts <- ifelse(curved, k_pts - 1L, 1L)
      idx <- rep(seq_len(n_pairs), seg_counts)
      base_d <- prune_pair_data_for_expansion(d, mapping, n_pairs)
      out <- base_d[idx, , drop = FALSE]
      out$.pair_x <- numeric(length(idx))
      out$.pair_y <- numeric(length(idx))
      out$.pair_xend <- numeric(length(idx))
      out$.pair_yend <- numeric(length(idx))
      is_st <- pair_mode == 'straight'
      if (any(is_st)) {
        m_st <- is_st[idx]
        out$.pair_x[m_st] <- a$x[is_st]
        out$.pair_y[m_st] <- a$y[is_st]
        out$.pair_xend[m_st] <- b$x[is_st]
        out$.pair_yend[m_st] <- b$y[is_st]
      }
      is_q <- pair_mode == 'quadratic'
      if (any(is_q)) {
        normal <- linear_pair_curve_normal(layout, source_d$.pair_chr1, first, a, b)
        auto_corridor <- isTRUE(layout$microsynteny) && is.null(object$gap) &&
          object$side1 == 'auto' && object$side2 == 'auto'
        adjacent <- if (auto_corridor) {
          abs(layout$index[source_d$.pair_chr1] - layout$index[source_d$.pair_chr2]) == 1L
        } else rep(FALSE, n_pairs)
        qi <- which(is_q)
        c_mat <- sample_linear_curves(a$x[qi], a$y[qi], b$x[qi], b$y[qi],
          source_d$.pair_curvature[qi], k_pts, normal[qi, , drop = FALSE],
          adjacent[qi], a$y[qi], b$y[qi])
        m_q <- is_q[idx]
        out$.pair_x[m_q] <- as.vector(c_mat$x[-k_pts, , drop = FALSE])
        out$.pair_y[m_q] <- as.vector(c_mat$y[-k_pts, , drop = FALSE])
        out$.pair_xend[m_q] <- as.vector(c_mat$x[-1L, , drop = FALSE])
        out$.pair_yend[m_q] <- as.vector(c_mat$y[-1L, , drop = FALSE])
      }
      is_bz <- pair_mode == 'bezier'
      if (any(is_bz)) {
        bi <- which(is_bz)
        b_mat <- sample_cubic_s_curves(a$x[bi], a$y[bi], b$x[bi], b$y[bi],
          source_d$.pair_curvature[bi], k_pts, identical(layout$orientation, 'vertical'))
        m_bz <- is_bz[idx]
        out$.pair_x[m_bz] <- as.vector(b_mat$x[-k_pts, , drop = FALSE])
        out$.pair_y[m_bz] <- as.vector(b_mat$y[-k_pts, , drop = FALSE])
        out$.pair_xend[m_bz] <- as.vector(b_mat$x[-1L, , drop = FALSE])
        out$.pair_yend[m_bz] <- as.vector(b_mat$y[-1L, , drop = FALSE])
      }
      d <- out
    }
    if (n_pairs > 200L) attr(d, 'source_pairs') <- source_d
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
    d <- linear_pair_polygons(d, ends, object, layout, first, a, b, pair_mode, mapping)
    if (n_pairs > 200L) attr(d, 'source_pairs') <- source_d
    mapping <- combine_track_mapping(mapping, ggplot2::aes(x = .data$.pair_x,
      y = .data$.pair_y, group = .data$.pair_id))
    geom <- ggplot2::geom_polygon
  }
  layout$bounds$x <- range(layout$bounds$x, d$.pair_x, d$.pair_xend)
  layout$bounds$y <- range(layout$bounds$y, d$.pair_y, d$.pair_yend)
  plot <- update_plot_ideogram_layout(plot, layout)
  add_chr_pair_native_layer(plot, object, mapping, d, geom, 'identity')
}

linear_pair_polygons <- function(data, ends, object, layout, center, a, b,
                                 pair_mode, mapping) {
  n_pairs <- nrow(data)
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

  k_pts <- resolve_adaptive_curve_points(object, n_pairs)
  curved <- pair_mode != 'straight'
  counts <- ifelse(curved, 2L * k_pts, nrow(x))
  index <- rep(seq_len(n_pairs), counts)
  base_d <- prune_pair_data_for_expansion(data, mapping, n_pairs)
  result <- base_d[index, , drop = FALSE]
  result$.pair_x <- numeric(length(index))
  result$.pair_y <- numeric(length(index))

  is_st <- pair_mode == 'straight'
  if (any(is_st)) {
    m_st <- is_st[index]
    result$.pair_x[m_st] <- as.vector(x[, is_st, drop = FALSE])
    result$.pair_y[m_st] <- as.vector(y[, is_st, drop = FALSE])
  }
  is_q <- pair_mode == 'quadratic'
  if (any(is_q)) {
    normal <- linear_pair_curve_normal(layout, data$.pair_chr1, center, a, b)
    auto_corridor <- isTRUE(layout$microsynteny) && is.null(object$gap) &&
      object$side1 == 'auto' && object$side2 == 'auto'
    adjacent <- if (auto_corridor) {
      abs(layout$index[data$.pair_chr1] - layout$index[data$.pair_chr2]) == 1L
    } else rep(FALSE, n_pairs)
    ci <- which(is_q)
    from_mat <- sample_linear_curves(x[1, ci], y[1, ci], x[2, ci], y[2, ci],
      data$.pair_curvature[ci], k_pts, normal[ci, , drop = FALSE],
      adjacent[ci], y[1, ci], y[2, ci])
    to_mat <- sample_linear_curves(x[3, ci], y[3, ci], x[4, ci], y[4, ci],
      data$.pair_curvature[ci], k_pts, normal[ci, , drop = FALSE],
      adjacent[ci], y[1, ci], y[2, ci])
    m_q <- is_q[index]
    result$.pair_x[m_q] <- as.vector(rbind(from_mat$x, to_mat$x))
    result$.pair_y[m_q] <- as.vector(rbind(from_mat$y, to_mat$y))
  }
  is_bz <- pair_mode == 'bezier'
  if (any(is_bz)) {
    bi <- which(is_bz)
    vertical <- identical(layout$orientation, 'vertical')
    from_mat <- sample_cubic_s_curves(x[1, bi], y[1, bi], x[2, bi], y[2, bi],
      data$.pair_curvature[bi], k_pts, vertical)
    to_mat <- sample_cubic_s_curves(x[3, bi], y[3, bi], x[4, bi], y[4, bi],
      data$.pair_curvature[bi], k_pts, vertical)
    m_bz <- is_bz[index]
    result$.pair_x[m_bz] <- as.vector(rbind(from_mat$x, to_mat$x))
    result$.pair_y[m_bz] <- as.vector(rbind(from_mat$y, to_mat$y))
  }
  result
}

# Anonymous native Geom subclasses that skip redundant coord_munch() on already
# Cartesian-sampled circular links while preserving class(geom)[1] and all
# standard ggplot2 aesthetics, guides, scales and drawing behavior.
cartesian_linear_coord <- function(coord) {
  ggplot2::ggproto(NULL, coord, is_linear = function(self) TRUE)
}

cartesian_geom_polygon <- function() {
  ggplot2::ggproto(NULL, ggplot2::GeomPolygon,
    parameters = function(self, extra = FALSE) {
      ggplot2::GeomPolygon$parameters(extra = extra)
    },
    draw_panel = function(self, data, panel_params, coord, ...) {
      if ("ideogram_cartesian" %in% names(data) &&
          all(data$ideogram_cartesian %in% TRUE) && !coord$is_linear()) {
        coord <- cartesian_linear_coord(coord)
      }
      ggplot2::GeomPolygon$draw_panel(data, panel_params, coord, ...)
    })
}

cartesian_geom_segment <- function() {
  ggplot2::ggproto(NULL, ggplot2::GeomSegment,
    parameters = function(self, extra = FALSE) {
      ggplot2::GeomSegment$parameters(extra = extra)
    },
    draw_panel = function(self, data, panel_params, coord, ...) {
      if ("ideogram_cartesian" %in% names(data) &&
          all(data$ideogram_cartesian %in% TRUE) && !coord$is_linear()) {
        coord <- cartesian_linear_coord(coord)
      }
      ggplot2::GeomSegment$draw_panel(data, panel_params, coord, ...)
    })
}

add_chr_pair_native_layer <- function(plot, object, mapping, data, geom, stat,
                                      cartesian_fast = FALSE) {
  if (object$type == "interval") {
    defaults <- list(fill = ideogram_colour("synteny"), alpha = 0.35)
    supplied <- union(names(object$params), names(mapping))
    if (isTRUE(object$inherit.aes)) supplied <- union(supplied, names(plot$mapping))
    object$params <- c(defaults[!names(defaults) %in% supplied], object$params)
  }
  add <- function(rows, params) {
    if (nrow(rows)) {
      lyr <- do.call(geom, c(list(mapping = mapping, data = rows,
        stat = stat, position = object$position, na.rm = object$na.rm,
        show.legend = object$show.legend, inherit.aes = object$inherit.aes), params))
      if (isTRUE(cartesian_fast)) {
        lyr$geom <- if (object$type == 'point') cartesian_geom_segment() else cartesian_geom_polygon()
      }
      plot <<- plot + lyr
    }
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

# Cubic S-curve Bezier leaving both chromosome boundaries perpendicularly:
# horizontal rows bend vertically (control points share endpoint x);
# vertical columns bend horizontally (control points share endpoint y).
sample_cubic_s_curves <- function(ax, ay, bx, by, tension, points, vertical = FALSE) {
  t <- seq(0, 1, length.out = points)
  u <- 1 - t
  w0 <- u^3
  w1 <- 3 * u^2 * t
  w2 <- 3 * u * t^2
  w3 <- t^3
  if (vertical) {
    dx <- bx - ax
    c1x <- ax + tension * dx; c1y <- ay
    c2x <- bx - tension * dx; c2y <- by
  } else {
    dy <- by - ay
    c1x <- ax; c1y <- ay + tension * dy
    c2x <- bx; c2y <- by - tension * dy
  }
  list(
    x = outer(w0, ax) + outer(w1, c1x) + outer(w2, c2x) + outer(w3, bx),
    y = outer(w0, ay) + outer(w1, c1y) + outer(w2, c2y) + outer(w3, by)
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
  curvature <- if (!is.null(object$curvature)) object$curvature else
    if (identical(object$curve_style, 'straight')) 0 else
    if (isTRUE(object$bend_explicit)) 0 else 1
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
  source_d <- d
  base_d <- prune_pair_data_for_expansion(d, mapping, n_pairs)
  k_curve <- resolve_adaptive_curve_points(object, n_pairs)

  if (object$type == 'point') {
    c_mat <- sample_circular_curves(xy$start1$x, xy$start1$y, xy$start2$x, xy$start2$y,
      curvature, k_curve, object$bend, nx, ny)
    k_pts <- nrow(c_mat$x)
    d <- base_d[rep(seq_len(n_pairs), each = k_pts - 1L), , drop = FALSE]
    d$.pair_x <- as.vector(c_mat$x[-k_pts, , drop = FALSE])
    d$.pair_y <- as.vector(c_mat$y[-k_pts, , drop = FALSE])
    d$.pair_xend <- as.vector(c_mat$x[-1L, , drop = FALSE])
    d$.pair_yend <- as.vector(c_mat$y[-1L, , drop = FALSE])
  } else {
    if (isTRUE(object$curve_points_explicit) || n_pairs <= 200L) {
      k_arc <- object$curve_points
    } else {
      max_arc_deg <- max(
        abs(endpoints$to1$x - endpoints$from1$x),
        abs(endpoints$to2$x - endpoints$from2$x)
      ) / layout$circular$radius * 180 / pi
      k_arc <- max(2L, min(k_curve, as.integer(ceiling(max_arc_deg / 2)) + 1L))
    }
    c1 <- sample_circular_curves(xy$from1$x, xy$from1$y, xy$from2$x, xy$from2$y,
      curvature, k_curve, object$bend, nx, ny)
    a1 <- sample_circular_arcs(layout, endpoints$from2$x, endpoints$to2$x,
      endpoints$from2$y, k_arc)
    c2 <- sample_circular_curves(xy$to2$x, xy$to2$y, xy$to1$x, xy$to1$y,
      curvature, k_curve, object$bend, nx, ny)
    a2 <- sample_circular_arcs(layout, endpoints$to1$x, endpoints$from1$x,
      endpoints$to1$y, k_arc)
    poly_x <- rbind(c1$x, a1$x, c2$x, a2$x)
    poly_y <- rbind(c1$y, a1$y, c2$y, a2$y)
    d <- base_d[rep(seq_len(n_pairs), each = nrow(poly_x)), , drop = FALSE]
    d$.pair_x <- as.vector(poly_x)
    d$.pair_y <- as.vector(poly_y)
  }
  if (n_pairs > 200L) attr(d, 'source_pairs') <- source_d
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
  add_chr_pair_native_layer(plot, object, mapping, d, geom, StatCartesianPair,
    cartesian_fast = TRUE)
}
