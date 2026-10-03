#' @noRd
geom_chr_text_repel <- function(mapping = NULL, data = NULL, track = NULL,
    side = c("right", "left", "outer", "inner"), width = NULL, gap = 0.2,
    columns = 1, max_labels = Inf, overflow = c("warn", "omit"),
    seed = 1, max.iter = 10000, box.padding = 0.4, force = 1, force_pull = 1,
    segment_colour = "#626A73", segment_linewidth = 0.22, background = "white",
    stat = "identity", position = "identity", parse = FALSE, ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) stopf("Install the optional `ggrepel` package to arrange physical labels.")
  if (!identical(stat, "identity")) stopf("Physical labels require stat = 'identity'.")
  if (!is.null(width)) check_positive_layout(width, "width")
  check_nonnegative_layout(gap, "gap")
  check_positive_layout(columns, "columns")
  if (columns != floor(columns)) stopf("`columns` must be a positive integer.")
  if (!is.numeric(max_labels) || length(max_labels) != 1L || is.na(max_labels) ||
      max_labels < 0 || (is.finite(max_labels) && max_labels != floor(max_labels))) {
    stopf("`max_labels` must be a non-negative integer or Inf.")
  }
  check_nonnegative_layout(box.padding, "box.padding")
  check_positive_layout(max.iter, "max.iter")
  check_nonnegative_layout(seed, "seed")
  structure(list(mapping = mapping, data = data, track = track, columns = columns,
    side = canonical_chr_side(match.arg(side)), width = width, gap = gap,
    max_labels = max_labels, overflow = match.arg(overflow), position = position,
    params = list(...), na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes,
    controls = list(seed = seed, max.iter = max.iter, max.time = Inf,
      max.overlaps = Inf, box.padding = grid::unit(box.padding, "mm"),
      point.padding = grid::unit(0, "mm"), min.segment.length = Inf,
      force = force, force_pull = force_pull, parse = parse),
    segment_colour = segment_colour, segment_linewidth = segment_linewidth,
    background = background), class = "ggideogram_repel_component")
}

prepare_repel_labels <- function(object, layout) {
  d <- object$data
  if (!is.data.frame(d)) stopf("Physical labels require data-frame `data`.")
  f <- mapped_fields(d, object$mapping, c("chr", "position", "label"), what = "mapping")
  d$.label_chr <- validate_chr(f$chr, "chr")
  d$.label_position <- validate_coordinate(f$position, "position")
  d$.label_text <- as.character(f$label)
  d$.label_id <- paste0("L", seq_len(nrow(d)))
  d$.label_priority <- if ("priority" %in% names(object$mapping))
    recycle_semantic(rlang::eval_tidy(object$mapping$priority, d), nrow(d), "priority") else rep(0, nrow(d))
  if (!is.numeric(d$.label_priority) || anyNA(d$.label_priority) || any(!is.finite(d$.label_priority))) {
    stopf("Label `priority` must be finite numeric values.")
  }
  source <- layout$data$karyotype %||% layout$chrom
  i <- match(d$.label_chr, source$.chr)
  if (anyNA(i) || any(d$.label_position < source$.start[i] | d$.label_position > source$.end[i])) {
    stopf("Label positions are outside the source chromosome bounds.")
  }
  i <- match(d$.label_chr, layout$chrom$.chr)
  keep <- !is.na(i) & d$.label_position >= layout$chrom$.start[i] &
    d$.label_position <= layout$chrom$.end[i] & !is.na(d$.label_text) & nzchar(d$.label_text)
  d <- d[keep, , drop = FALSE]
  ids <- object$track
  if (!is.character(ids) || !length(ids) %in% 1:2 || anyNA(ids)) stopf("`track` must name one or two beside tracks.")
  specs <- lapply(ids, function(id) track_table_row(layout, id))
  sides <- unname(vapply(specs, function(s) s$side, character(1)))
  if (any(sides == "overlay") || anyDuplicated(sides)) stopf("Labels require distinct beside track sides.")
  if (length(ids) == 2 && (is.null(names(ids)) || !identical(canonical_chr_side(names(ids)), sides))) {
    stopf("Name the two tracks according to their `left`/`right` sides.")
  }
  names(ids) <- sides
  if ("side" %in% names(object$mapping)) {
    side <- canonical_chr_side(recycle_semantic(rlang::eval_tidy(object$mapping$side, object$data), nrow(object$data), "side")[keep])
    if (anyNA(side) || any(!side %in% sides)) stopf("Mapped label `side` must name an available track side.")
  } else {
    side <- rep(sides[1], nrow(d))
    for (chr in unique(d$.label_chr)) {
      rows <- which(d$.label_chr == chr)
      rows <- rows[order(d$.label_position[rows], d$.label_id[rows])]
      side[rows] <- rep(sides, length.out = length(rows))
    }
  }
  d$.label_track <- unname(ids[side])
  groups <- interaction(d$.label_chr, d$.label_track, drop = TRUE)
  d$.label_column <- 1L
  retain <- rep(TRUE, nrow(d))
  for (rows in split(seq_len(nrow(d)), groups)) {
    rows <- rows[order(-d$.label_priority[rows], d$.label_position[rows], d$.label_id[rows])]
    if (is.finite(object$max_labels) && length(rows) > object$max_labels) {
      omitted <- rows[seq_along(rows) > object$max_labels]
      retain[omitted] <- FALSE
      warning(sprintf("Label limit omitted: %s.", paste(d$.label_text[omitted], collapse = ", ")), call. = FALSE)
    }
    rows <- rows[retain[rows]]
    rows <- rows[order(d$.label_position[rows], d$.label_id[rows])]
    d$.label_column[rows] <- rep(seq_len(object$columns), length.out = length(rows))
  }
  d[retain, , drop = FALSE]
}

repel_column_scope <- function(layout, chr, spec, column, columns) {
  g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
  offsets <- spec$low_offset + c(column - 1, column) / columns * (spec$high_offset - spec$low_offset)
  if (is_circular_layout(layout)) {
    positions <- seq(g$.start, g$.end, length.out = layout$curve_points)
    points <- project_positions_checked(layout, rep(chr, length(positions) * 2),
      rep(positions, 2))
    p <- offset_chr_points_signed(layout, rep(chr, length(points$x)), points,
      rep(offsets, each = length(positions)))
    xy <- as.data.frame(circular_xy(layout, p$x, p$y))
    polygon <- rbind(xy[seq_along(positions), ], xy[rev(seq_along(positions) + length(positions)), ])
    polygon <- rbind(polygon, polygon[1, ])
    return(list(xlim = range(xy$x), ylim = range(xy$y), offset = mean(offsets), polygon = polygon))
  }
  points <- project_positions_checked(layout, rep(chr, 4), rep(c(g$.start, g$.end), each = 2))
  corners <- offset_chr_points_signed(layout, rep(chr, 4), points, rep(offsets, 2))
  list(xlim = range(corners$x), ylim = range(corners$y), offset = offsets[1])
}

chr_repel_geom <- function(kind) {
  parent <- if (kind == "leaders") ggplot2::GeomSegment else ggrepel::GeomTextRepel
  defaults <- utils::modifyList(ggplot2::GeomSegment$default_aes, ggrepel::GeomTextRepel$default_aes)
  ggplot2::ggproto(paste0("GeomChrRepel", if (kind == "leaders") "Segment" else "Text"), parent,
    required_aes = unique(c(parent$required_aes, "label", "x_orig", "y_orig", "ideogram_label_id")),
    default_aes = defaults,
    draw_panel = function(data, panel_params, coord, label_key, options, na.rm = FALSE) {
      registry <- panel_params$ideogram_labels
      if (kind == "leaders") {
        state <- new.env(parent = emptyenv())
        state$data <- data; state$panel <- panel_params
        state$coord <- if (isTRUE(options$circular)) ggplot2::ggproto(NULL, ggplot2::CoordCartesian) else coord
        state$options <- options; state$cache <- NULL
        registry[[label_key]] <- state
      } else state <- registry[[label_key]]
      grid::gTree(state = state, kind = kind, cl = "chr_repel_render")
    })
}

#' @method ggplot_add ggideogram_repel_component
#' @export
ggplot_add.ggideogram_repel_component <- function(object, plot, ...) {
  if (is.null(object$track)) {
    declared <- declare_chr_labels_track(plot, object)
    plot <- declared$plot; object$track <- declared$track
  }
  layout <- ideogram_plot_layout(plot)
  d <- prepare_repel_labels(object, layout)
  if (!nrow(d)) return(plot)
  groups <- split(d, interaction(d$.label_chr, d$.label_track, d$.label_column, drop = TRUE))
  geoms <- list(leaders = chr_repel_geom("leaders"), labels = chr_repel_geom("labels"))
  layers <- list(leaders = list(), labels = list())
  for (i in seq_along(groups)) {
    rows <- groups[[i]]
    chr <- rows$.label_chr[1]
    spec <- track_table_row(layout, rows$.label_track[1])
    scope <- repel_column_scope(layout, chr, spec, rows$.label_column[1], object$columns)
    p <- project_positions_checked(layout, rows$.label_chr, rows$.label_position)
    anchor <- offset_chr_points_signed(layout, rows$.label_chr, p,
      rep(if (spec$side == "right") layout$chromosome_width / 2 else -layout$chromosome_width / 2, nrow(rows)))
    target <- offset_chr_points_signed(layout, rows$.label_chr, p, rep(scope$offset, nrow(rows)))
    if (is_circular_layout(layout)) {
      anchor <- circular_xy(layout, anchor$x, anchor$y)
      target <- circular_xy(layout, target$x, target$y)
    }
    rows$.anchor_x <- anchor$x; rows$.anchor_y <- anchor$y
    rows$.label_x <- target$x; rows$.label_y <- target$y
    mapping <- track_mapping_without(object$mapping, c("chr", "position", "side", "priority"))
    mapping <- combine_track_mapping(mapping, ggplot2::aes(
      x = .data$.label_x, y = .data$.label_y, x_orig = .data$.anchor_x,
      y_orig = .data$.anchor_y, ideogram_label_id = .data$.label_id))
    defaults <- if (is.null(mapping$size)) list(size = ideogram_text_size("label")) else list()
    params <- utils::modifyList(defaults, object$params)
    align <- if (layout$orientation == "vertical") "hjust" else "vjust"
    if (is_circular_layout(layout)) {
      params <- utils::modifyList(list(hjust = 0.5, vjust = 0.5), params)
    } else if (!align %in% names(params) && !align %in% names(mapping)) {
      params[[align]] <- if (layout$orientation == "vertical") as.numeric(spec$side == "left") else as.numeric(spec$side == "right")
    }
    key_glyph <- params$key_glyph; params$key_glyph <- NULL
    options <- list(controls = utils::modifyList(object$controls,
      list(direction = if (is_circular_layout(layout)) "both" else if (layout$orientation == "vertical") "y" else "x",
        xlim = scope$xlim, ylim = scope$ylim)),
      source = rows[c(".label_id", ".label_text", ".label_position", ".label_priority")],
      overflow = object$overflow, background = object$background,
      circular = is_circular_layout(layout), scope_polygon = scope$polygon,
      size.unit = object$text_params$size.unit %||% "mm",
      segment_colour = object$segment_colour, segment_linewidth = object$segment_linewidth)
    for (kind in names(geoms)) {
      m <- mapping
      if (kind == "leaders") m <- combine_track_mapping(m, ggplot2::aes(xend = .data$.label_x, yend = .data$.label_y))
      layers[[kind]][[i]] <- ggplot2::layer(data = rows, mapping = m, stat = "identity",
        geom = geoms[[kind]], position = object$position,
        show.legend = if (kind == "leaders") FALSE else object$show.legend,
        inherit.aes = object$inherit.aes, key_glyph = key_glyph,
        params = c(list(na.rm = object$na.rm, label_key = paste0("labels_", length(plot$layers), "_", i),
          options = options), params))
    }
  }
  plot + c(layers$leaders, layers$labels)
}

chr_repel_text_grobs <- function(grob) {
  if (inherits(grob, "text") && grepl("^textrepelgrob[0-9]+$", grob$name)) return(list(grob))
  if (!length(grob$children)) return(list())
  unlist(lapply(grob$children, chr_repel_text_grobs), recursive = FALSE)
}

chr_repel_boxes <- function(engine) {
  text <- chr_repel_text_grobs(engine)
  if (!length(text)) return(data.frame())
  do.call(rbind, lapply(text, function(tg) data.frame(
    index = as.integer(sub("textrepelgrob", "", tg$name)),
    left = grid::convertX(grid::grobX(tg, "west"), "native", TRUE),
    right = grid::convertX(grid::grobX(tg, "east"), "native", TRUE),
    bottom = grid::convertY(grid::grobY(tg, "south"), "native", TRUE),
    top = grid::convertY(grid::grobY(tg, "north"), "native", TRUE),
    fontsize = tg$gp$fontsize)))
}

chr_repel_conflicts <- function(boxes, bounds) {
  if (!nrow(boxes)) return(integer())
  tolerance <- 1e-6
  bad <- which(boxes$left < bounds$xlim[1] - tolerance | boxes$right > bounds$xlim[2] + tolerance |
    boxes$bottom < bounds$ylim[1] - tolerance | boxes$top > bounds$ylim[2] + tolerance)
  if (nrow(boxes) > 1L) for (i in seq_len(nrow(boxes) - 1L)) {
    j <- seq.int(i + 1L, nrow(boxes))
    collision <- pmin(boxes$right[i], boxes$right[j]) - pmax(boxes$left[i], boxes$left[j]) > tolerance &
      pmin(boxes$top[i], boxes$top[j]) - pmax(boxes$bottom[i], boxes$bottom[j]) > tolerance
    if (any(collision)) bad <- c(bad, i, j[collision])
  }
  unique(bad)
}

rectangle_inside_polygon <- function(box, polygon) {
  corners <- data.frame(x = c(box$left, box$right, box$right, box$left),
    y = c(box$bottom, box$bottom, box$top, box$top))
  if (!all(points_inside_polygon(corners$x, corners$y, polygon))) return(FALSE)
  x1 <- utils::head(polygon$x, -1L); x2 <- utils::tail(polygon$x, -1L)
  y1 <- utils::head(polygon$y, -1L); y2 <- utils::tail(polygon$y, -1L)
  tolerance <- 1e-6
  left <- box$left + tolerance; right <- box$right - tolerance
  bottom <- box$bottom + tolerance; top <- box$top - tolerance
  if (left >= right || bottom >= top) return(TRUE)
  dx <- x2 - x1; dy <- y2 - y1
  tx1 <- (left - x1) / dx; tx2 <- (right - x1) / dx
  ty1 <- (bottom - y1) / dy; ty2 <- (top - y1) / dy
  near <- pmax(0, pmin(tx1, tx2), pmin(ty1, ty2))
  far <- pmin(1, pmax(tx1, tx2), pmax(ty1, ty2))
  stationary_outside <- (dx == 0 & (x1 <= left | x1 >= right)) |
    (dy == 0 & (y1 <= bottom | y1 >= top))
  !any(near < far & !stationary_outside, na.rm = TRUE)
}

chr_repel_sector_candidates <- function(boxes, polygon) {
  n <- (nrow(polygon) - 1L) / 2L
  edges <- list(c(1L, 2L * n), c(n, n + 1L))
  inside <- list(c(2L, 2L * n - 1L), c(n - 1L, n + 2L))
  corners <- cbind(x = c(boxes$left, boxes$right, boxes$right, boxes$left),
    y = c(boxes$bottom, boxes$bottom, boxes$top, boxes$top))
  normals <- matrix(0, 2L, 2L)
  distance <- numeric(2L)
  for (i in seq_along(edges)) {
    edge <- as.matrix(polygon[edges[[i]], c("x", "y")])
    centre <- colMeans(edge)
    delta <- edge[2L, ] - edge[1L, ]
    normal <- c(-delta[2L], delta[1L]) / sqrt(sum(delta^2))
    if (sum(normal * (colMeans(polygon[inside[[i]], c("x", "y")]) - centre)) < 0) {
      normal <- -normal
    }
    normals[i, ] <- normal
    distance[i] <- 1e-6 - min(sweep(corners, 2L, centre) %*% normal)
  }
  candidates <- lapply(which(distance > 0), function(i) normals[i, ] * distance[i])
  if (length(candidates) && abs(det(normals)) > 1e-8) {
    candidates <- c(candidates, list(as.numeric(solve(normals, distance))))
  }
  candidates[order(vapply(candidates, function(x) sum(x^2), numeric(1)))]
}

chr_repel_fit_sector <- function(engine, boxes, polygon, bounds) {
  if (all(vapply(seq_len(nrow(boxes)), function(i)
      rectangle_inside_polygon(boxes[i, ], polygon), logical(1))) ||
      length(chr_repel_conflicts(boxes, list(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))))) {
    return(list(engine = engine, boxes = boxes))
  }
  mm <- boxes
  mm$left <- grid::convertX(grid::unit(boxes$left, "native"), "mm", TRUE)
  mm$right <- grid::convertX(grid::unit(boxes$right, "native"), "mm", TRUE)
  mm$bottom <- grid::convertY(grid::unit(boxes$bottom, "native"), "mm", TRUE)
  mm$top <- grid::convertY(grid::unit(boxes$top, "native"), "mm", TRUE)
  polygon_mm <- data.frame(
    x = grid::convertX(grid::unit(polygon$x, "native"), "mm", TRUE),
    y = grid::convertY(grid::unit(polygon$y, "native"), "mm", TRUE))
  text <- chr_repel_text_grobs(engine)
  for (shift in chr_repel_sector_candidates(mm, polygon_mm)) {
    candidate <- engine
    for (tg in text) {
      candidate <- grid::editGrob(candidate, grid::gPath(tg$name),
        x = tg$x + grid::unit(shift[1L], "mm"),
        y = tg$y + grid::unit(shift[2L], "mm"))
    }
    moved <- chr_repel_boxes(candidate)
    if (!length(chr_repel_conflicts(moved, bounds)) &&
        all(vapply(seq_len(nrow(moved)), function(i)
          rectangle_inside_polygon(moved[i, ], polygon), logical(1)))) {
      return(list(engine = candidate, boxes = moved))
    }
  }
  list(engine = engine, boxes = boxes)
}

resolve_chr_repel <- function(state) {
  device <- c(grDevices::dev.cur(), grDevices::dev.size("px"),
    grid::convertWidth(grid::unit(1, "npc"), "mm", TRUE), grid::convertHeight(grid::unit(1, "npc"), "mm", TRUE))
  if (!is.null(state$cache) && identical(device, state$device)) return(state$cache)
  data <- locus_text_mm(state$data, state$options$size.unit)
  omitted <- character()
  repeat {
    engine <- do.call(ggrepel::GeomTextRepel$draw_panel,
      c(list(data = data, panel_scales = state$panel, coord = state$coord), state$options$controls))
    engine <- grid::forceGrob(engine)
    boxes <- chr_repel_boxes(engine)
    scope <- state$coord$transform(data.frame(x = state$options$controls$xlim,
      y = state$options$controls$ylim), state$panel)
    bounds <- list(xlim = range(scope$x), ylim = range(scope$y))
    if (!is.null(state$options$scope_polygon) && nrow(boxes)) {
      polygon <- state$coord$transform(state$options$scope_polygon, state$panel)
      fitted <- chr_repel_fit_sector(engine, boxes, polygon, bounds)
      engine <- fitted$engine; boxes <- fitted$boxes
    }
    bad <- chr_repel_conflicts(boxes, bounds)
    if (!is.null(state$options$scope_polygon) && nrow(boxes)) {
      outside <- vapply(seq_len(nrow(boxes)), function(i) {
        !rectangle_inside_polygon(boxes[i, ], polygon)
      }, logical(1))
      bad <- unique(c(bad, which(outside)))
    }
    if (!length(bad) || state$options$overflow == "warn" || !nrow(boxes)) break
    ids <- as.character(engine$data$ideogram_label_id[boxes$index[bad]])
    source <- state$options$source
    priority <- source$.label_priority[match(ids, source$.label_id)]
    drop <- ids[order(priority, ids)][1]
    omitted <- c(omitted, drop)
    data <- data[as.character(data$ideogram_label_id) != drop, , drop = FALSE]
    if (!nrow(data)) { boxes <- data.frame(); engine <- grid::nullGrob(); break }
  }
  source <- state$options$source
  if (length(bad) && state$options$overflow == "warn") {
    ids <- as.character(engine$data$ideogram_label_id[boxes$index[bad]])
    warning(sprintf("Labels do not fit their column: %s.",
      paste(source$.label_text[match(ids, source$.label_id)], collapse = ", ")), call. = FALSE)
  }
  if (length(omitted)) warning(sprintf("Physical label layout omitted: %s.",
    paste(source$.label_text[match(omitted, source$.label_id)], collapse = ", ")), call. = FALSE)
  positions <- boxes
  if (nrow(boxes)) {
    positions$id <- as.character(engine$data$ideogram_label_id[boxes$index])
    positions$anchor_x <- engine$data$x[boxes$index]
    positions$anchor_y <- engine$data$y[boxes$index]
    positions$position <- source$.label_position[match(positions$id, source$.label_id)]
  }
  state$device <- device
  state$cache <- list(engine = engine, positions = positions, omitted = omitted)
  state$cache
}

#' @method makeContent chr_repel_render
#' @importFrom grid makeContent
#' @export
makeContent.chr_repel_render <- function(x) {
  solved <- resolve_chr_repel(x$state)
  boxes <- solved$positions
  if (!nrow(boxes)) return(grid::setChildren(x, grid::gList()))
  identity <- ggplot2::ggproto(NULL, ggplot2::CoordCartesian,
    transform = function(data, panel_params) data)
  if (x$kind == "leaders") {
    d <- data.frame(x = boxes$anchor_x, y = boxes$anchor_y,
      xend = pmin(pmax(boxes$anchor_x, boxes$left), boxes$right),
      yend = pmin(pmax(boxes$anchor_y, boxes$bottom), boxes$top),
      colour = x$state$options$segment_colour, linewidth = x$state$options$segment_linewidth,
      alpha = NA_real_, linetype = 1, group = 1)
    grob <- ggplot2::GeomSegment$draw_panel(d, x$state$panel, identity)
  } else {
    background <- grid::nullGrob()
    if (!is.na(x$state$options$background)) {
      d <- data.frame(xmin = boxes$left, xmax = boxes$right,
        ymin = boxes$bottom, ymax = boxes$top, fill = x$state$options$background,
        colour = NA_character_, linewidth = 0, linetype = 1, alpha = NA_real_, group = 1)
      background <- ggplot2::GeomRect$draw_panel(d, x$state$panel, identity)
    }
    grob <- grid::grobTree(background, solved$engine)
  }
  x$positions <- boxes
  grid::setChildren(x, grid::gList(grob))
}
