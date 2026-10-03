GeomCircularGapText <- ggplot2::ggproto("GeomCircularGapText", GeomIdeogramAxisText,
  default_aes = local({
    mapping <- ggplot2::GeomText$default_aes
    mapping$ideogram_cartesian <- TRUE
    mapping
  }),
  draw_panel = function(data, panel_params, coord, tick_length = 1.5,
      label_gap = 0.25, gap_track = NULL, gap_title = FALSE, na.rm = FALSE) {
    text <- GeomIdeogramAxisText$draw_panel(data, panel_params, coord,
      tick_length = tick_length, label_gap = label_gap, na.rm = na.rm)
    xy <- coord$transform(data, panel_params)
    registry <- panel_params$ideogram_labels
    if (!is.environment(registry)) registry <- new.env(parent = emptyenv())
    key <- ".circular_gap_text"
    if (!exists(key, registry, inherits = FALSE)) {
      state <- new.env(parent = emptyenv())
      state$entries <- list()
      state$panel <- panel_params
      state$coord <- coord
      assign(key, state, registry)
    }
    state <- get(key, registry, inherits = FALSE)
    id <- length(state$entries) + 1L
    state$entries[[id]] <- list(text = text, data = data, xy = xy,
      track = gap_track, title = isTRUE(gap_title),
      tick_length = tick_length, label_gap = label_gap)
    state$cache <- NULL
    grid::gTree(children = grid::gList(text), state = state, id = id,
      cl = "circular_gap_text")
  })

# The templates are native text grobs from the final build. All gap layers
# register before drawing, so packing sees the actual scales and font settings.
resolve_circular_gap_text <- function(state) {
  mm <- c(grid::convertWidth(grid::unit(1, "npc"), "mm", TRUE),
    grid::convertHeight(grid::unit(1, "npc"), "mm", TRUE))
  device <- c(grDevices::dev.cur(), grDevices::dev.size("in"), mm)
  if (!is.null(state$cache) && identical(state$device, device)) return(state$cache)
  rows <- lapply(seq_along(state$entries), function(id) {
    e <- state$entries[[id]]
    text <- e$text
    dimensions <- chromosome_text_dimensions(e$data)
    x <- grid::convertX(text$x, "mm", TRUE)
    y <- grid::convertY(text$y, "mm", TRUE)
    axis <- state$coord$layout$circular_gap_axes[[e$track]]
    linewidth <- if (is.null(axis)) 0.3 else axis$linewidth
    data.frame(id = id, row = seq_len(nrow(e$data)), track = e$track,
      title = e$title, label = as.character(text$label), x = x, y = y,
      original_x = x, original_y = y,
      anchor_x = e$xy$x * mm[1], anchor_y = e$xy$y * mm[2],
      tick_x = e$xy$x * mm[1] + e$xy$nx * e$tick_length,
      tick_y = e$xy$y * mm[2] + e$xy$ny * e$tick_length,
      nx = e$xy$nx, ny = e$xy$ny,
      width = grid::convertWidth(dimensions$width, "mm", TRUE),
      height = grid::convertHeight(dimensions$height, "mm", TRUE),
      angle = rep(text$rot, length.out = nrow(e$data)),
      hjust = rep(text$hjust, length.out = nrow(e$data)),
      vjust = rep(text$vjust, length.out = nrow(e$data)),
      size = e$data$size, padding = e$label_gap * e$data$size,
      linewidth = linewidth,
      shift = 0, radial_shift = 0)
  })
  positions <- do.call(rbind, rows)
  # Keep enough room for the native leader stroke between neighbouring labels.
  positions$packing_padding <- pmax(positions$padding, max(positions$linewidth))
  boxes <- lapply(seq_len(nrow(positions)), function(i)
    circular_gap_text_corners(positions[i, ]))
  placed <- integer()
  # Numeric labels retain their radial/value location. Only the tangent into
  # the gap changes, including for explicitly rotated closing boundaries.
  radial <- positions$x * -positions$ny + positions$y * positions$nx
  numeric <- which(!positions$title)
  titles <- which(positions$title)
  order <- c(numeric[order(radial[numeric], positions$id[numeric], positions$row[numeric])],
    titles[order(radial[titles], positions$id[titles], positions$row[titles])])
  for (i in order) {
    normal <- c(positions$nx[i], positions$ny[i])
    minimum <- if (positions$title[i]) 0 else max(0,
      positions$linewidth[i] / 2 - positions$padding[i] + 1e-6)
    values <- numeric[positions$track[numeric] == positions$track[i]]
    if (positions$title[i] && length(values)) {
      outer <- max(vapply(boxes[values], function(b) max(b %*% normal), numeric(1)))
      minimum <- max(0, outer + positions$packing_padding[i] - min(boxes[[i]] %*% normal))
    }
    shift <- min(2 * positions$size[i], circular_gap_clearance(
      boxes[[i]], normal, boxes[placed], positions$packing_padding[i], minimum))
    positions$shift[i] <- shift
    positions$x[i] <- positions$x[i] + normal[1] * shift
    positions$y[i] <- positions$y[i] + normal[2] * shift
    boxes[[i]] <- sweep(boxes[[i]], 2, normal * shift, "+")
    placed <- c(placed, i)
  }
  center <- state$coord$transform(data.frame(x = 0, y = 0,
    ideogram_cartesian = TRUE), state$panel)
  center <- c(center$x * mm[1], center$y * mm[2])
  gap <- state$coord$layout$circular$label_gap
  outside <- vapply(boxes, function(b) {
    b <- sweep(b, 2, center, "-")
    delta <- (atan2(b[, 1], b[, 2]) * 180 / pi - gap$left_angle) %% 360
    any(delta > gap$angle + 1e-7 & delta < 360 - 1e-7)
  }, logical(1))
  crowded <- rep(FALSE, nrow(positions))
  if (length(boxes) > 1L) for (i in seq_len(length(boxes) - 1L)) {
    for (j in seq.int(i + 1L, length(boxes))) {
      if (circular_gap_clearance(boxes[[i]], c(1, 0), boxes[j], 0) > 1e-7) {
        crowded[c(i, j)] <- TRUE
      }
    }
  }
  texts <- segments <- leaders <- vector("list", length(state$entries))
  identity <- ggplot2::ggproto(NULL, ggplot2::CoordCartesian,
    transform = function(data, panel_params) data)
  for (id in seq_along(state$entries)) {
    e <- state$entries[[id]]
    z <- positions[positions$id == id, , drop = FALSE]
    text <- e$text
    text$x <- grid::unit(z$x / mm[1], "native")
    text$y <- grid::unit(z$y / mm[2], "native")
    texts[[id]] <- text
    moved <- which(!z$title & z$shift > 1e-7)
    if (length(moved)) {
      linewidth <- z$linewidth[1]
      paths <- lapply(moved, function(i) {
        end <- c(z$x[i], z$y[i]) - c(z$nx[i], z$ny[i]) *
          max(z$padding[i] / 2, linewidth / 2)
        path <- rbind(c(z$tick_x[i], z$tick_y[i]), end)
        data.frame(x = utils::head(path[, 1], -1) / mm[1],
          y = utils::head(path[, 2], -1) / mm[2],
          xend = utils::tail(path[, 1], -1) / mm[1],
          yend = utils::tail(path[, 2], -1) / mm[2], row = i,
          colour = e$data$colour[i], alpha = e$data$alpha[i],
          linewidth = linewidth, linetype = 1)
      })
      leader <- do.call(rbind, paths)
      leaders[[id]] <- leader
      segments[[id]] <- ggplot2::GeomSegment$draw_panel(leader, state$panel, identity)
    }
  }
  state$device <- device
  state$cache <- list(texts = texts, segments = segments, leaders = leaders,
    positions = positions, boxes = boxes, fits = !any(outside | crowded))
  if (any(outside | crowded)) warning(sprintf(
    "Circular track labels exceed the closing gap or overlap for %s; increase `opening_angle`, track spacing or the output size.",
    paste(sprintf("`%s`", unique(positions$track[outside | crowded])), collapse = ", ")), call. = FALSE)
  state$cache
}

circular_gap_text_corners <- function(z) {
  corners <- rbind(c(-z$hjust * z$width, -z$vjust * z$height),
    c((1 - z$hjust) * z$width, -z$vjust * z$height),
    c((1 - z$hjust) * z$width, (1 - z$vjust) * z$height),
    c(-z$hjust * z$width, (1 - z$vjust) * z$height))
  angle <- z$angle * pi / 180
  rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
  sweep(corners %*% rotation, 2, c(z$x, z$y), "+")
}

# Each obstacle gives a finite interval of displacements that intersects its
# measured rectangle. The first free displacement preserves the normal axis.
circular_gap_clearance <- function(box, normal, obstacles, padding, minimum = 0) {
  intervals <- lapply(obstacles, function(other) {
    edges <- rbind(box[c(2, 3, 4, 1), ] - box,
      other[c(2, 3, 4, 1), ] - other)
    axes <- cbind(-edges[, 2], edges[, 1])
    axes <- axes / sqrt(rowSums(axes^2))
    lo <- -Inf; hi <- Inf
    for (j in seq_len(nrow(axes))) {
      a <- range(box %*% axes[j, ]); b <- range(other %*% axes[j, ])
      velocity <- sum(normal * axes[j, ])
      if (abs(velocity) < 1e-12) {
        if (a[2] <= b[1] - padding || a[1] >= b[2] + padding) return(NULL)
      } else {
        interval <- sort(c((b[1] - padding - a[2]) / velocity,
          (b[2] + padding - a[1]) / velocity))
        lo <- max(lo, interval[1]); hi <- min(hi, interval[2])
      }
      if (lo >= hi) return(NULL)
    }
    c(lo, hi)
  })
  intervals <- Filter(Negate(is.null), intervals)
  if (!length(intervals)) return(minimum)
  intervals <- do.call(rbind, intervals)
  intervals <- intervals[order(intervals[, 1]), , drop = FALSE]
  shift <- minimum
  for (j in seq_len(nrow(intervals)))
    if (shift >= intervals[j, 1] - 1e-8 && shift <= intervals[j, 2] + 1e-8)
      shift <- intervals[j, 2] + 1e-8
  shift
}

#' @method makeContent circular_gap_text
#' @importFrom grid makeContent
#' @export
makeContent.circular_gap_text <- function(x) {
  solved <- resolve_circular_gap_text(x$state)
  x$positions <- solved$positions[solved$positions$id == x$id, , drop = FALSE]
  children <- grid::gList(solved$texts[[x$id]])
  if (!is.null(solved$segments[[x$id]])) children <-
    grid::gList(solved$segments[[x$id]], solved$texts[[x$id]])
  grid::setChildren(x, children)
}
