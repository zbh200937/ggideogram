#' Internal ordered label component using native text and segment geoms
#' @noRd
geom_chr_labels <- function(mapping = NULL, data = NULL, track = NULL,
    side = c("right", "left", "outer", "inner"), width = NULL, gap = 0.2,
    connection_height = 5, label_gap = 0, padding = 0.6,
    segment_colour = ideogram_colour("annotation"), segment_linewidth = ideogram_linewidth("annotation"),
    stat = "identity", position = "identity", parse = FALSE, ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  if (!identical(stat, "identity")) stopf("Ordered labels require stat = 'identity'.")
  validate_optional_marker_track(track)
  if (!is.null(width)) check_positive_layout(width, "width")
  check_nonnegative_layout(gap, "gap")
  check_nonnegative_layout(connection_height, "connection_height")
  check_nonnegative_layout(label_gap, "label_gap")
  check_nonnegative_layout(padding, "padding")
  check_nonnegative_layout(segment_linewidth, "segment_linewidth")
  if (!is.null(mapping$x) && is.null(mapping$position)) mapping$position <- mapping$x
  structure(list(mapping = mapping, data = data, track = track,
    side = canonical_chr_side(match.arg(side)), width = width, gap = gap,
    connection_height = connection_height, label_gap = label_gap,
    padding = padding, position = position, parse = parse,
    params = list(...), segment_colour = segment_colour,
    segment_linewidth = segment_linewidth, na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes),
    class = "ggideogram_labels_component")
}

declare_chr_labels_track <- function(plot, object) {
  layout <- ideogram_plot_layout(plot)
  if (!is.null(object$track)) {
    spec <- track_table_row(layout, object$track)
    if (spec$side == "overlay") stopf("Ordered labels require a beside track.")
    return(list(plot = plot, track = object$track))
  }
  source <- list(data = object$data, mapping = object$mapping,
    side = object$side, width = object$width, gap = object$gap)
  existing <- names(Filter(function(spec) identical(spec$.labels_source, source),
    layout$base_spec$tracks %||% list()))
  if (length(existing)) return(list(plot = plot, track = existing[1]))
  id <- "chr_labels"
  while (id %in% layout$tracks$id) id <- paste0(id, "_")
  width <- object$width %||% if (is_circular_layout(layout) && object$side == "left") 14 else 32
  spec <- geom_track(track = id, side = object$side, width = width, gap = object$gap)
  spec$.labels_source <- source
  plot <- ggplot_add.ggideogram_track_scope(spec, plot)
  list(plot = plot, track = id)
}

StatChrLabels <- ggplot2::ggproto("StatChrLabels", ggplot2::StatIdentity,
  required_aes = c("chr", "position", "label", "ideogram_label_id"),
  extra_params = c("na.rm", "track"),
  compute_layer = function(self, data, params, layout) {
    object <- layout$coord$layout
    spec <- track_table_row(object, params$track)
    p <- project_positions_checked(object, as.character(data$chr), data$position)
    anchor <- offset_chr_points(object, as.character(data$chr), p,
      side = spec$side, distance = object$chromosome_width / 2)
    target <- offset_chr_points_signed(object, as.character(data$chr), p,
      rep(spec$low_offset, nrow(data)))
    data$x_orig <- anchor$x; data$y_orig <- anchor$y
    data$x <- target$x; data$y <- target$y
    data
  })

chr_spread_geom <- function(kind, text_geom = ggplot2::GeomText) {
  parent <- if (kind == "leaders") ggplot2::GeomSegment else text_geom
  ggplot2::ggproto(paste0("GeomChrLabels", if (kind == "leaders") "Segment" else "Text"), parent,
    required_aes = c("x", "y", "label", "chr", "position", "ideogram_label_id"),
    default_aes = utils::modifyList(ggplot2::GeomSegment$default_aes, text_geom$default_aes),
    draw_panel = function(data, panel_params, coord, label_key, options, na.rm = FALSE) {
      registry <- panel_params$ideogram_labels
      if (kind == "leaders") {
        state <- new.env(parent = emptyenv())
        state$data <- data; state$panel <- panel_params
        state$coord <- coord; state$options <- options; state$cache <- NULL
        registry[[label_key]] <- state
      } else state <- registry[[label_key]]
      grid::gTree(state = state, kind = kind, cl = "chr_labels_render")
    })
}

#' @method ggplot_add ggideogram_labels_component
#' @export
ggplot_add.ggideogram_labels_component <- function(object, plot, ...) {
  declared <- declare_chr_labels_track(plot, object)
  plot <- declared$plot; object$track <- declared$track
  layout <- ideogram_plot_layout(plot)
  prepare <- object
  prepare$columns <- 1; prepare$max_labels <- Inf
  d <- prepare_repel_labels(prepare, layout)
  if (!nrow(d)) return(plot)
  mapping <- track_mapping_without(object$mapping, c("chr", "position", "x", "side", "priority"))
  mapping <- combine_track_mapping(mapping, ggplot2::aes(chr = .data$.label_chr,
    position = .data$.label_position, label = .data$.label_text,
    ideogram_label_id = .data$.label_id))
  defaults <- if (is.null(mapping$size)) list(size = ideogram_text_size("label")) else list()
  params <- utils::modifyList(defaults, object$params)
  key_glyph <- params$key_glyph; params$key_glyph <- NULL
  options <- list(track = object$track, source = d,
    connection_height = object$connection_height,
    label_gap = object$label_gap, padding = object$padding, parse = object$parse,
    text_geom = object$text_geom %||% ggplot2::GeomText,
    text_params = object$text_params %||% list(),
    segment_colour = object$segment_colour, segment_linewidth = object$segment_linewidth)
  layers <- lapply(c("leaders", "labels"), function(kind) ggplot2::layer(
    data = d, mapping = mapping, stat = StatChrLabels,
    geom = chr_spread_geom(kind, options$text_geom), position = object$position,
    inherit.aes = object$inherit.aes,
    show.legend = if (kind == "leaders") FALSE else object$show.legend,
    key_glyph = key_glyph,
    params = c(list(track = object$track, na.rm = object$na.rm,
      label_key = paste0("ordered_labels_", length(plot$layers)), options = options), params)))
  plot <- plot + layers
  reserve_chromosome_axis_space(plot, list())
}

ordered_label_positions <- function(target, halfwidth, lo, hi) {
  i <- order(target, seq_along(target))
  position <- target[i]; half <- halfwidth[i]
  fit <- 2 * sum(half) <= hi - lo
  if (!length(i)) return(list(position = numeric(), fits = TRUE))
  position[1] <- max(position[1], lo + half[1])
  if (length(i) > 1L) for (j in 2:length(i))
    position[j] <- max(position[j], position[j - 1L] + half[j - 1L] + half[j])
  if (fit) {
    position[length(i)] <- min(position[length(i)], hi - half[length(i)])
    if (length(i) > 1L) for (j in rev(seq_len(length(i) - 1L)))
      position[j] <- min(position[j], position[j + 1L] - half[j + 1L] - half[j])
  } else {
    position <- lo + cumsum(2 * half) - half
    position <- position + (hi - lo - 2 * sum(half)) / 2
  }
  result <- numeric(length(i)); result[i] <- position
  list(position = result, fits = fit)
}

resolve_chr_labels <- function(state) {
  mm <- c(grid::convertWidth(grid::unit(1, "npc"), "mm", TRUE),
    grid::convertHeight(grid::unit(1, "npc"), "mm", TRUE))
  device <- c(grDevices::dev.cur(), grDevices::dev.size("in"), mm)
  if (!is.null(state$cache) && identical(state$device, device)) return(state$cache)
  data <- locus_text_mm(state$data, state$options$text_params$size.unit %||% "mm")
  coord <- state$coord; layout <- coord$layout
  spec <- track_table_row(layout, state$options$track)
  dimensions <- chromosome_text_dimensions(data, state$options$parse)
  width <- grid::convertWidth(dimensions$width, "mm", TRUE)
  height <- grid::convertHeight(dimensions$height, "mm", TRUE)
  base <- coord$transform(data, state$panel)
  anchors <- coord$transform(data.frame(x = data$x_orig, y = data$y_orig), state$panel)
  xy <- cbind(base$x * mm[1], base$y * mm[2])
  anchor_mm <- cbind(anchors$x * mm[1], anchors$y * mm[2])
  baseline <- centres <- xy
  angle <- rep(0, nrow(data)); hjust <- rep(0, nrow(data)); bad <- character()
  connection <- state$options$connection_height + state$options$label_gap * max(data$size)
  label_theta <- NULL
  origin_mm <- NULL
  far <- offset_chr_points_signed(layout, as.character(data$chr),
    project_positions_checked(layout, as.character(data$chr), data$position),
    rep(spec$high_offset, nrow(data)))
  far <- coord$transform(as.data.frame(far), state$panel)
  far_mm <- cbind(far$x * mm[1], far$y * mm[2])
  if (is_circular_layout(layout)) {
    origin <- coord$transform(data.frame(x = 0, y = 0), state$panel)
    centre <- c(origin$x * mm[1], origin$y * mm[2])
    origin_mm <- centre
    radial <- xy - rep(centre, each = nrow(data))
    radius <- sqrt(rowSums(radial^2))
    sign <- if (spec$side == "right") 1 else -1
    label_radius <- radius[1] + sign * connection
    if (label_radius <= 0 || (sign < 0 && label_radius <= max(width)))
      stopf("Physical labels reach the circle center. Increase the radius or reduce text size.")
    far_radius <- sqrt(rowSums(sweep(far_mm, 2, centre, "-")^2))
    outer_corner <- sqrt((label_radius + width)^2 + (height / 2)^2)
    inner_edge <- label_radius - width
    outside <- if (sign > 0) outer_corner > far_radius else inner_edge < far_radius
    bad <- c(bad, as.character(data$label[outside]))
    label_theta <- numeric(nrow(data))
    for (chr in unique(as.character(data$chr))) {
      rows <- which(as.character(data$chr) == chr)
      g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
      packing_radius <- if (sign > 0) label_radius else label_radius - width[rows]
      half <- label_radius * atan2((height[rows] + state$options$padding) / 2, packing_radius)
      packed <- ordered_label_positions(data$x[rows] / layout$circular$radius * label_radius,
        half, g$.arc_start / layout$circular$radius * label_radius,
        g$.arc_end / layout$circular$radius * label_radius)
      if (!packed$fits) bad <- c(bad, as.character(data$label[rows]))
      theta <- layout$circular$start_angle * pi / 180 +
        (if (layout$circular$clockwise) 1 else -1) * packed$position / label_radius
      baseline[rows, ] <- cbind(centre[1] + label_radius * sin(theta),
        centre[2] + label_radius * cos(theta))
      glyph_radius <- label_radius + sign * width[rows] / 2
      centres[rows, ] <- cbind(centre[1] + glyph_radius * sin(theta),
        centre[2] + glyph_radius * cos(theta))
      label_theta[rows] <- theta
      radial_angle <- 90 - theta * 180 / pi
      angle[rows] <- (radial_angle + 90) %% 180 - 90
      hjust[rows] <- as.numeric(sign * cos((angle[rows] - radial_angle) * pi / 180) < 0)
    }
  } else {
    vertical <- layout$orientation == "vertical"
    axis <- if (vertical) 2L else 1L
    normal_axis <- 3L - axis
    sign <- if (spec$side == "right") 1 else -1
    if (!vertical) sign <- -sign
    baseline[, normal_axis] <- baseline[, normal_axis] + sign * connection
    centres[, normal_axis] <- baseline[, normal_axis] + sign * width / 2
    angle[] <- if (vertical) 0 else 90
    hjust[] <- as.numeric(sign < 0)
    available <- abs(far_mm[, normal_axis] - xy[, normal_axis])
    bad <- c(bad, as.character(data$label[available < connection + width]))
    for (chr in unique(as.character(data$chr))) {
      rows <- which(as.character(data$chr) == chr)
      g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
      ends <- coord$transform(data.frame(x = c(g$.axis_start_x, g$.axis_end_x),
        y = c(g$.axis_start_y, g$.axis_end_y)), state$panel)
      extent <- range(if (vertical) ends$y * mm[2] else ends$x * mm[1])
      half <- (height[rows] + state$options$padding) / 2
      packed <- ordered_label_positions(xy[rows, axis], half, extent[1], extent[2])
      baseline[rows, axis] <- packed$position
      centres[rows, axis] <- packed$position
      if (!packed$fits) bad <- c(bad, as.character(data$label[rows]))
    }
  }
  if (length(bad)) warning(sprintf("Labels exceed the available display span or track width: %s.",
    paste(unique(bad), collapse = ", ")), call. = FALSE)
  data$x <- baseline[, 1] / mm[1]; data$y <- baseline[, 2] / mm[2]
  data$angle <- angle; data$hjust <- hjust; data$vjust <- 0.5
  identity <- ggplot2::ggproto(NULL, ggplot2::CoordCartesian,
    transform = function(data, panel_params) data)
  text_params <- state$options$text_params
  text_params$parse <- state$options$parse
  if ("size.unit" %in% names(text_params)) text_params$size.unit <- "mm"
  text <- do.call(state$options$text_geom$draw_panel,
    c(list(data = data, panel_params = state$panel, coord = identity), text_params))
  leaders <- chr_label_leader_data(anchor_mm, baseline, layout, mm, origin_mm,
    if (is_circular_layout(layout)) circular_theta(layout, data$x_orig) else NULL,
    label_theta)
  leaders$colour <- state$options$segment_colour
  leaders$linewidth <- state$options$segment_linewidth
  leaders$linetype <- 1; leaders$alpha <- NA_real_
  segments <- ggplot2::GeomSegment$draw_panel(leaders, state$panel, identity)
  positions <- data.frame(id = as.character(data$ideogram_label_id), chr = as.character(data$chr),
    position = data$position, x = centres[, 1], y = centres[, 2],
    baseline_x = baseline[, 1], baseline_y = baseline[, 2],
    near_x = xy[, 1], near_y = xy[, 2],
    anchor_x = anchor_mm[, 1], anchor_y = anchor_mm[, 2],
    width = width, height = height, angle = angle, hjust = hjust, vjust = 0.5,
    size = data$size)
  state$device <- device
  state$cache <- list(text = text, segments = segments, leaders = leaders,
    positions = positions, fits = !length(bad))
  state$cache
}

locus_text_mm <- function(data, size.unit = "mm") {
  factor <- switch(size.unit, mm = 1, cm = 10, `in` = 25.4,
    pt = 1 / ggplot2::.pt, pc = 12 / ggplot2::.pt,
    stopf("Text `size.unit` must be mm, cm, in, pt or pc."))
  data$size <- data$size * factor
  data
}

chr_label_leader_data <- function(anchor, end, layout, mm, origin = NULL,
    theta1 = NULL, theta2 = NULL) {
  if (is_circular_layout(layout)) {
    a <- sweep(anchor, 2, origin, "-")
    b <- sweep(end, 2, origin, "-")
    delta <- theta2 - theta1
    r1 <- sqrt(rowSums(a^2)); r2 <- sqrt(rowSums(b^2))
    return(do.call(rbind, lapply(seq_len(nrow(anchor)), function(i) {
      t <- seq(0, 1, length.out = layout$curve_points)
      theta <- theta1[i] + delta[i] * t^2
      radius <- r1[i] + (r2[i] - r1[i]) * t
      x <- (origin[1] + radius * sin(theta)) / mm[1]
      y <- (origin[2] + radius * cos(theta)) / mm[2]
      data.frame(x = utils::head(x, -1L), y = utils::head(y, -1L),
        xend = utils::tail(x, -1L), yend = utils::tail(y, -1L), group = i)
    })))
  }
  data.frame(x = anchor[, 1] / mm[1], y = anchor[, 2] / mm[2],
    xend = end[, 1] / mm[1], yend = end[, 2] / mm[2], group = seq_len(nrow(anchor)))
}

#' @method makeContent chr_labels_render
#' @importFrom grid makeContent
#' @export
makeContent.chr_labels_render <- function(x) {
  solved <- resolve_chr_labels(x$state)
  x$positions <- solved$positions
  grid::setChildren(x, grid::gList(if (x$kind == "leaders") solved$segments else solved$text))
}
