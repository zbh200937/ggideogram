# Value-axis components in declared chromosome tracks.

track_axis_gap <- function(gap = NULL) {
  gap <- gap %||% 0.2
  check_nonnegative_layout(gap, "axis gap")
  gap
}

circular_gap_axis_arc <- function(layout, gap) {
  circle <- layout$circular
  room <- circle$opening_angle * pi / 180 * circle$radius
  direction <- if (circle$clockwise) 1 else -1
  circle$label_gap$left_arc + direction * min(gap, room / 2)
}

#' @noRd
geom_track_axis <- function(
    track,
    chr = TRUE,
    position = c("auto", "end", "start", "gap"),
    breaks = NULL,
    n = 1,
    labels = scales::label_number(),
    tick_length = 1.5,
    label_gap = 0.25,
    size = NULL,
    colour = ideogram_colour("axis"),
    linewidth = ideogram_linewidth("axis"),
    family = NULL, gap = NULL) {
  position <- match.arg(position)
  if (!isTRUE(chr) && (!is.character(chr) || !length(chr))) {
    stopf("`chr` must be `TRUE` or a non-empty chromosome vector.")
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 1) {
    stopf("`n` must be one positive number.")
  }
  check_nonnegative_layout(tick_length, "tick_length")
  check_nonnegative_layout(label_gap, "label_gap")
  if (!is.null(size)) check_positive_layout(size, "size")
  check_nonnegative_layout(linewidth, "linewidth")
  structure(
    list(
      track = track, chr = chr, position = position, breaks = breaks,
      n = n, labels = labels, tick_length = tick_length,
      label_gap = label_gap, size = size, colour = colour,
      linewidth = linewidth, family = family, gap = track_axis_gap(gap)
    ),
    class = "ggideogram_track_axis_component"
  )
}

#' @method ggplot_add ggideogram_track_axis_component
#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.ggideogram_track_axis_component <- function(
    object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  object_name <- ggplot_add_object_name(..., fallback = "track_axis")
  track_table_row(layout, object$track)
  chromosomes <- if (isTRUE(object$chr)) {
    track_chromosomes(layout, object$track)
  } else {
    as.character(object$chr)
  }
  unknown <- setdiff(chromosomes, layout$chrom$.chr)
  if (length(unknown)) {
    stopf("Unknown track-axis chromosome%s: %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  check_track_chromosomes(layout, object$track, chromosomes)
  if (!length(chromosomes)) return(plot)
  layers <- build_track_axis_layers(layout, object, chromosomes)
  if (length(layers) && inherits(layers[[3]]$geom, "GeomCircularGapText")) {
    layout$circular_gap_axes[[object$track]] <- object
    plot <- update_plot_ideogram_layout(plot, layout)
  }
  for (index in seq_along(layers)) {
    plot <- ggplot2::ggplot_add(
      layers[[index]], plot, paste0(object_name, "[[", index, "]]"))
  }
  reserve_chromosome_axis_space(plot, layers)
}

# Resolve axes against the final coordinate layout, after all track layers have
# registered their values. Keep the ordinary segment/tick/text Geoms intact.
StatTrackAxis <- ggplot2::ggproto(
  "StatTrackAxis", ggplot2::StatIdentity,
  extra_params = c("na.rm", "axis_spec", "chromosomes", "axis_part"),
  compute_layer = function(self, data, params, layout) {
    layers <- build_track_axis_layers(
      layout$coord$layout, params$axis_spec, params$chromosomes)
    if (!length(layers)) return(data[FALSE, , drop = FALSE])
    result <- layers[[params$axis_part]]$data
    panels <- unique(data$PANEL)
    do.call(rbind, lapply(panels, function(panel) {
      result$PANEL <- rep(panel, nrow(result))
      result$group <- rep(-1L, nrow(result))
      result
    }))
  }
)

auto_track_axis_breaks <- function(limits, n) {
  breaks <- pretty(limits, n = n)
  breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]
  if (length(breaks) >= 2L || limits[2] <= limits[1]) return(breaks)
  for (cand_n in unique(c(max(2L, ceiling(n) + 1L), 3L, 4L, 5L, 6L, 8L, 10L))) {
    cand <- pretty(limits, n = cand_n)
    cand <- cand[cand >= limits[1] & cand <= limits[2]]
    if (length(cand) >= 2L) {
      return(if (n <= 1 && length(cand) > 2L) c(cand[1L], cand[length(cand)]) else cand)
    }
  }
  limits
}

build_track_axis_layers <- function(layout, object, chromosomes) {
  object$family <- object$family %||% layout$base_spec$base_family %||% ""
  spec <- track_table_row(layout, object$track)
  raw_position <- object$position
  position <- raw_position
  gap <- track_axis_gap(object$gap)
  if (position == "auto") position <- if (is_circular_layout(layout) &&
      spec$value_scale != "per_chr" && isTRUE(object$chr)) "gap" else "end"
  in_gap <- position == "gap"
  if (in_gap && !is_circular_layout(layout)) stopf("Track-axis `position = \"gap\"` requires a circular layout.")
  if (in_gap && spec$value_scale == "per_chr" && length(chromosomes) > 1L) {
    stopf("A gap axis with `value_scale = \"per_chr\"` must select one chromosome.")
  }
  if (in_gap) chromosomes <- if (isTRUE(object$chr) && layout$circular$clockwise)
    utils::tail(chromosomes, 1L) else chromosomes[1L]
  spine <- list()
  ticks <- list()
  text <- list()
  for (index in seq_along(chromosomes)) {
    chr <- chromosomes[index]
    key <- track_key_for_rows(spec, chr)
    limits <- layout$track_ranges[[key]]$limits
    if (is.null(limits)) {
      stopf(paste0(
        "Track `%s` has no resolved range for chromosome `%s`.\n",
        "  Add its data layer before an automatic-range axis, or declare ",
        "explicit track limits."), object$track, chr)
    }
    if (is.null(object$breaks)) {
      breaks <- auto_track_axis_breaks(limits, object$n)
      fixed_labels <- NULL
    } else {
      raw_breaks <- if (is.function(object$breaks)) object$breaks(limits) else object$breaks
      if (!is.numeric(raw_breaks) || anyNA(raw_breaks) || any(!is.finite(raw_breaks))) {
        stopf("Track-axis `breaks` must resolve to finite numbers.")
      }
      keep <- raw_breaks >= limits[1] & raw_breaks <= limits[2]
      if (!is.function(object$labels) && !is.null(object$labels)) {
        if (length(object$labels) != length(raw_breaks)) {
          stopf("Track-axis `labels` must match the number of breaks.")
        }
        fixed_labels <- object$labels[keep]
      } else {
        fixed_labels <- NULL
      }
      breaks <- raw_breaks[keep]
    }
    if (!length(breaks)) next
    label <- if (!is.null(fixed_labels)) {
      fixed_labels
    } else if (is.function(object$labels)) {
      object$labels(breaks)
    } else {
      NULL
    }
    if (is.null(label)) label <- rep("", length(breaks))
    if (length(label) != length(breaks)) {
      stopf("Track-axis `labels` must match the number of breaks.")
    }

    g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
    eff_pos <- if (raw_position == "auto" && !in_gap && !is_circular_layout(layout) &&
        g$.direction < 0) "start" else position
    genomic_position <- if (eff_pos == "end") g$.end else g$.start
    ends <- project_track_values_raw(
      layout, object$track, rep(chr, 2), rep(genomic_position, 2), limits)
    at <- project_track_values_raw(
      layout, object$track, rep(chr, length(breaks)),
      rep(genomic_position, length(breaks)), breaks)
    dx <- g$.axis_end_x - g$.axis_start_x
    dy <- g$.axis_end_y - g$.axis_start_y
    axis_length <- sqrt(dx^2 + dy^2)
    direction <- if (eff_pos == "end") 1 else -1
    if (in_gap) {
      arc <- circular_gap_axis_arc(layout, gap)
      ends <- circular_xy(layout, rep(arc, 2), ends$y)
      at <- circular_xy(layout, rep(arc, length(breaks)), at$y)
      theta <- circular_theta(layout, arc)
      dx <- cos(theta); dy <- -sin(theta); axis_length <- direction <- 1
      if (abs(dx) < 1e-12) dx <- 0
      if (abs(dy) < 1e-12) dy <- 0
    } else {
      if (is_circular_layout(layout)) {
        # Source start/end can reverse; choose the sector gap on the display side.
        after <- direction * dx > 0
        keys <- names(layout$circular$gap_angle)
        i <- match(chr, keys)
        adjacent <- if (after) i else if (i == 1L) length(keys) else i - 1L
        room <- layout$circular$gap_angle[adjacent] * pi / 180 * layout$circular$radius
        end_gap <- min(gap, room / 2)
      } else end_gap <- gap
      sx <- direction * dx / axis_length * end_gap
      sy <- direction * dy / axis_length * end_gap
      ends$x <- ends$x + sx; ends$y <- ends$y + sy
      at$x <- at$x + sx; at$y <- at$y + sy
    }
    spine[[index]] <- data.frame(
      x = ends$x[1], y = ends$y[1],
      xend = ends$x[2], yend = ends$y[2]
    )
    ticks[[index]] <- data.frame(
      x = at$x, y = at$y,
      nx = rep(direction * dx / axis_length, length(breaks)),
      ny = rep(direction * dy / axis_length, length(breaks))
    )
    # An empty label suppresses only the text, leaving its value tick intact.
    label <- as.character(label)
    show <- !is.na(label) & nzchar(trimws(label))
    text[[index]] <- cbind(ticks[[index]][show, , drop = FALSE], label = label[show])
    if (is.null(object$size)) text[[index]]$track_axis_index <- rep(index, sum(show))
  }
  spine <- do.call(rbind, spine)
  ticks <- do.call(rbind, ticks)
  text <- do.call(rbind, text)
  if (is.null(spine) || !nrow(spine)) return(list())
  if (in_gap) {
    spine$ideogram_cartesian <- ticks$ideogram_cartesian <- TRUE
    text$ideogram_cartesian <- rep(TRUE, nrow(text))
  }
  layers <- list(
    ggplot2::geom_segment(
      data = spine,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      colour = object$colour, linewidth = object$linewidth,
      inherit.aes = FALSE, show.legend = FALSE
    ),
    ggplot2::layer(
      data = ticks,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity", geom = GeomIdeogramTick,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = object$tick_length, colour = object$colour,
        linewidth = object$linewidth, na.rm = FALSE)
    ),
    ggplot2::layer(
      data = text,
      mapping = ggplot2::aes(
        x = .data$x, y = .data$y, label = .data$label,
        nx = .data$nx, ny = .data$ny),
      stat = "identity", position = "identity",
      geom = if (in_gap) GeomCircularGapText else GeomIdeogramAxisText,
      inherit.aes = FALSE, show.legend = FALSE,
      params = list(
        tick_length = object$tick_length,
        label_gap = object$label_gap,
        size = object$size %||% ideogram_text_size("value"), colour = object$colour,
        family = object$family, na.rm = FALSE)
    )
  )
  if (is.null(object$size)) layers[[3]]$geom_params$track_axis_width <- spec$width
  if (in_gap) layers[[3]]$geom_params$gap_track <- object$track
  for (part in seq_along(layers)) {
    layers[[part]]$stat <- StatTrackAxis
    layers[[part]]$stat_params <- list(
      na.rm = FALSE, axis_spec = object,
      chromosomes = chromosomes, axis_part = part)
  }
  layers
}

adaptive_track_axis_grob <- function(text, coordinates, coord, panel, width, label_gap) {
  corners <- coord$transform(data.frame(x = c(0, width), y = c(0, width),
    ideogram_cartesian = TRUE), panel)
  grid::gTree(text = text, coordinates = coordinates, label_gap = label_gap,
    span = c(diff(corners$x), diff(corners$y)), cl = "ideogram_track_axis_text")
}

#' @method makeContent ideogram_track_axis_text
#' @export
makeContent.ideogram_track_axis_text <- function(x) {
  text <- x$text
  data <- x$coordinates
  panel_mm <- c(grid::convertWidth(grid::unit(1, "npc"), "mm", TRUE),
    grid::convertHeight(grid::unit(1, "npc"), "mm", TRUE))
  glyphs <- chromosome_text_dimensions(data)
  width <- grid::convertWidth(glyphs$width, "mm", TRUE)
  height <- grid::convertHeight(glyphs$height, "mm", TRUE)
  min_size <- 1.8
  size <- rep(min_size, nrow(data))
  for (rows in split(seq_len(nrow(data)), data$track_axis_index)) {
    first <- rows[1]
    # The value spine runs perpendicular to the outward tick direction.
    spine <- c(-data$ny[first], data$nx[first]) * x$span * panel_mm
    length_mm <- sqrt(sum(spine^2))
    if (length_mm == 0) next
    limit <- length_mm / 3
    if (length(rows) > 1L) {
      along <- spine / length_mm
      position <- data$x[rows] * panel_mm[1] * along[1] +
        data$y[rows] * panel_mm[2] * along[2]
      extent <- width[rows] * abs(along[1]) + height[rows] * abs(along[2])
      extent <- extent[order(position)] / data$size[first]
      pitch <- diff(sort(position))
      pair_extent <- (utils::head(extent, -1) + utils::tail(extent, -1)) / 2
      limit <- min(limit, 0.8 * pitch / pmax(pair_extent, .Machine$double.eps))
    }
    size[rows] <- max(min_size, min(data$size[first], limit))
  }
  delta <- grid::unit(x$label_gap * (size - data$size), "mm")
  text$gp$fontsize <- size * ggplot2::.pt
  text$x <- text$x + data$nx * delta
  text$y <- text$y + data$ny * delta
  grid::setChildren(x, grid::gList(text))
}
