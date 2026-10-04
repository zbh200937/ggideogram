# Value-axis components in declared chromosome tracks.

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
    size = ideogram_text_size("value"),
    colour = ideogram_colour("axis"),
    linewidth = ideogram_linewidth("axis"),
    family = NULL) {
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
  check_positive_layout(size, "size")
  check_nonnegative_layout(linewidth, "linewidth")
  structure(
    list(
      track = track, chr = chr, position = position, breaks = breaks,
      n = n, labels = labels, tick_length = tick_length,
      label_gap = label_gap, size = size, colour = colour,
      linewidth = linewidth, family = family
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
    layout$chrom$.chr
  } else {
    as.character(object$chr)
  }
  unknown <- setdiff(chromosomes, layout$chrom$.chr)
  if (length(unknown)) {
    stopf("Unknown track-axis chromosome%s: %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
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
      result$PANEL <- panel
      result$group <- -1L
      result
    }))
  }
)

build_track_axis_layers <- function(layout, object, chromosomes) {
  object$family <- object$family %||% layout$base_spec$base_family %||% ""
  spec <- track_table_row(layout, object$track)
  position <- object$position
  if (position == "auto") position <- if (is_circular_layout(layout) &&
      spec$value_scale != "per_chr" && isTRUE(object$chr)) "gap" else "end"
  in_gap <- position == "gap"
  if (in_gap && !is_circular_layout(layout)) stopf("Track-axis `position = \"gap\"` requires a circular layout.")
  if (in_gap && spec$value_scale == "per_chr" && length(chromosomes) > 1L) {
    stopf("A gap axis with `value_scale = \"per_chr\"` must select one chromosome.")
  }
  if (in_gap) chromosomes <- if (isTRUE(object$chr))
    layout$chrom$.chr[if (layout$circular$clockwise) nrow(layout$chrom) else 1L] else chromosomes[1L]
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
    breaks <- if (is.function(object$breaks)) {
      object$breaks(limits)
    } else if (is.null(object$breaks)) {
      pretty(limits, n = object$n)
    } else {
      object$breaks
    }
    if (!is.numeric(breaks) || anyNA(breaks) || any(!is.finite(breaks))) {
      stopf("Track-axis `breaks` must resolve to finite numbers.")
    }
    breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]
    if (!length(breaks)) next
    label <- if (is.function(object$labels)) {
      object$labels(breaks)
    } else {
      object$labels
    }
    if (length(label) != length(breaks)) {
      stopf("Track-axis `labels` must match the number of breaks.")
    }

    g <- layout$chrom[layout$index[[chr]], , drop = FALSE]
    genomic_position <- if (position == "end") g$.end else g$.start
    ends <- project_track_values_raw(
      layout, object$track, rep(chr, 2), rep(genomic_position, 2), limits)
    at <- project_track_values_raw(
      layout, object$track, rep(chr, length(breaks)),
      rep(genomic_position, length(breaks)), breaks)
    dx <- g$.axis_end_x - g$.axis_start_x
    dy <- g$.axis_end_y - g$.axis_start_y
    axis_length <- sqrt(dx^2 + dy^2)
    direction <- if (position == "end") 1 else -1
    if (in_gap) {
      arc <- layout$circular$label_gap$left_arc
      ends <- circular_xy(layout, rep(arc, 2), ends$y)
      at <- circular_xy(layout, rep(arc, length(breaks)), at$y)
      theta <- layout$circular$label_gap$left_angle * pi / 180
      dx <- cos(theta); dy <- -sin(theta); axis_length <- direction <- 1
      if (abs(dx) < 1e-12) dx <- 0
      if (abs(dy) < 1e-12) dy <- 0
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
    text[[index]] <- data.frame(
      x = at$x, y = at$y, label = as.character(label),
      nx = rep(direction * dx / axis_length, length(breaks)),
      ny = rep(direction * dy / axis_length, length(breaks))
    )
  }
  spine <- do.call(rbind, spine)
  ticks <- do.call(rbind, ticks)
  text <- do.call(rbind, text)
  if (is.null(spine) || !nrow(spine)) return(list())
  if (in_gap) {
    spine$ideogram_cartesian <- ticks$ideogram_cartesian <- text$ideogram_cartesian <- TRUE
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
        size = object$size, colour = object$colour,
        family = object$family, na.rm = FALSE)
    )
  )
  if (in_gap) layers[[3]]$geom_params$gap_track <- object$track
  for (part in seq_along(layers)) {
    layers[[part]]$stat <- StatTrackAxis
    layers[[part]]$stat_params <- list(
      na.rm = FALSE, axis_spec = object,
      chromosomes = chromosomes, axis_part = part)
  }
  layers
}
