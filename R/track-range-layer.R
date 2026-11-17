# Native range layers finish their statistics and Position before projection.

clamp_smooth_to_limits <- function(layout, track_id, chr, vals) {
  if (is.null(layout) || is.null(track_id) || !length(vals)) return(vals)
  spec <- track_table_row(layout, track_id)
  keys <- track_key_for_rows(spec, as.character(chr))
  out <- vals
  for (key in unique(keys)) {
    entry <- layout$track_ranges[[key]]
    if (is.null(entry) || is.null(entry$limits)) next
    idx <- which(keys == key & is.finite(out))
    if (length(idx)) {
      out[idx] <- pmin(pmax(out[idx], entry$limits[1]), entry$limits[2])
    }
  }
  out
}

native_range_track_values <- function(prototype, data, mapping, layout = NULL, track = NULL) {
  chr <- as.character(mapped_fields(data, mapping, "chr", what = "mapping")$chr)
  raw_val <- if ("value" %in% names(mapping)) {
    mapped_fields(data, mapping, "value", what = "mapping")$value
  } else if ("y" %in% names(mapping)) {
    mapped_fields(data, mapping, "y", what = "mapping")$y
  } else NULL
  is_smooth <- inherits(prototype$geom, "GeomSmooth")
  explicit_lim <- if (!is.null(layout) && !is.null(track)) {
    track_table_row(layout, track)$limits[[1]]
  } else NULL
  base_z <- if (!is.null(layout) && !is.null(track) && inherits(prototype$geom, "GeomArea")) {
    track_zero_baseline(layout, track, raw_val)
  } else 0
  mapping <- track_mapping_without(mapping, c("chr", "position", "value"))
  pieces <- lapply(unique(chr), function(key) {
    layer <- ggplot2::ggproto(NULL, prototype, data = data[chr == key, , drop = FALSE],
      mapping = mapping, inherit.aes = FALSE)
    layer$geom_params <- prototype$geom_params
    layer$stat_params <- prototype$stat_params
    layer$geom_params$orientation <- "x"
    layer$stat_params$orientation <- "x"
    computed <- suppressMessages(ggplot2::ggplot_build(ggplot2::ggplot() + layer))$data[[1]]
    values <- unlist(computed[intersect(c("y", "ymin", "ymax"), names(computed))],
      use.names = FALSE)
    values <- values[is.finite(values)]
    if (is_smooth && !is.null(explicit_lim) && length(values)) {
      raw_local <- raw_val[chr == key]
      raw_local <- raw_local[is.finite(raw_local)]
      if (!any(raw_local < explicit_lim[1] | raw_local > explicit_lim[2])) {
        values <- pmin(pmax(values, explicit_lim[1]), explicit_lim[2])
      } else {
        values <- c(raw_local, pmin(pmax(values, explicit_lim[1]), explicit_lim[2]))
      }
    }
    if (inherits(prototype$geom, "GeomArea") && is.finite(base_z) && base_z != 0 &&
        !any(raw_val[chr == key] == 0, na.rm = TRUE)) {
      values[values == 0] <- base_z
    }
    list(chr = rep(key, length(values)), value = values)
  })
  list(chr = unlist(lapply(pieces, `[[`, "chr"), use.names = FALSE),
    value = unlist(lapply(pieces, `[[`, "value"), use.names = FALSE))
}

native_range_track_constructor <- function(prototype) {
  force(prototype)
  function(mapping = NULL, data = NULL, stat = prototype$stat,
      position = prototype$position, na.rm = FALSE,
      show.legend = prototype$show.legend, inherit.aes = FALSE, ...) {
    supplied <- list(...)
    track <- supplied$track
    if (is.character(position)) {
      position <- ggplot2::layer(geom = prototype$geom, stat = "identity",
        position = position, params = list())$position
    }
    layer <- ggplot2::ggproto(NULL, prototype, mapping = mapping, data = data,
      stat = native_range_track_stat(prototype$stat, track),
      geom = native_range_track_geom(prototype$geom),
      position = native_range_track_position(position, prototype$geom),
      show.legend = show.legend, inherit.aes = inherit.aes)
    layer$geom_params <- prototype$geom_params
    layer$aes_params <- prototype$aes_params
    layer$stat_params <- prototype$stat_params
    layer$geom_params$na.rm <- na.rm
    layer$stat_params$na.rm <- na.rm
    if (inherits(prototype$stat, c("StatSmooth", "StatAlign"))) {
      layer$stat_params$orientation <- "x"
      layer$geom_params$orientation <- "x"
    }
    layer
  }
}

native_range_track_stat <- function(native, track) {
  # StatAlign / StatSmooth deliberately create new observations along x = bp.
  # Semantic aliases are reconstructed from the resulting bp/value.
  delegate <- ggplot2::ggproto(NULL, native,
    dropped_aes = union(native$dropped_aes, c("position", "value")))
  ggplot2::ggproto("StatChrTrackRange", native,
    required_aes = union(native$required_aes, "chr"),
    extra_params = union(native$extra_params,
      c("track", "chromosome_separation", "circular", "track_geometry")),
    setup_params = function(self, data, params) {
      if ("chr" %in% names(data) && nrow(data)) {
        data$group <- as.integer(interaction(as.character(data$chr), data$group,
          drop = TRUE, lex.order = TRUE))
      }
      params$orientation <- "x"
      delegate$setup_params(data, params)
    },
    compute_layer = function(self, data, params, layout) {
      if (!nrow(data)) return(data)
      chr <- as.character(data$chr)
      data$group <- as.integer(interaction(chr, data$group,
        drop = TRUE, lex.order = TRUE))
      params[c("track", "chromosome_separation", "circular", "track_geometry")] <- NULL
      pieces <- lapply(unique(chr), function(key) {
        rows <- data[chr == key, , drop = FALSE]
        computed <- delegate$compute_layer(rows, params, layout)
        if (!nrow(computed)) return(computed)
        computed$chr <- key
        computed$ideogram_chr <- key
        computed$ideogram_track <- track
        computed$ideogram_position_shift <- 0
        computed$ideogram_track_geometry <- TRUE
        if (inherits(native, c("StatAlign", "StatSmooth")) || is.null(computed$position)) {
          computed$position <- computed$x
        }
        if (inherits(native, c("StatAlign", "StatSmooth")) || is.null(computed$value)) {
          computed$value <- computed$y
        }
        computed
      })
      do.call(rbind, pieces)
    })
}

native_range_track_position <- function(native, prototype_geom = NULL) {
  ggplot2::ggproto("PositionChrTrackRange", native,
    setup_params = function(data) list(),
    setup_data = function(data, params) data,
    compute_layer = function(self, data, params, layout) {
      if (!nrow(data)) return(data)
      pieces <- lapply(unique(data$ideogram_chr), function(key) {
        rows <- data[data$ideogram_chr == key, , drop = FALSE]
        settings <- native$setup_params(rows)
        rows <- native$setup_data(rows, settings)
        native$compute_layer(rows, settings, layout)
      })
      data <- do.call(rbind, pieces)
      if (!nrow(data)) return(data)
      raw <- data
      bounds <- all(c("ymin", "ymax") %in% names(raw))
      if (bounds) {
        data$xmin <- data$x
        data$xmax <- data$x
      }
      chr_layout <- layout$coord$layout
      is_area <- inherits(prototype_geom, "GeomArea")
      is_smooth <- inherits(prototype_geom, "GeomSmooth")
      if (is_area && bounds) {
        raw$ymin <- clamp_track_baseline(chr_layout, raw$ideogram_track[1], raw$ideogram_chr, raw$ymin)
        raw$ymax <- clamp_track_baseline(chr_layout, raw$ideogram_track[1], raw$ideogram_chr, raw$ymax)
        data$ymin <- raw$ymin
        data$ymax <- raw$ymax
      }
      if (is_area && "y" %in% names(data)) {
        data$y <- clamp_track_baseline(chr_layout, data$ideogram_track[1], data$ideogram_chr, data$y)
      }
      if (is_smooth && bounds) {
        raw$ymin <- clamp_smooth_to_limits(chr_layout, raw$ideogram_track[1], raw$ideogram_chr, raw$ymin)
        raw$ymax <- clamp_smooth_to_limits(chr_layout, raw$ideogram_track[1], raw$ideogram_chr, raw$ymax)
        data$ymin <- raw$ymin
        data$ymax <- raw$ymax
      }
      if (is_smooth && "y" %in% names(data)) {
        data$y <- clamp_smooth_to_limits(chr_layout, data$ideogram_track[1], data$ideogram_chr, data$y)
      }
      data <- transform_ideogram_track(chr_layout, data)
      if (bounds) {
        lower <- upper <- raw
        lower$y <- raw$ymin
        upper$y <- raw$ymax
        lower[c("xmin", "xmax", "ymin", "ymax")] <- NULL
        upper[c("xmin", "xmax", "ymin", "ymax")] <- NULL
        lower <- transform_ideogram_track(chr_layout, lower)
        upper <- transform_ideogram_track(chr_layout, upper)
        if (chr_layout$orientation == "vertical") {
          data$xmin <- lower$x
          data$xmax <- upper$x
        } else {
          data$ymin <- lower$y
          data$ymax <- upper$y
        }
      }
      data$flipped_aes <- chr_layout$orientation == "vertical"
      data$ideogram_projected <- TRUE
      data
    })
}

native_range_track_geom <- function(native) {
  ggplot2::ggproto(NULL, native,
    parameters = function(extra = FALSE) native$parameters(extra),
    draw_group = function(self, data, panel_params, coord, ..., flipped_aes = FALSE) {
      if (inherits(coord$layout, "ideogram_layout_v2")) {
        flipped_aes <- coord$layout$orientation == "vertical"
      }
      native$draw_group(data, panel_params, coord, ..., flipped_aes = flipped_aes)
    })
}
