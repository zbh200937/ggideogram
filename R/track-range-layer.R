# Native range layers finish their statistics and Position before projection.

native_range_track_values <- function(prototype, data, mapping) {
  chr <- as.character(mapped_fields(data, mapping, "chr", what = "mapping")$chr)
  mapping <- track_mapping_without(mapping, c("chr", "position", "value"))
  pieces <- lapply(unique(chr), function(key) {
    layer <- ggplot2::ggproto(NULL, prototype, data = data[chr == key, , drop = FALSE],
      mapping = mapping, inherit.aes = FALSE)
    layer$geom_params <- prototype$geom_params
    layer$stat_params <- prototype$stat_params
    layer$geom_params$orientation <- "x"
    layer$stat_params$orientation <- "x"
    computed <- ggplot2::ggplot_build(ggplot2::ggplot() + layer)$data[[1]]
    values <- unlist(computed[intersect(c("y", "ymin", "ymax"), names(computed))],
      use.names = FALSE)
    values <- values[is.finite(values)]
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
      position = native_range_track_position(position),
      show.legend = show.legend, inherit.aes = inherit.aes)
    layer$geom_params <- prototype$geom_params
    layer$aes_params <- prototype$aes_params
    layer$stat_params <- prototype$stat_params
    layer$geom_params$na.rm <- na.rm
    layer$stat_params$na.rm <- na.rm
    layer
  }
}

native_range_track_stat <- function(native, track) {
  # StatAlign deliberately creates new observations. Semantic aliases are
  # reconstructed from its resulting bp/value rather than treated as aesthetics.
  delegate <- ggplot2::ggproto(NULL, native,
    dropped_aes = union(native$dropped_aes, c("position", "value")))
  ggplot2::ggproto("StatChrTrackRange", native,
    required_aes = union(native$required_aes, "chr"),
    extra_params = union(native$extra_params,
      c("track", "chromosome_separation", "circular", "track_geometry")),
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
        if (inherits(native, "StatAlign") || is.null(computed$position)) {
          computed$position <- computed$x
        }
        if (inherits(native, "StatAlign") || is.null(computed$value)) {
          computed$value <- computed$y
        }
        computed
      })
      do.call(rbind, pieces)
    })
}

native_range_track_position <- function(native) {
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
    draw_group = function(self, data, panel_params, coord, lineend = "butt",
        linejoin = "round", linemitre = 10, na.rm = FALSE,
        flipped_aes = FALSE, outline.type = "both") {
      if (inherits(coord$layout, "ideogram_layout_v2")) {
        flipped_aes <- coord$layout$orientation == "vertical"
      }
      native$draw_group(data, panel_params, coord, lineend = lineend,
        linejoin = linejoin, linemitre = linemitre, na.rm = na.rm,
        flipped_aes = flipped_aes, outline.type = outline.type)
    })
}
