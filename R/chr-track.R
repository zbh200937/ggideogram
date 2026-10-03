# Track scopes keep ordinary layers and their data together before projection.

ggideogram_track_spec <- function(...) track(...)

#' Draw ordinary ggplot layers in a named chromosome track
#'
#' `geom_track()` declares a track and its contents together. Add it with `+`
#' and a unique `track` identifier, or put it in a named `ggideogram(tracks =
#' list(...))`. The same function handles left/right and inner/outer tracks.
#' Within a track, `x` is original genomic bp and `y` is the measured value,
#' including in vertical and circular layouts. Ordinary layers inherit the
#' track's data and mapping; their own data and mappings can override them.
#'
#' Pass an ordinary layer to `geom`, or several layers to `layers`. Compatible
#' identity-stat geoms keep their original Geom, Position, parameters and keys.
#' Area, ribbon, boxplot and violin layers finish their native statistics and
#' Position in raw bp/value coordinates before projection. Widths use bp;
#' violin bandwidths use the original value units. Automatic value limits
#' include statistical geometry such as density tails and boxplot notches.
#' All contents train the shared value range before drawing.
#' Adding another track recomputes existing components from source data.
#' Reuse a track identifier without `geom`/`layers` to adjust it. Omitted
#' settings remain unchanged; `axis = FALSE` closes its axis and `label = NULL`
#' removes its title. Lists update title or axis settings individually.
#'
#' @param mapping A mapping such as `aes(chr = Chr, x = Pos, y = Value)`.
#' @param data Track data. Full source data can be supplied to a local view.
#' @param geom An ordinary layer or geom constructor. `NULL` for `layers` or
#'   an empty track declaration.
#' @param track A unique track name when adding with `+`. In a named list the
#'   list name supplies the identifier.
#' @param side Track side: right/left, outer/inner for circles, or overlay
#'   inside the chromosome body.
#' @param width,gap Track width and gap in chromosome-body-width units. The
#'   default width is 2 beside the body and 1 inside it; gap is 0.2 beside
#'   and 0 inside. Native geom widths
#'   belong to the layer, e.g. geom = ggplot2::geom_col(width = 1e6).
#' @param value_scale Auto range: per_track (default), per_chr, or global.
#' @param limits Optional increasing raw-value limits.
#' @param transform A scales transformation name or object.
#' @param reverse Reverse the transverse value direction.
#' @param offset Signed overlay displacement; the full width must fit in the
#'   chromosome body. Beside tracks use zero.
#' @param layers A layer or list of layers, including [geom_genemodel()] and
#'   [geom_locus()]. Nested domain layers inherit this track identifier.
#' @param replace Replace this track's contents instead of appending layers.
#'   Geometry and omitted track settings remain unchanged. Use with new
#'   native layers to edit their point, line, column or text styling.
#' @param label Optional title string or a list with text, size (mm), colour,
#'   family, fontface and gap (em). A list can adjust the existing title's
#'   style without repeating its text. Titles sit beside the start endpoint
#'   in linear layouts and inside the closing gap in circular layouts.
#' @param axis `FALSE`, `TRUE`, or a list of local value-axis arguments.
#'   `position` accepts `"auto"` (default), `"start"`, `"end"`, or `"gap"`.
#'   Auto places shared circular value axes in the closing gap; linear and
#'   explicitly chromosome-specific axes keep their end placement. A gap
#'   axis with per-chromosome ranges must select one chromosome with `chr`.
#' @param chr Optional chromosome identifiers to select from the source data.
#' @param clip Tile clipping: `"auto"` clips overlay tiles to the chromosome
#'   silhouette, `"on"` requires an overlay track, `"off"` leaves tiles intact.
#' @param ... Parameters for `geom` when it is a constructor function.
#' @return A named-track component accepted by `ggideogram()` and `+`.
#' @examples
#' k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
#' d <- data.frame(Chr = rep(c("A", "B"), each = 3),
#'   Pos = rep(c(20, 50, 80), 2), Value = c(2, 5, 3, 4, 2, 6))
#' ggideogram(k, orientation = "horizontal", axis = TRUE) +
#'   geom_track(data = d, ggplot2::aes(chr = Chr, x = Pos, y = Value),
#'     track = "signal", side = "right", label = "Signal", axis = TRUE,
#'     layers = list(ggplot2::geom_line(), ggplot2::geom_point(size = 1)))
#' @export
geom_track <- function(mapping = NULL, data = NULL, geom = NULL, track = NULL,
    side = c("right", "left", "overlay", "outer", "inner"), width = 2,
    gap = NULL, value_scale = c("per_track", "global", "per_chr"),
    limits = NULL, transform = "identity", reverse = FALSE, offset = 0,
    layers = NULL, label = NULL, axis = FALSE, chr = NULL,
    clip = c("auto", "on", "off"), replace = FALSE, ...) {
  side <- match.arg(side)
  if (missing(width) && side == "overlay") width <- 1
  spec <- ggideogram_track_spec(side, width, gap, value_scale, limits,
    transform, reverse, offset, .partial = !is.null(track))
  if (!is.null(geom)) {
    if (is.function(geom)) geom <- geom(...)
    else check_unused_args(...)
    layers <- c(list(geom), chr_track_layers(layers))
  } else {
    check_unused_args(...)
    layers <- chr_track_layers(layers)
  }
  if (!is.null(data) && !is.data.frame(data)) stopf("Track `data` must be a data frame.")
  if (!is.null(mapping) && !inherits(mapping, "uneval") && !inherits(mapping, "ggplot2::mapping")) {
    stopf("Track `mapping` must be created with aes().")
  }
  label <- check_track_label(label)
  if ((!is.logical(axis) || length(axis) != 1L || is.na(axis)) && !is.list(axis)) {
    stopf("Track `axis` must be FALSE, TRUE, or an argument list.")
  }
  validate_optional_marker_track(track)
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) stopf("`replace` must be TRUE or FALSE.")
  spec$data <- data; spec$mapping <- mapping; spec$layers <- layers
  spec$label <- label; spec$axis <- axis; spec$chr <- chr; spec$track <- track
  spec$clip <- match.arg(clip)
  spec$replace <- replace
  spec$supplied_args <- names(match.call())[-1]
  spec$geometry_args <- intersect(spec$supplied_args,
    c("side", "width", "gap", "value_scale", "limits", "transform", "reverse", "offset"))
  class(spec) <- c("ggideogram_track_scope", class(spec))
  spec
}

check_track_label <- function(label) {
  if (is.null(label)) return(NULL)
  if (is.character(label) && length(label) == 1L && !is.na(label)) return(label)
  allowed <- c("text", "size", "colour", "family", "fontface", "gap")
  if (!is.list(label) || (length(label) && (is.null(names(label)) ||
      any(!nzchar(names(label))) || any(!names(label) %in% c(allowed, "color"))))) {
    stopf("Track `label` must be NULL, one string, or a list of text/size/colour/family/fontface/gap.")
  }
  if ("color" %in% names(label)) {
    label$colour <- label$colour %||% label$color
    label$color <- NULL
  }
  if (!is.null(label$text) && (!is.character(label$text) ||
      length(label$text) != 1L || is.na(label$text))) stopf("Track title `text` must be one string.")
  if (!is.null(label$size)) check_positive_layout(label$size, "label size")
  if (!is.null(label$gap)) check_nonnegative_layout(label$gap, "label gap")
  label
}

track_label_settings <- function(label) {
  if (is.character(label)) list(text = label) else label %||% list()
}

chr_track_layers <- function(layers) {
  if (is.null(layers)) return(list())
  if (inherits(layers, "ggplot") || inherits(layers, "grob")) {
    stopf("Track contents must be layers; use geom_locus_inset() for a complete plot or grob.")
  }
  if (inherits(layers, "Layer") || inherits(layers, "ggideogram_object") ||
      any(grepl("^ggideogram_.*component$", class(layers)))) return(list(layers))
  if (!is.list(layers)) stopf("Track `layers` must be a layer or a list of layers.")
  unlist(lapply(layers, chr_track_layers), recursive = FALSE)
}

normalise_chr_tracks <- function(tracks) {
  if (is.null(tracks)) return(NULL)
  if (!inherits(tracks, "ideogram_track_layout")) {
    if (!is.list(tracks)) stopf("`tracks` must be a named list of geom_track() declarations.")
    tracks <- do.call(track_layout, tracks)
  }
  for (id in names(tracks)) {
    specified <- tracks[[id]]$track
    if (!is.null(specified) && specified != id) stopf("Track `%s` has a different identifier `%s`.", id, specified)
    tracks[[id]] <- validate_track_spec(tracks[[id]])
  }
  tracks
}

scope_mapping <- function(parent, own, inherit = TRUE) {
  if (!inherit) return(own)
  if (is.null(own)) return(parent)
  if (is.null(parent)) return(own)
  combine_track_mapping(parent, own)
}

scope_layer_data <- function(own, parent) {
  if (is.null(own) || inherits(own, "waiver")) return(parent)
  if (is.function(own)) return(own(parent))
  own
}

scope_track_mapping <- function(mapping) {
  if (is.null(mapping)) return(mapping)
  if (!"position" %in% names(mapping) && "x" %in% names(mapping)) mapping$position <- mapping$x
  if (!"value" %in% names(mapping) && "y" %in% names(mapping)) mapping$value <- mapping$y
  track_mapping_without(mapping, c("x", "y"))
}

scope_point_data <- function(data, mapping, layout, selected = NULL) {
  if (!is.data.frame(data)) stopf("Track layers require data-frame `data`.")
  fields <- mapped_fields(data, mapping, c("chr", "position"), what = "mapping")
  chr <- validate_chr(fields$chr, "chr")
  bp <- validate_coordinate(fields$position, "position")
  source <- layout$data$karyotype
  i <- match(chr, source$.chr)
  if (anyNA(i) || any(bp < source$.start[i] | bp > source$.end[i])) {
    stopf("Track positions are outside the source chromosome bounds.")
  }
  i <- match(chr, layout$chrom$.chr)
  keep <- !is.na(i) & bp >= layout$chrom$.start[i] & bp <= layout$chrom$.end[i]
  if (!is.null(selected)) keep <- keep & chr %in% selected
  data[keep, , drop = FALSE]
}

track_layer_constructor <- function(prototype, ribbon = FALSE) {
  force(prototype); force(ribbon)
  function(mapping = NULL, data = NULL, stat = prototype$stat,
      position = prototype$position, na.rm = FALSE, show.legend = prototype$show.legend,
      inherit.aes = FALSE, ...) {
    supplied <- list(...)
    if (is.character(stat)) {
      if (!identical(stat, "identity")) stopf("Unsupported projected track Stat `%s`.", stat)
      stat <- ggplot2::StatIdentity
    }
    layer <- ggplot2::ggproto(NULL, prototype, data = data, mapping = mapping,
      stat = stat, position = position, show.legend = show.legend, inherit.aes = inherit.aes)
    if (ribbon) layer$geom <- ggplot2::GeomRibbon
    layer$geom_params <- prototype$geom_params
    layer$aes_params <- prototype$aes_params
    layer$stat_params <- prototype$stat_params
    if (identical(stat, StatChrTrack)) {
      layer$stat_params <- supplied[intersect(names(supplied),
        c("track", "chromosome_separation", "circular", "track_geometry"))]
    }
    layer$stat_params$na.rm <- na.rm
    layer$geom_params$na.rm <- na.rm
    if (!is.null(supplied$orientation)) {
      layer$geom_params$orientation <- supplied$orientation
      if ("orientation" %in% prototype$stat$parameters()) layer$stat_params$orientation <- supplied$orientation
    }
    for (aesthetic in intersect(c("xend", "xmin", "xmax"), names(mapping))) {
      if (!aesthetic %in% names(supplied)) layer$aes_params[[aesthetic]] <- NULL
    }
    layer
  }
}

native_track_component <- function(layer, spec, id, layout) {
  domain <- inherits(layer$stat, "StatChrMarker") || inherits(layer$stat, "StatChrInterval")
  mapping <- scope_mapping(spec$mapping, layer$mapping, domain || isTRUE(layer$inherit.aes))
  data <- scope_layer_data(layer$data, spec$data)
  if (domain) {
    prototype <- layer
    layer <- ggplot2::ggproto(NULL, prototype, data = data, mapping = mapping, inherit.aes = FALSE)
    layer$stat_params <- utils::modifyList(layer$stat_params, list(track = id))
    if (inherits(layer$stat, "StatChrInterval")) layer$stat_params$placement <- "beside"
    if (!is.null(layer$ideogram_marker_spec)) layer$ideogram_marker_spec$track <- id
    if ("position" %in% names(mapping)) layer$data <- scope_point_data(data, mapping, layout, spec$chr)
    if (inherits(layer$stat, "StatChrInterval")) {
      f <- mapped_fields(data, mapping, c("chr", "start", "end"), what = "mapping")
      chr <- validate_chr(f$chr, "chr")
      first <- validate_coordinate(f$start, "start")
      last <- validate_coordinate(f$end, "end")
      source <- layout$data$karyotype
      i <- match(chr, source$.chr)
      if (anyNA(i) || any(first < 1 | first > last | first != floor(first) |
          last != floor(last) | first - 1 < source$.start[i] | last > source$.end[i])) {
        stopf("Track intervals are outside the source chromosome bounds.")
      }
      i <- match(chr, layout$chrom$.chr)
      keep <- !is.na(i) & last > layout$chrom$.start[i] & first - 1 < layout$chrom$.end[i]
      if (!is.null(spec$chr)) keep <- keep & chr %in% spec$chr
      data$.source_start <- data$.source_start %||% first
      data$.source_end <- data$.source_end %||% last
      data$.ideogram_start <- pmax(first, floor(layout$chrom$.start[i]) + 1)
      data$.ideogram_end <- pmin(last, ceiling(layout$chrom$.end[i]))
      layer$data <- data[keep, , drop = FALSE]
      layer$mapping$start <- ggplot2::aes(start = .data$.ideogram_start)$start
      layer$mapping$end <- ggplot2::aes(end = .data$.ideogram_end)$end
    }
    return(layer)
  }
  kind <- if (inherits(layer$geom, "GeomArea") &&
      (inherits(layer$stat, "StatAlign") || inherits(layer$stat, "StatIdentity"))) "area" else
    if (inherits(layer$geom, "GeomRibbon") && inherits(layer$stat, "StatIdentity")) "ribbon" else
    if (inherits(layer$geom, "GeomBoxplot") && inherits(layer$stat, "StatBoxplot")) "boxplot" else
    if (inherits(layer$geom, "GeomViolin") && inherits(layer$stat, "StatYdensity")) "violin" else
    if (inherits(layer$stat, "StatIdentity")) "dynamic" else NA_character_
  if (is.na(kind)) stopf("Track `%s` does not support the `%s` Stat; supply prepared values with an identity layer.",
    id, class(layer$stat)[1])
  range_layer <- kind %in% c("area", "ribbon")
  area <- identical(kind, "area")
  mapping <- scope_track_mapping(mapping)
  if (range_layer) {
    mapping$x <- mapping$position
    if (is.null(mapping$value)) {
      lower <- mapping$ymin %||% layer$aes_params$ymin
      upper <- mapping$ymax %||% layer$aes_params$ymax
      if (is.null(lower) || is.null(upper)) stopf("A ribbon track requires ymin and ymax.")
      mapping$value <- rlang::new_quosure(rlang::expr((!!lower + !!upper) / 2))
    }
    mapping$y <- mapping$value
    kind <- "dynamic"
  }
  data <- scope_point_data(data, mapping, layout, spec$chr)
  if (inherits(layer$geom, "GeomTile")) {
    clip <- spec$clip %||% "auto"
    if (clip == "on" && spec$side != "overlay") stopf("Tile clipping requires an overlay track.")
    if (clip == "on" || (clip == "auto" && spec$side == "overlay")) layer <- clipped_native_tile(layer)
  }
  params <- utils::modifyList(layer$geom_params, layer$aes_params)
  params[c("orientation", "na.rm", "position", "show.legend", "inherit.aes")] <- NULL
  constructor <- if (range_layer) native_range_track_constructor(layer) else track_layer_constructor(layer)
  object <- new_track_component(kind, mapping, data, constructor,
    id, layer$position, params, layer$geom_params$na.rm %||% FALSE,
    layer$show.legend, FALSE, include_zero = inherits(layer$geom, "GeomCol") || area)
  object$native_geom <- layer$geom
  if (range_layer) object$range_prototype <- layer
  object
}

clipped_native_tile <- function(layer) {
  parent <- layer$geom
  prototype <- layer
  layer <- ggplot2::ggproto(NULL, prototype)
  layer$geom <- ggplot2::ggproto("GeomChrTrackTile", parent,
    parameters = function(extra = FALSE) parent$parameters(extra),
    draw_panel = function(data, panel_params, coord, ...) {
      params <- list(...)
      mask_chr_geometry(data, panel_params, coord, function(rows) {
        do.call(parent$draw_panel, c(list(data = rows, panel_params = panel_params,
          coord = coord), params))
      })
    })
  layer
}

prepare_track_scope <- function(spec, id, layout) {
  parts <- spec$parts %||% list(spec)
  unlist(lapply(parts, function(part) {
    if (is.null(part$inherit_data) || isTRUE(part$inherit_data)) part$data <- spec$data
    if (is.null(part$inherit_mapping) || isTRUE(part$inherit_mapping)) part$mapping <- spec$mapping
    part$chr <- spec$chr
    part$clip <- spec$clip
    lapply(part$layers %||% list(), function(layer) {
      if (inherits(layer, "ggideogram_object")) layer <- layer$component
      if (inherits(layer, "ggideogram_chr_component") && any(layer$components != "fill")) {
        stopf("Add chromosome bodies, names and bp axes to the main plot with geom_chr().")
      }
      if (inherits(layer, "ggideogram_pair_component")) {
        stopf("Add double-ended chromosome relations to the main plot with geom_chrlink().")
      }
      if (inherits(layer, "Layer")) return(native_track_component(layer, part, id, layout))
      if (!is.list(layer) || !any(grepl("^ggideogram_.*component$", class(layer)))) {
        stopf("Track `%s` contains an unsupported object; supply an ordinary geom layer.", id)
      }
      if (!is.null(layer$track) && (!is.character(layer$track) || !identical(layer$track, id))) {
        stopf("A layer inside track `%s` targets another track.", id)
      }
      layer$track <- id
      layer$data <- scope_layer_data(layer$data, part$data)
      layer$mapping <- scope_mapping(part$mapping, layer$mapping)
      if (inherits(layer, "ggideogram_track_component")) layer$mapping <- scope_track_mapping(layer$mapping)
      if (!is.null(part$chr)) {
        interval <- all(c("start", "end") %in% names(layer$mapping))
        f <- mapped_fields(layer$data, layer$mapping,
          if (interval) c("chr", "start", "end") else "chr", what = "mapping")
        if (interval) {
          first <- validate_coordinate(f$start, "start")
          last <- validate_coordinate(f$end, "end")
          if (any(first < 1 | first > last | first != floor(first) | last != floor(last))) {
            stopf("Track annotations require positive integer start <= end.")
          }
          source <- layout$data$karyotype
          index <- match(validate_chr(f$chr, "chr"), source$.chr)
          if (anyNA(index) || any(first - 1 < source$.start[index] | last > source$.end[index])) {
            stopf("Track annotations are outside the source chromosome bounds.")
          }
        } else if ("position" %in% names(layer$mapping)) {
          scope_point_data(layer$data, layer$mapping, layout, part$chr)
        }
        layer$data <- layer$data[as.character(f$chr) %in% part$chr, , drop = FALSE]
      }
      if (inherits(layer, "ggideogram_track_component")) {
        layer$data <- scope_point_data(layer$data, layer$mapping, layout)
      }
      layer
    })
  }), recursive = FALSE)
}

train_track_scopes <- function(plot, components) {
  layout <- ideogram_plot_layout(plot)
  frozen <- character()
  for (spec in layout$base_spec$tracks %||% list()) {
    for (key in names(spec$frozen_ranges)) layout$track_ranges[[key]] <- spec$frozen_ranges[[key]]
  }
  for (object in components) {
    if (!inherits(object, "ggideogram_track_component") || !nrow(object$data)) next
    if (object$kind == "dynamic") {
      layout <- register_dynamic_track_range(object, layout)$layout
    } else {
      fields <- mapped_fields(object$data, object$mapping,
        if (object$kind == "ribbon") c("chr", "position", "ymin", "ymax") else c("chr", "position", "value"), what = "mapping")
      values <- if (object$kind == "ribbon") c(fields$ymin, fields$ymax) else fields$value
      chr <- rep(as.character(fields$chr), length.out = length(values))
      if (object$kind == "area") { values <- c(values, 0); chr <- c(chr, chr[1]) }
      if (object$kind == "area" && !inherits(object$position, "PositionIdentity")) {
        stacked <- track_position_value_range(object$position, fields$chr, fields$position, fields$value)
        if (!is.null(stacked)) { values <- c(values, stacked$value); chr <- c(chr, stacked$chr) }
      }
      if (object$kind %in% c("boxplot", "violin")) {
        computed <- distribution_track_values(object, fields)
        values <- c(values, computed$value)
        chr <- c(chr, computed$chr)
      }
      layout <- register_track_values(layout, object$track, chr, values)
      frozen <- union(frozen, object$track)
    }
  }
  for (id in frozen) {
    keys <- track_range_keys(track_table_row(layout, id), layout$chrom$.chr)
    for (key in keys) if (!is.null(layout$track_ranges[[key]]$limits)) layout$track_ranges[[key]]$fixed <- TRUE
    layout$base_spec$tracks[[id]]$frozen_ranges <- layout$track_ranges[intersect(keys, names(layout$track_ranges))]
  }
  layout$scope_fixed_tracks <- union(layout$scope_fixed_tracks %||% character(), frozen)
  update_plot_ideogram_layout(plot, layout)
}

track_scope_owners <- function(layout, spec, contents) {
  owners <- unique(unlist(lapply(contents, function(object) {
    if (is.null(object$data) || is.null(object$mapping) || !nrow(object$data)) return(NULL)
    as.character(mapped_fields(object$data, object$mapping, "chr", what = "mapping")$chr)
  }), use.names = FALSE))
  if (!length(owners)) owners <- spec$chr %||% layout$chrom$.chr
  layout$chrom$.chr[layout$chrom$.chr %in% owners]
}

add_track_scope_title <- function(plot, id, title, owners) {
  settings <- utils::modifyList(list(size = ideogram_text_size("track"),
    colour = "#202020", family = "", fontface = 1, gap = 0.3),
    track_label_settings(title))
  if (is.null(settings$text) || !nzchar(settings$text) || !length(owners)) return(plot)
  layout <- ideogram_plot_layout(plot)
  selected <- match(owners[1], layout$chrom$.chr)
  g <- layout$chrom[selected, , drop = FALSE]
  spec <- track_table_row(layout, id)
  reverse <- g$.direction < 0
  bp <- if (reverse) g$.end else g$.start
  p <- offset_chr_points_signed(layout, g$.chr,
    project_positions_checked(layout, g$.chr, bp), mean(c(spec$low_offset, spec$high_offset)))
  data <- data.frame(x = p$x, y = p$y, label = settings$text,
    nx = if (layout$orientation == "vertical") 0 else -1,
    ny = if (layout$orientation == "vertical") 1 else 0)
  if (is_circular_layout(layout)) {
    point <- circular_xy(layout, layout$circular$label_gap$left_arc, p$y)
    theta <- layout$circular$label_gap$left_angle * pi / 180
    data$x <- point$x; data$y <- point$y
    data$nx <- cos(theta); data$ny <- -sin(theta)
    if (abs(data$nx) < 1e-12) data$nx <- 0
    if (abs(data$ny) < 1e-12) data$ny <- 0
    data$ideogram_cartesian <- TRUE
  }
  layer <- ggplot2::layer(data = data, mapping = ggplot2::aes(
    x = .data$x, y = .data$y, label = .data$label, nx = .data$nx, ny = .data$ny),
    stat = "identity", geom = if (is_circular_layout(layout)) GeomCircularGapText else
      GeomIdeogramAxisText, position = "identity",
    inherit.aes = FALSE, show.legend = FALSE,
    params = list(tick_length = 1, label_gap = settings$gap,
      size = settings$size, colour = settings$colour, family = settings$family,
      fontface = settings$fontface, na.rm = FALSE))
  layer$ideogram_track_title <- list(chr = g$.chr)
  if (is_circular_layout(layout)) {
    layer$geom_params$gap_track <- id
    layer$geom_params$gap_title <- TRUE
    layout$circular_gap_titles[[id]] <- list(data = data,
      size = settings$size, family = settings$family, tick_length = 1, label_gap = settings$gap)
    plot <- update_plot_ideogram_layout(plot, layout)
  }
  reserve_chromosome_axis_space(plot + layer, list(layer))
}

add_chr_track_contents <- function(plot, tracks) {
  if (is.null(tracks)) return(plot)
  components <- lapply(names(tracks), function(id) prepare_track_scope(tracks[[id]], id, ideogram_plot_layout(plot)))
  plot <- train_track_scopes(plot, unlist(components, recursive = FALSE))
  for (i in seq_along(tracks)) {
    id <- names(tracks)[i]; start <- length(plot$layers)
    for (content in seq_along(components[[i]])) {
      component <- components[[i]][[content]]
      if (is.data.frame(component$data) && !nrow(component$data)) next
      before <- length(plot$layers)
      plot <- ggplot2::ggplot_add(component, plot, id)
      if (length(plot$layers) > before) for (j in seq.int(before + 1L, length(plot$layers))) {
        plot$layers[[j]]$ideogram_scope_slot <- paste("content", content, j - before, sep = ":")
      }
    }
    layout <- ideogram_plot_layout(plot)
    owners <- track_scope_owners(layout, tracks[[i]], components[[i]])
    before <- length(plot$layers)
    plot <- add_track_scope_title(plot, id, tracks[[i]]$label, owners)
    if (length(plot$layers) > before) plot$layers[[length(plot$layers)]]$ideogram_scope_slot <- "title"
    axis <- tracks[[i]]$axis
    if (isTRUE(axis) || is.list(axis)) {
      args <- if (is.list(axis)) axis else list()
      if (is.null(args$chr)) {
        if (is_circular_layout(layout) && tracks[[i]]$value_scale != "per_chr" &&
            (args$position %||% "auto") == "auto") args$position <- "gap"
        args$chr <- if (identical(args$position, "gap")) owners[1] else owners
      }
      before <- length(plot$layers)
      if (length(owners)) plot <- plot + do.call(geom_track_axis, c(list(track = id), args))
      if (length(plot$layers) > before) for (j in seq.int(before + 1L, length(plot$layers))) {
        plot$layers[[j]]$ideogram_scope_slot <- paste0("axis:", j - before)
      }
    }
    if (length(plot$layers) > start) for (j in seq.int(start + 1L, length(plot$layers))) plot$layers[[j]]$ideogram_scope <- id
  }
  plot
}
