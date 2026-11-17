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
#' Locus annotations can inherit `x` as their bp position. Interval annotations
#' use `start` and `end`, and neither annotation uses the track's `y` mapping.
#' Put layers that should overlap in the same track. Two overlay tracks must
#' occupy separate bands on any chromosome they share. In horizontal local
#' synteny plots, tracks restricted to different chromosomes can reuse space.
#'
#' Pass an ordinary layer to `geom`, or several layers to `layers`. Compatible
#' identity-stat geoms keep their original Geom, Position, parameters and keys,
#' and `geom_hline()`, `geom_vline()` and `geom_rug()` project against each
#' chromosome's track boundaries. Discrete `y` (`factor` or `character`) maps
#' categories to equal sub-channels `1..K` with automatic category axis labels.
#' Area, ribbon, smooth (`geom_smooth()`), boxplot and violin layers finish
#' their native statistics and Position independently per chromosome in raw
#' bp/value coordinates before projection. Widths use bp; violin bandwidths use
#' the original value units. Automatic value limits include statistical geometry
#' such as confidence bands, density tails and boxplot notches.
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
#'   default width is 6 for gene-model tracks, 3 for other beside tracks,
#'   and 1 inside the body; gap is 0.2 beside and 0 inside. Use `NULL`
#'   to restore either default. Automatic dimensions follow the final side
#'   and contents. Explicit dimensions and offsets are retained on edits;
#'   adjust them if they conflict with the new side. Native geom widths
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
#'   In linear layouts, titles at a shared start align at their track-facing
#'   edge and clear any value axes at that end.
#' @param axis `FALSE`, `TRUE`, or a list of local value-axis arguments.
#'   `position` accepts `"auto"` (default), `"start"`, `"end"`, or `"gap"`.
#'   Auto places shared circular value axes in the closing gap; linear and
#'   explicitly chromosome-specific axes keep their end placement. A gap
#'   axis with per-chromosome ranges must select one chromosome with `chr`.
#'   Automatic axes skip chromosomes without a resolved value range. Use
#'   `labels = NULL` to hide all tick text. Empty strings or `NA` hide individual
#'   labels. Ticks remain visible in both cases.
#'   `gap` sets the distance from the track endpoint to the value-axis spine,
#'   in body-width units (default 0.2). Circular axes use at most half the
#'   adjacent sector gap, measured along the chromosome centerline.
#'   With `size = NULL` (default), tick text uses 1.8–2.8 mm based on the
#'   rendered axis length and label spacing. A numeric `size` fixes the
#'   physical text size. Set `family`, `colour`, `tick_length` (mm) and
#'   `label_gap` (em) independently.
#' @param chr Optional chromosome keys from the source karyotype. Known keys
#'   outside the current view are not drawn.
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
    side = c("right", "left", "overlay", "outer", "inner"), width = NULL,
    gap = NULL, value_scale = c("per_track", "global", "per_chr"),
    limits = NULL, transform = "identity", reverse = FALSE, offset = 0,
    layers = NULL, label = NULL, axis = FALSE, chr = NULL,
    clip = c("auto", "on", "off"), replace = FALSE, ...) {
  side <- match.arg(side)
  if (!is.null(geom)) {
    if (is.function(geom)) geom <- geom(...)
    else check_unused_args(...)
    layers <- c(list(geom), chr_track_layers(layers))
  } else {
    check_unused_args(...)
    layers <- chr_track_layers(layers)
  }
  spec <- ggideogram_track_spec(side, width %||% default_track_width(side, layers), gap,
    value_scale, limits, transform, reverse, offset, .partial = !is.null(track))
  spec$auto_width <- is.null(width)
  spec$auto_gap <- is.null(gap)
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
  spec$label <- label; spec$axis <- axis; spec$track <- track
  spec$chr <- if (!is.null(chr)) validate_chr(chr, "chr")
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
  m_ext <- attr(tracks, "marker_extent")
  m_ext_chr <- attr(tracks, "marker_extent_by_chr")
  if (!inherits(tracks, "ideogram_track_layout")) {
    if (!is.list(tracks)) stopf("`tracks` must be a named list of geom_track() declarations.")
    tracks <- do.call(track_layout, tracks)
  }
  for (id in names(tracks)) {
    specified <- tracks[[id]]$track
    if (!is.null(specified) && specified != id) stopf("Track `%s` has a different identifier `%s`.", id, specified)
    tracks[[id]] <- validate_track_spec(tracks[[id]])
  }
  attr(tracks, "marker_extent") <- m_ext
  attr(tracks, "marker_extent_by_chr") <- m_ext_chr
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

scope_annotation_mapping <- function(mapping, point = FALSE) {
  if (point && is.null(mapping$position) && !is.null(mapping$x)) mapping$position <- mapping$x
  track_mapping_without(mapping, c("x", "y", "value"))
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
        c("track", "chromosome_separation", "circular", "track_geometry", "zero_baseline"))]
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

build_track_hline_component <- function(layer, spec, id, layout) {
  own_data <- !is.null(layer$data) && !inherits(layer$data, "waiver")
  hdata <- if (own_data) {
    if (is.function(layer$data)) layer$data(spec$data) else layer$data
  } else spec$data
  hmapping <- if (isTRUE(layer$inherit.aes)) {
    scope_mapping(spec$mapping, layer$mapping, TRUE)
  } else layer$mapping
  yint_param <- layer$aes_params$yintercept %||% layer$geom_params$yintercept
  if (!is.data.frame(hdata)) {
    if (is.null(yint_param)) stopf("`geom_hline()` in track `%s` requires `yintercept`.", id)
    hdata <- data.frame(.hline_y = if (is_discrete_track_vec(yint_param)) as.character(yint_param) else as.numeric(yint_param),
      stringsAsFactors = FALSE)
    hmapping <- ggplot2::aes(yintercept = .data$.hline_y)
  }
  if (!nrow(hdata)) {
    empty_df <- data.frame(.hline_chr = character(0), .hline_pos = numeric(0),
      .hline_y = numeric(0), .hline_grp = character(0))
    proto <- ggplot2::geom_line(data = empty_df, show.legend = layer$show.legend, inherit.aes = FALSE)
    proto$mapping <- ggplot2::aes(chr = .data$.hline_chr, x = .data$.hline_pos,
      y = .data$.hline_y, group = .data$.hline_grp)
    comp <- native_track_component(proto, spec, id, layout)
    comp$from_hline <- TRUE
    return(comp)
  }
  yint <- if (!is.null(hmapping$yintercept)) {
    recycle_semantic(rlang::eval_tidy(hmapping$yintercept, data = hdata), nrow(hdata), "yintercept")
  } else if (!is.null(yint_param)) {
    recycle_semantic(yint_param, nrow(hdata), "yintercept")
  } else stopf("`geom_hline()` in track `%s` requires `yintercept`.", id)
  if (!is_discrete_track_vec(yint) && (!is.numeric(yint) || any(!is.finite(yint[!is.na(yint)])))) {
    stopf("Track `yintercept` must contain finite numbers or discrete track levels.")
  }
  owned_chrs <- track_chromosomes(layout, id)
  if (!is.null(hmapping) && "chr" %in% names(hmapping)) {
    row_chr <- validate_chr(recycle_semantic(rlang::eval_tidy(hmapping$chr, data = hdata), nrow(hdata), "chr"), "chr")
    source <- layout$data$karyotype
    if (any(!row_chr %in% source$.chr)) stopf("Track positions are outside the source chromosome bounds.")
    keep <- row_chr %in% owned_chrs & !is.na(yint)
    idx <- which(keep)
    pairs_r <- idx
    pairs_c <- row_chr[idx]
  } else {
    valid_r <- which(!is.na(yint))
    grid_rc <- expand.grid(r = valid_r, c = owned_chrs, stringsAsFactors = FALSE)
    pairs_r <- grid_rc$r
    pairs_c <- grid_rc$c
  }
  if (!length(pairs_r)) {
    empty_df <- data.frame(.hline_chr = character(0), .hline_pos = numeric(0),
      .hline_y = numeric(0), .hline_grp = character(0))
      proto <- ggplot2::geom_line(data = empty_df, show.legend = layer$show.legend, inherit.aes = FALSE)
    proto$mapping <- ggplot2::aes(chr = .data$.hline_chr, x = .data$.hline_pos,
      y = .data$.hline_y, group = .data$.hline_grp)
    comp <- native_track_component(proto, spec, id, layout)
    comp$from_hline <- TRUE
    return(comp)
  }
  g <- layout$chrom[match(pairs_c, layout$chrom$.chr), , drop = FALSE]
  rep2 <- rep(seq_along(pairs_r), each = 2L)
  out <- hdata[pairs_r[rep2], , drop = FALSE]
  out$.hline_chr <- pairs_c[rep2]
  out$.hline_pos <- as.vector(rbind(g$.start, g$.end))
  out$.hline_y <- if (is.factor(yint)) yint[pairs_r[rep2]] else yint[pairs_r[rep2]]
  out$.hline_grp <- paste0(pairs_r[rep2], ":", pairs_c[rep2])
  keep_aes <- setdiff(names(hmapping), c("yintercept", "chr", "x", "y", "position", "value", "group"))
  line_map <- combine_track_mapping(hmapping[keep_aes], ggplot2::aes(
    chr = .data$.hline_chr, x = .data$.hline_pos,
    y = .data$.hline_y, group = .data$.hline_grp))
  aes_p <- layer$aes_params
  aes_p$yintercept <- NULL
  geom_p <- layer$geom_params[intersect(names(layer$geom_params), c("na.rm", "lineend", "linejoin", "linemitre"))]
  proto <- do.call(ggplot2::geom_line, c(list(data = out,
    show.legend = layer$show.legend, inherit.aes = FALSE), geom_p, aes_p))
  proto$mapping <- line_map
  comp <- native_track_component(proto, spec, id, layout)
  comp$from_hline <- TRUE
  comp
}

build_track_vline_layer <- function(layer, spec, id, layout) {
  own_data <- !is.null(layer$data) && !inherits(layer$data, "waiver")
  vdata <- if (own_data) {
    if (is.function(layer$data)) layer$data(spec$data) else layer$data
  } else spec$data
  vmapping <- if (isTRUE(layer$inherit.aes)) {
    scope_mapping(spec$mapping, layer$mapping, TRUE)
  } else layer$mapping
  xint_param <- layer$aes_params$xintercept %||% layer$geom_params$xintercept
  if (!is.data.frame(vdata)) {
    if (is.null(xint_param)) stopf("`geom_vline()` in track `%s` requires `xintercept`.", id)
    vdata <- data.frame(.vline_x = as.numeric(xint_param))
    vmapping <- ggplot2::aes(xintercept = .data$.vline_x)
  }
  trow <- track_table_row(layout, id)
  if (!nrow(vdata)) {
    empty_seg <- data.frame(.vx = numeric(0), .vy = numeric(0), .vxend = numeric(0), .vyend = numeric(0))
    return(ggplot2::geom_segment(data = empty_seg,
      mapping = ggplot2::aes(x = .data$.vx, y = .data$.vy, xend = .data$.vxend, yend = .data$.vyend),
      inherit.aes = FALSE, show.legend = layer$show.legend))
  }
  xint <- if (!is.null(vmapping$xintercept)) {
    recycle_semantic(rlang::eval_tidy(vmapping$xintercept, data = vdata), nrow(vdata), "xintercept")
  } else if (!is.null(xint_param)) {
    recycle_semantic(xint_param, nrow(vdata), "xintercept")
  } else stopf("`geom_vline()` in track `%s` requires `xintercept`.", id)
  if (!is.numeric(xint) || anyNA(xint) || any(!is.finite(xint))) {
    stopf("Track `xintercept` must contain finite numeric coordinates.")
  }
  owned_chrs <- track_chromosomes(layout, id)
  if (!is.null(vmapping) && "chr" %in% names(vmapping)) {
    row_chr <- validate_chr(recycle_semantic(rlang::eval_tidy(vmapping$chr, data = vdata), nrow(vdata), "chr"), "chr")
    source <- layout$data$karyotype
    src_i <- match(row_chr, source$.chr)
    if (anyNA(src_i) || any(xint < source$.start[src_i] | xint > source$.end[src_i])) {
      stopf("Track positions are outside the source chromosome bounds.")
    }
    vis_i <- match(row_chr, layout$chrom$.chr)
    keep <- row_chr %in% owned_chrs & !is.na(vis_i) &
      xint >= layout$chrom$.start[vis_i] & xint <= layout$chrom$.end[vis_i]
    pairs_r <- which(keep)
    pairs_c <- row_chr[keep]
  } else {
    grid_rc <- expand.grid(r = seq_len(nrow(vdata)), c = owned_chrs, stringsAsFactors = FALSE)
    vis_i <- match(grid_rc$c, layout$chrom$.chr)
    keep <- !is.na(vis_i) & xint[grid_rc$r] >= layout$chrom$.start[vis_i] &
      xint[grid_rc$r] <= layout$chrom$.end[vis_i]
    pairs_r <- grid_rc$r[keep]
    pairs_c <- grid_rc$c[keep]
  }
  if (!length(pairs_r)) {
    empty_seg <- data.frame(.vx = numeric(0), .vy = numeric(0), .vxend = numeric(0), .vyend = numeric(0))
    return(ggplot2::geom_segment(data = empty_seg,
      mapping = ggplot2::aes(x = .data$.vx, y = .data$.vy, xend = .data$.vxend, yend = .data$.vyend),
      inherit.aes = FALSE, show.legend = layer$show.legend))
  }
  out <- vdata[pairs_r, , drop = FALSE]
  pts <- project_positions_checked(layout, pairs_c, xint[pairs_r])
  lo_pt <- offset_chr_points_signed(layout, pairs_c, pts, rep(trow$low_offset, length(pairs_c)))
  hi_pt <- offset_chr_points_signed(layout, pairs_c, pts, rep(trow$high_offset, length(pairs_c)))
  out$.vx <- lo_pt$x; out$.vy <- lo_pt$y
  out$.vxend <- hi_pt$x; out$.vyend <- hi_pt$y
  keep_aes <- setdiff(names(vmapping), c("xintercept", "chr", "x", "y", "xend", "yend", "position", "value"))
  seg_map <- combine_track_mapping(vmapping[keep_aes], ggplot2::aes(
    x = .data$.vx, y = .data$.vy, xend = .data$.vxend, yend = .data$.vyend))
  aes_p <- layer$aes_params
  aes_p$xintercept <- NULL
  geom_p <- layer$geom_params[intersect(names(layer$geom_params), c("na.rm", "lineend", "linejoin", "arrow", "arrow.fill"))]
  seg_layer <- do.call(ggplot2::geom_segment, c(list(mapping = seg_map, data = out,
    show.legend = layer$show.legend, inherit.aes = FALSE), geom_p, aes_p))
  seg_layer$ideogram_track_chrs <- unique(pairs_c)
  seg_layer
}

build_track_rug_layer <- function(layer, spec, id, layout) {
  own_data <- !is.null(layer$data) && !inherits(layer$data, "waiver")
  rdata <- if (own_data) {
    if (is.function(layer$data)) layer$data(spec$data) else layer$data
  } else spec$data
  rmapping <- if (isTRUE(layer$inherit.aes)) {
    scope_mapping(spec$mapping, layer$mapping, TRUE)
  } else layer$mapping
  rmapping <- scope_track_mapping(rmapping)
  if (!is.data.frame(rdata)) stopf("`geom_rug()` in track `%s` requires data-frame `data`.", id)
  if (is.null(rmapping) || !all(c("chr", "position") %in% names(rmapping))) {
    stopf("`geom_rug()` in track `%s` requires `chr` and `x` (or `position`) mappings.", id)
  }
  sides <- as.character(layer$geom_params$sides %||% layer$aes_params$sides %||% "b")[1]
  has_b <- grepl("b", sides, fixed = TRUE)
  has_t <- grepl("t", sides, fixed = TRUE)
  if (!has_b && !has_t) {
    stopf("`geom_rug()` in track `%s` requires `sides` containing `\"b\"` and/or `\"t\"`.", id)
  }
  raw_len <- layer$geom_params$length %||% layer$aes_params$length
  frac <- if (is.null(raw_len) || identical(as.character(raw_len), "0.03npc")) {
    0.15
  } else if (grid::is.unit(raw_len)) {
    val_u <- as.numeric(raw_len)[1]
    if (!is.finite(val_u) || val_u <= 0) stopf("`geom_rug()` `length` must be positive.")
    min(val_u, 1)
  } else if (is.numeric(raw_len) && length(raw_len) == 1L && is.finite(raw_len) && raw_len > 0) {
    min(raw_len, 1)
  } else {
    stopf("`geom_rug()` `length` in track `%s` must be a positive fraction or grid unit.", id)
  }
  rdata <- scope_point_data(rdata, rmapping, layout, spec$chr)
  empty_seg <- function() {
    empty_df <- data.frame(.rx = numeric(0), .ry = numeric(0), .rxend = numeric(0), .ryend = numeric(0))
    sl <- ggplot2::geom_segment(data = empty_df,
      mapping = ggplot2::aes(x = .data$.rx, y = .data$.ry, xend = .data$.rxend, yend = .data$.ryend),
      inherit.aes = FALSE, show.legend = layer$show.legend)
    sl$ideogram_track_chrs <- character(0)
    sl
  }
  if (!nrow(rdata)) return(empty_seg())
  f <- mapped_fields(rdata, rmapping, c("chr", "position"), what = "mapping")
  chr_v <- validate_chr(f$chr, "chr")
  pos_v <- validate_coordinate(f$position, "position")
  owned_chrs <- track_chromosomes(layout, id)
  keep <- chr_v %in% owned_chrs
  if (!any(keep)) return(empty_seg())
  rdata <- rdata[keep, , drop = FALSE]
  chr_v <- chr_v[keep]
  pos_v <- pos_v[keep]

  trow <- track_table_row(layout, id)
  bot_edge <- if (isTRUE(trow$reverse)) trow$high_offset else trow$low_offset
  top_edge <- if (isTRUE(trow$reverse)) trow$low_offset else trow$high_offset
  tick_delta <- (top_edge - bot_edge) * frac

  pieces <- list()
  if (has_b) {
    b_df <- rdata
    pts <- project_positions_checked(layout, chr_v, pos_v)
    p0 <- offset_chr_points_signed(layout, chr_v, pts, rep(bot_edge, length(chr_v)))
    p1 <- offset_chr_points_signed(layout, chr_v, pts, rep(bot_edge + tick_delta, length(chr_v)))
    b_df$.rx <- p0$x; b_df$.ry <- p0$y
    b_df$.rxend <- p1$x; b_df$.ryend <- p1$y
    pieces[[length(pieces) + 1L]] <- b_df
  }
  if (has_t) {
    t_df <- rdata
    pts <- project_positions_checked(layout, chr_v, pos_v)
    p0 <- offset_chr_points_signed(layout, chr_v, pts, rep(top_edge, length(chr_v)))
    p1 <- offset_chr_points_signed(layout, chr_v, pts, rep(top_edge - tick_delta, length(chr_v)))
    t_df$.rx <- p0$x; t_df$.ry <- p0$y
    t_df$.rxend <- p1$x; t_df$.ryend <- p1$y
    pieces[[length(pieces) + 1L]] <- t_df
  }
  out <- do.call(rbind, pieces)
  keep_aes <- setdiff(names(rmapping), c("chr", "x", "y", "xend", "yend", "position", "value"))
  seg_map <- combine_track_mapping(rmapping[keep_aes], ggplot2::aes(
    x = .data$.rx, y = .data$.ry, xend = .data$.rxend, yend = .data$.ryend))
  aes_p <- layer$aes_params
  aes_p[c("sides", "length", "outside")] <- NULL
  geom_p <- layer$geom_params[intersect(names(layer$geom_params), c("na.rm", "lineend"))]
  seg_layer <- do.call(ggplot2::geom_segment, c(list(mapping = seg_map, data = out,
    show.legend = layer$show.legend, inherit.aes = FALSE), geom_p, aes_p))
  seg_layer$ideogram_track_chrs <- unique(chr_v)
  seg_layer
}

native_track_component <- function(layer, spec, id, layout) {
  if (inherits(layer$geom, "GeomHline")) {
    return(build_track_hline_component(layer, spec, id, layout))
  }
  if (inherits(layer$geom, "GeomVline")) {
    return(build_track_vline_layer(layer, spec, id, layout))
  }
  if (inherits(layer$geom, "GeomRug")) {
    return(build_track_rug_layer(layer, spec, id, layout))
  }
  domain <- inherits(layer$stat, "StatChrMarker") || inherits(layer$stat, "StatChrInterval")
  if (domain && !is.null(layer$stat_params$track) && !identical(layer$stat_params$track, id)) {
    stopf("A layer inside track `%s` targets another track.", id)
  }
  inherit_scope <- if (domain) !isTRUE(layer$ideogram_no_scope_inherit) else isTRUE(layer$inherit.aes)
  mapping <- scope_mapping(spec$mapping, layer$mapping, inherit_scope)
  data <- scope_layer_data(layer$data, spec$data)
  if (domain) {
    mapping <- scope_annotation_mapping(mapping, point = inherits(layer$stat, "StatChrMarker"))
    prototype <- layer
    layer <- ggplot2::ggproto(NULL, prototype, data = data, mapping = mapping, inherit.aes = FALSE)
    layer$stat_params <- utils::modifyList(layer$stat_params, list(track = id))
    if (inherits(layer$stat, "StatChrInterval")) layer$stat_params$placement <- "beside"
    if (!is.null(layer$ideogram_marker_spec)) layer$ideogram_marker_spec$track <- id
    if ("position" %in% names(mapping)) layer$data <- scope_point_data(data, mapping, layout, spec$chr)
    if (inherits(layer$stat, "StatChrInterval")) {
      if (!is.data.frame(data)) stopf("Track layers require data-frame `data`.")
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
      if (!".source_start" %in% names(data)) data$.source_start <- first
      if (!".source_end" %in% names(data)) data$.source_end <- last
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
    if (inherits(layer$geom, "GeomSmooth") &&
      (inherits(layer$stat, "StatSmooth") || inherits(layer$stat, "StatIdentity"))) "smooth" else
    if (inherits(layer$geom, "GeomBoxplot") && inherits(layer$stat, "StatBoxplot")) "boxplot" else
    if (inherits(layer$geom, "GeomViolin") && inherits(layer$stat, "StatYdensity")) "violin" else
    if (inherits(layer$stat, "StatIdentity")) "dynamic" else NA_character_
  if (is.na(kind)) stopf("Track `%s` does not support the `%s` Stat; supply prepared values with an identity layer.",
    id, class(layer$stat)[1])
  range_layer <- kind %in% c("area", "ribbon", "smooth")
  area <- identical(kind, "area")
  mapping <- scope_track_mapping(mapping)
  is_rect <- inherits(layer$geom, "GeomRect")
  if (is.null(mapping$position)) {
    xmin_expr <- mapping$xmin %||% layer$aes_params$xmin %||% layer$geom_params$xmin
    xmax_expr <- mapping$xmax %||% layer$aes_params$xmax %||% layer$geom_params$xmax
    if (!is.null(xmin_expr) && !is.null(xmax_expr)) {
      mapping$position <- rlang::new_quosure(rlang::expr((!!xmin_expr + !!xmax_expr) / 2))
    }
  }
  if (is.null(mapping$value) && !range_layer) {
    ymin_expr <- mapping$ymin %||% layer$aes_params$ymin %||% layer$geom_params$ymin
    ymax_expr <- mapping$ymax %||% layer$aes_params$ymax %||% layer$geom_params$ymax
    if (!is.null(ymin_expr) && !is.null(ymax_expr)) {
      mapping$value <- rlang::new_quosure(rlang::expr((!!ymin_expr + !!ymax_expr) / 2))
    }
  }
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
  if (is_rect && all(c("xmin", "xmax") %in% unique(c(names(mapping), names(layer$aes_params), names(layer$geom_params))))) {
    if (!is.data.frame(data)) stopf("Track layers require data-frame `data`.")
    f_chr <- validate_chr(mapped_fields(data, mapping, "chr", what = "mapping")$chr, "chr")
    xmin_v <- if ("xmin" %in% names(mapping)) recycle_semantic(rlang::eval_tidy(mapping$xmin, data = data), nrow(data), "xmin") else
      recycle_semantic(layer$aes_params$xmin %||% layer$geom_params$xmin, nrow(data), "xmin")
    xmax_v <- if ("xmax" %in% names(mapping)) recycle_semantic(rlang::eval_tidy(mapping$xmax, data = data), nrow(data), "xmax") else
      recycle_semantic(layer$aes_params$xmax %||% layer$geom_params$xmax, nrow(data), "xmax")
    if (!is.numeric(xmin_v) || !is.numeric(xmax_v) || anyNA(xmin_v) || anyNA(xmax_v) ||
        any(!is.finite(xmin_v)) || any(!is.finite(xmax_v)) || any(xmin_v > xmax_v)) {
      stopf("Track rectangle `xmin` and `xmax` must be finite numbers with `xmin <= xmax`.")
    }
    source <- layout$data$karyotype
    src_i <- match(f_chr, source$.chr)
    if (anyNA(src_i) || any(xmin_v < source$.start[src_i] | xmax_v > source$.end[src_i])) {
      stopf("Track positions are outside the source chromosome bounds.")
    }
    vis_i <- match(f_chr, layout$chrom$.chr)
    keep <- !is.na(vis_i) & xmax_v > layout$chrom$.start[vis_i] & xmin_v < layout$chrom$.end[vis_i]
    if (!is.null(spec$chr)) keep <- keep & f_chr %in% spec$chr
    data <- data[keep, , drop = FALSE]
    if (nrow(data)) {
      vis_k <- vis_i[keep]
      c_lo <- pmax(xmin_v[keep], layout$chrom$.start[vis_k])
      c_hi <- pmin(xmax_v[keep], layout$chrom$.end[vis_k])
      data$.ideogram_rect_mid <- (c_lo + c_hi) / 2
      mapping$position <- ggplot2::aes(position = .data$.ideogram_rect_mid)$position
    }
  } else {
    data <- scope_point_data(data, mapping, layout, spec$chr)
  }
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
  object$native_layer <- layer
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
    part$side <- spec$side
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
      layer$mapping <- scope_mapping(part$mapping, layer$mapping, !isTRUE(layer$ideogram_no_scope_inherit))
      if (inherits(layer, "ggideogram_track_component")) {
        layer$mapping <- scope_track_mapping(layer$mapping)
      } else if (inherits(layer, c("ggideogram_labels_component", "ggideogram_repel_component"))) {
        layer$mapping <- scope_annotation_mapping(layer$mapping, point = TRUE)
      } else if (inherits(layer, c("ggideogram_gene_component", "ggideogram_chr_component"))) {
        layer$mapping <- scope_annotation_mapping(layer$mapping)
      }
      if (!is.null(part$chr)) {
        if (!is.data.frame(layer$data)) stopf("Track layers require data-frame `data`.")
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
      if (inherits(layer, c("ggideogram_track_component", "ggideogram_labels_component", "ggideogram_repel_component"))) {
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
    if (!is.null(spec$limits)) next
    for (key in names(spec$frozen_ranges)) {
      layout$track_ranges[[key]] <- spec$frozen_ranges[[key]]
    }
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
      if (object$kind == "area") {
        present_chr <- unique(chr)
        base_z <- track_zero_baseline(layout, object$track, values)
        if (is.finite(base_z) && length(present_chr)) {
          values <- c(values, rep(base_z, length(present_chr)))
          chr <- c(chr, present_chr)
        }
      }
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

track_scope_owners <- function(layout, id, contents) {
  if (!length(contents)) return(track_chromosomes(layout, id))
  owners <- unique(unlist(lapply(contents, function(object) {
    if (!is.null(object$ideogram_track_chrs)) return(object$ideogram_track_chrs)
    if (is.null(object$data) || !is.data.frame(object$data) || !nrow(object$data)) return(NULL)
    if (is.null(object$mapping) || !"chr" %in% names(object$mapping)) {
      return(track_chromosomes(layout, id))
    }
    chr_col <- as.character(mapped_fields(object$data, object$mapping, "chr", what = "mapping")$chr)
    if (all(c("start", "end") %in% names(object$mapping))) {
      f <- mapped_fields(object$data, object$mapping, c("start", "end"), what = "mapping")
      vis_i <- match(chr_col, layout$chrom$.chr)
      keep <- !is.na(vis_i) & f$end > layout$chrom$.start[vis_i] & f$start - 1 < layout$chrom$.end[vis_i]
      chr_col <- chr_col[keep]
    }
    chr_col
  }), use.names = FALSE))
  assigned <- track_chromosomes(layout, id)
  intersect(assigned, owners)
}

add_track_scope_title <- function(plot, id, title, owners, circular_gap = 0) {
  layout <- ideogram_plot_layout(plot)
  settings <- utils::modifyList(list(size = ideogram_text_size("track"),
    colour = ideogram_colour("text"), family = layout$base_spec$base_family %||% "",
    fontface = 1, gap = 0.3),
    Filter(Negate(is.null), track_label_settings(title)))
  if (is.null(settings$text) || !nzchar(settings$text) || !length(owners)) return(plot)
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
    arc <- circular_gap_axis_arc(layout, circular_gap)
    point <- circular_xy(layout, arc, p$y)
    theta <- circular_theta(layout, arc)
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
  layer$ideogram_track_title <- list(chr = g$.chr, track = id, anchor = data[c("x", "y")])
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
  encoded <- encode_discrete_track_scopes(plot, tracks, components)
  plot <- encoded$plot
  components <- encoded$components
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
    owners <- track_scope_owners(layout, id, components[[i]])
    axis <- tracks[[i]]$axis
    axis_owners <- owners
    args <- if (is.list(axis)) axis else if (isTRUE(axis)) list()
    circular_gap <- 0
    if (!is.null(args)) {
      if (is.null(args$chr)) {
        spec <- track_table_row(layout, id)
        if (identical(args$position, "gap") && spec$value_scale == "per_chr") {
          stopf("A gap axis with `value_scale = \"per_chr\"` must select one chromosome with `chr`.")
        }
        keys <- track_key_for_rows(spec, owners)
        axis_owners <- owners[vapply(keys, function(key)
          !is.null(layout$track_ranges[[key]]$limits), logical(1))]
        if (is_circular_layout(layout) && tracks[[i]]$value_scale != "per_chr" &&
            (args$position %||% "auto") == "auto") args$position <- "gap"
        args$chr <- if (identical(args$position, "gap")) axis_owners[1] else axis_owners
      }
      position <- args$position %||% "auto"
      if (is_circular_layout(layout) && length(axis_owners) && (position == "gap" ||
          (position == "auto" && isTRUE(args$chr) && tracks[[i]]$value_scale != "per_chr"))) {
        circular_gap <- track_axis_gap(args$gap)
      }
    }
    before <- length(plot$layers)
    plot <- add_track_scope_title(plot, id, tracks[[i]]$label, owners, circular_gap)
    if (length(plot$layers) > before) plot$layers[[length(plot$layers)]]$ideogram_scope_slot <- "title"
    if (!is.null(args) && length(axis_owners)) {
      before <- length(plot$layers)
      plot <- plot + do.call(geom_track_axis, c(list(track = id), args))
      if (length(plot$layers) > before) for (j in seq.int(before + 1L, length(plot$layers))) {
        plot$layers[[j]]$ideogram_scope_slot <- paste0("axis:", j - before)
      }
    }
    if (length(plot$layers) > start) for (j in seq.int(start + 1L, length(plot$layers))) plot$layers[[j]]$ideogram_scope <- id
  }
  plot
}
