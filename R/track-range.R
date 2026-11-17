# Training and validation of track value ranges.

is_discrete_track_vec <- function(x) {
  is.factor(x) || is.character(x) || (is.logical(x) && any(!is.na(x)))
}

extract_discrete_levels <- function(x, present_only = FALSE) {
  if (is.factor(x)) {
    lv <- as.character(levels(x))
    lv <- lv[!is.na(lv)]
    if (present_only) {
      obs <- unique(as.character(x[!is.na(x)]))
      lv <- intersect(lv, obs)
    }
    return(lv)
  }
  unique(as.character(x[!is.na(x)]))
}

encode_discrete_track_scopes <- function(plot, tracks, components) {
  layout <- ideogram_plot_layout(plot)
  track_ids <- names(tracks)
  if (!length(track_ids)) return(list(plot = plot, components = components))

  global_idx <- which(vapply(tracks, function(sp) identical(sp$value_scale, "global"), logical(1)))
  other_idx <- setdiff(seq_along(tracks), global_idx)
  groups <- c(if (length(global_idx)) list(global_idx) else list(), as.list(other_idx))

  for (grp in groups) {
    has_disc <- FALSE
    has_cont <- FALSE
    first_id <- track_ids[grp[1]]
    for (ti in grp) {
      for (obj in components[[ti]]) {
        if (!inherits(obj, "ggideogram_track_component")) next
        if ("value" %in% names(obj$mapping)) {
          val <- mapped_fields(obj$data, obj$mapping, "value", what = "mapping")$value
          if (is_discrete_track_vec(val)) {
            has_disc <- TRUE
          } else if (is.numeric(val) && !isTRUE(obj$from_hline)) {
            has_cont <- TRUE
          }
        }
        for (a in intersect(c("ymin", "ymax", "yend"), union(names(obj$mapping), names(obj$params)))) {
          av <- if (a %in% names(obj$params)) obj$params[[a]] else
            rlang::eval_tidy(obj$mapping[[a]], data = obj$data)
          if (is_discrete_track_vec(av)) has_disc <- TRUE
        }
      }
    }
    if (!has_disc) next
    if (has_cont) {
      stopf("Track `%s` cannot mix discrete and continuous `y` values.", first_id)
    }

    for (ti in grp) {
      tid <- track_ids[ti]
      trow <- track_table_row(layout, tid)
      if (!identical(trow$transform_name, "identity")) {
        stopf("Discrete `y` values in track `%s` require `transform = \"identity\"`.", tid)
      }
      for (obj in components[[ti]]) {
        if (!inherits(obj, "ggideogram_track_component")) next
        if (obj$kind != "dynamic" || !is.null(obj$range_prototype)) {
          stopf("Discrete `y` values in track `%s` require identity track geoms (such as geom_tile, geom_rect, geom_point or geom_text).", tid)
        }
      }
    }

    is_per_chr <- length(grp) == 1L && identical(tracks[[grp[1]]]$value_scale, "per_chr")
    if (!is_per_chr) {
      lvl_parts <- list()
      for (ti in grp) {
        for (obj in components[[ti]]) {
          if (!inherits(obj, "ggideogram_track_component")) next
          if ("value" %in% names(obj$mapping)) {
            val <- mapped_fields(obj$data, obj$mapping, "value", what = "mapping")$value
            if (is_discrete_track_vec(val)) lvl_parts[[length(lvl_parts) + 1L]] <- extract_discrete_levels(val)
          }
          for (a in intersect(c("ymin", "ymax", "yend"), union(names(obj$mapping), names(obj$params)))) {
            av <- if (a %in% names(obj$params)) obj$params[[a]] else
              rlang::eval_tidy(obj$mapping[[a]], data = obj$data)
            if (is_discrete_track_vec(av)) lvl_parts[[length(lvl_parts) + 1L]] <- extract_discrete_levels(av)
          }
        }
      }
      levels_all <- unique(unlist(lvl_parts, use.names = FALSE))
      for (ti in grp) {
        tid <- track_ids[ti]
        trow <- track_table_row(layout, tid)
        keys <- track_range_keys(trow, layout$chrom$.chr)
        if (length(levels_all)) {
          lim <- trow$limits[[1]] %||% c(0.5, length(levels_all) + 0.5)
          for (key in keys) {
            layout$track_ranges[[key]] <- list(limits = lim, fixed = TRUE, discrete_levels = levels_all)
          }
        }
        for (ci in seq_along(components[[ti]])) {
          obj <- components[[ti]][[ci]]
          if (!inherits(obj, "ggideogram_track_component")) next
          components[[ti]][[ci]] <- encode_one_discrete_component(obj, function(chr_v, v) {
            match(as.character(v), levels_all)
          })
        }
      }
    } else {
      tid <- track_ids[grp[1]]
      trow <- track_table_row(layout, tid)
      chr_levels_map <- list()
      for (ch in layout$chrom$.chr) {
        parts_c <- list()
        for (obj in components[[grp[1]]]) {
          if (!inherits(obj, "ggideogram_track_component") || !nrow(obj$data)) next
          f_c <- as.character(mapped_fields(obj$data, obj$mapping, "chr", what = "mapping")$chr)
          if (!any(f_c == ch)) next
          if ("value" %in% names(obj$mapping)) {
            val <- mapped_fields(obj$data, obj$mapping, "value", what = "mapping")$value
            if (is_discrete_track_vec(val)) {
              sub_v <- if (is.factor(val)) val[f_c == ch] else as.character(val[f_c == ch])
              parts_c[[length(parts_c) + 1L]] <- extract_discrete_levels(sub_v, present_only = TRUE)
            }
          }
        }
        lv_c <- unique(unlist(parts_c, use.names = FALSE))
        if (length(lv_c)) {
          chr_levels_map[[ch]] <- lv_c
          key <- paste0("track:", tid, ":chr:", ch)
          lim <- trow$limits[[1]] %||% c(0.5, length(lv_c) + 0.5)
          layout$track_ranges[[key]] <- list(limits = lim, fixed = TRUE, discrete_levels = lv_c)
        }
      }
      for (ci in seq_along(components[[grp[1]]])) {
        obj <- components[[grp[1]]][[ci]]
        if (!inherits(obj, "ggideogram_track_component")) next
        components[[grp[1]]][[ci]] <- encode_one_discrete_component(obj, function(chr_v, v) {
          out <- rep(NA_real_, length(v))
          for (ch in unique(chr_v)) {
            lv_c <- chr_levels_map[[ch]]
            if (!is.null(lv_c)) {
              idx <- which(chr_v == ch)
              out[idx] <- match(as.character(v[idx]), lv_c)
            }
          }
          out
        })
      }
    }
  }
  list(plot = update_plot_ideogram_layout(plot, layout), components = components)
}

encode_one_discrete_component <- function(obj, encoder) {
  n <- nrow(obj$data)
  chr_v <- if (n && "chr" %in% names(obj$mapping)) {
    as.character(mapped_fields(obj$data, obj$mapping, "chr", what = "mapping")$chr)
  } else character(0)
  if ("value" %in% names(obj$mapping)) {
    val <- if (n) mapped_fields(obj$data, obj$mapping, "value", what = "mapping")$value else character(0)
    if (is_discrete_track_vec(val) || !n) {
      codes <- if (n) as.numeric(encoder(chr_v, val)) else numeric(0)
      obj$data$.ideogram_disc_y <- codes
      q_y <- ggplot2::aes(value = .data$.ideogram_disc_y)$value
      obj$mapping$value <- q_y
      if ("y" %in% names(obj$mapping)) obj$mapping$y <- q_y
    }
  }
  for (a in intersect(c("ymin", "ymax", "yend"), union(names(obj$mapping), names(obj$params)))) {
    av <- if (a %in% names(obj$params)) obj$params[[a]] else
      if (n) rlang::eval_tidy(obj$mapping[[a]], data = obj$data) else character(0)
    if (is_discrete_track_vec(av)) {
      av <- recycle_semantic(av, n, a)
      col_nm <- paste0(".ideogram_disc_", a)
      obj$data[[col_nm]] <- if (n) as.numeric(encoder(chr_v, av)) else numeric(0)
      obj$mapping[[a]] <- rlang::new_quosure(rlang::expr(.data[[!!col_nm]]), env = rlang::base_env())
      obj$params[[a]] <- NULL
    }
  }
  if (inherits(obj$native_geom, "GeomTile") &&
      !"height" %in% union(names(obj$mapping), names(obj$params))) {
    obj$params$height <- 1
    if (!is.null(obj$native_layer)) obj$native_layer$aes_params$height <- 1
  }
  if (inherits(obj$native_geom, "GeomRect") &&
      !all(c("ymin", "ymax") %in% union(names(obj$mapping), names(obj$params))) &&
      ".ideogram_disc_y" %in% names(obj$data)) {
    h <- obj$params$height %||% 0.9
    obj$params$height <- NULL
    if (!is.null(obj$native_layer)) {
      obj$native_layer$geom_params$height <- NULL
      obj$native_layer$aes_params$height <- NULL
    }
    obj$data$.ideogram_disc_ymin <- obj$data$.ideogram_disc_y - h / 2
    obj$data$.ideogram_disc_ymax <- obj$data$.ideogram_disc_y + h / 2
    obj$mapping$ymin <- ggplot2::aes(ymin = .data$.ideogram_disc_ymin)$ymin
    obj$mapping$ymax <- ggplot2::aes(ymax = .data$.ideogram_disc_ymax)$ymax
  }
  obj$include_zero <- FALSE
  obj
}

register_dynamic_track_range <- function(object, layout) {
  fields <- mapped_fields(
    object$data, object$mapping, c("chr", "position", "value"),
    what = "mapping")
  validate_track_semantics(layout, fields$chr, fields$position, fields$value)
  chr <- as.character(fields$chr)
  value <- as.numeric(fields$value)
  native_range <- !is.null(object$native_layer) &&
    ((inherits(object$native_geom, "GeomTile") &&
        inherits(object$position, "PositionIdentity")) ||
      inherits(object$position, c("PositionStack", "PositionNudge")))
  if (!is.null(object$range_prototype)) {
    computed <- native_range_track_values(
      object$range_prototype, object$data, object$mapping, layout, object$track)
    chr <- computed$chr
    value <- computed$value
  } else if (native_range) {
    computed <- native_identity_track_values(object, layout, fields)
    chr <- computed$chr
    value <- computed$value
  } else {
    tile_range <- track_tile_value_range(object, value)
    if (!is.null(tile_range)) {
      chr <- rep(chr, 3L)
      value <- c(value, tile_range$lower, tile_range$upper)
    }
    stacked <- track_position_value_range(object$position, chr, fields$position, value)
    if (!is.null(stacked)) {
      chr <- c(chr, stacked$chr)
      value <- c(value, stacked$value)
    }
  }
  if (isTRUE(object$include_zero) && !native_range && is.null(object$range_prototype)) {
    present_chr <- unique(chr)
    base_z <- track_zero_baseline(layout, object$track, value)
    if (is.finite(base_z) && length(present_chr)) {
      chr <- c(chr, present_chr)
      value <- c(value, rep(base_z, length(present_chr)))
    }
  }
  coordinate_names <- intersect(c("xend", "xmin", "xmax", "yend", "ymin", "ymax"),
    union(names(object$mapping), names(object$params)))
  coordinates <- lapply(stats::setNames(coordinate_names, coordinate_names), function(aesthetic) {
    input <- if (aesthetic %in% names(object$params)) object$params[[aesthetic]] else
      rlang::eval_tidy(object$mapping[[aesthetic]], data = object$data)
    recycle_semantic(input, nrow(object$data), paste0("track ", aesthetic))
  })
  is_rect_geom <- inherits(object$native_geom, "GeomRect")
  for (aesthetic in intersect(c("xend", "xmin", "xmax"), coordinate_names)) {
    if (is_rect_geom && aesthetic %in% c("xmin", "xmax")) {
      coord_x <- coordinates[[aesthetic]]
      if (!is.numeric(coord_x) || anyNA(coord_x) || any(!is.finite(coord_x))) {
        stopf("Track `%s` must contain finite numeric coordinates.", aesthetic)
      }
      source <- layout$data$karyotype
      src_i <- match(as.character(fields$chr), source$.chr)
      if (anyNA(src_i) || any(coord_x < source$.start[src_i] | coord_x > source$.end[src_i])) {
        stopf("Track positions are outside the source chromosome bounds.")
      }
    } else {
      validate_track_semantics(layout, fields$chr, coordinates[[aesthetic]], fields$value)
    }
  }
  for (aesthetic in intersect(c("yend", "ymin", "ymax"), coordinate_names)) {
    validate_track_semantics(layout, fields$chr, fields$position, coordinates[[aesthetic]])
    if (!native_range) {
      chr <- c(chr, as.character(fields$chr))
      value <- c(value, coordinates[[aesthetic]])
    }
  }
  layout <- register_track_values(layout, object$track, chr, value)
  list(layout = layout, coordinates = coordinates, coordinate_names = coordinate_names)
}

track_tile_value_range <- function(object, value) {
  if (!identical(object$geom, ggplot2::geom_tile) &&
      !inherits(object$native_geom, "GeomTile")) return(NULL)
  height <- if ("height" %in% names(object$mapping)) {
    mapped <- rlang::eval_tidy(object$mapping$height, data = object$data)
    recycle_semantic(mapped, nrow(object$data), "mapping$height")
  } else if (!is.null(object$params$height)) {
    rep(object$params$height, nrow(object$data))
  } else {
    rep(ggplot2::resolution(value, zero = FALSE), length(value))
  }
  if (!is.numeric(height) || anyNA(height) || any(!is.finite(height)) ||
      any(height < 0)) {
    stopf("Tile `height` must contain finite non-negative numbers.")
  }
  list(lower = value - height / 2, upper = value + height / 2)
}

# Use the native Geom and Position to train their finished value geometry.
# The same bp separation as StatChrTrack keeps native widths and positions
# independent across chromosomes without changing the source observations.
native_identity_track_values <- function(object, layout, fields) {
  data <- object$data
  chr <- as.character(fields$chr)
  shift <- (match(chr, unique(chr)) - as.integer(!is_circular_layout(layout))) *
    track_chromosome_separation(layout, data, object$mapping, object$params)
  mapping <- track_mapping_without(object$mapping, c("position", "value"))
  mapping$x <- fields$position + shift
  mapping$y <- fields$value
  for (field in intersect(c("xend", "xmin", "xmax"),
      union(names(mapping), names(object$params)))) {
    value <- object$params[[field]] %||%
      rlang::eval_tidy(mapping[[field]], data = data)
    mapping[[field]] <- recycle_semantic(value, nrow(data), field) + shift
  }
  prototype <- object$native_layer
  layer <- ggplot2::ggproto(NULL, prototype, data = data,
    mapping = mapping, inherit.aes = FALSE)
  layer$aes_params <- prototype$aes_params
  layer$aes_params[c("xend", "xmin", "xmax")] <- NULL
  layer$stat <- ggplot2::ggproto(NULL, prototype$stat,
    required_aes = union(prototype$stat$required_aes, "chr"))
  computed <- ggplot2::ggplot_build(ggplot2::ggplot() + layer)$data[[1]]
  columns <- intersect(c("y", "ymin", "ymax", "yend"), names(computed))
  vals <- unlist(computed[columns], use.names = FALSE)
  if (isTRUE(object$include_zero)) {
    base_z <- track_zero_baseline(layout, object$track, fields$value)
    if (is.finite(base_z) && base_z != 0 && !any(fields$value == 0, na.rm = TRUE)) {
      vals[!is.na(vals) & vals == 0] <- base_z
    }
  }
  list(chr = rep(as.character(computed$chr), length(columns)),
    value = vals)
}

track_chromosome_separation <- function(layout, data, mapping, params) {
  genomic_span <- max(layout$chrom$.end) - min(layout$chrom$.start)
  width <- if ("width" %in% names(mapping)) {
    rlang::eval_tidy(mapping$width, data = data)
  } else {
    params$width %||% 0
  }
  width <- suppressWarnings(max(abs(width), na.rm = TRUE))
  if (!is.finite(width)) width <- 0
  separation <- genomic_span + max(width, genomic_span * sqrt(.Machine$double.eps))
  if (is_circular_layout(layout)) separation <- max(separation, layout$circular$semantic_separation)
  separation
}

track_position_value_range <- function(position, chr, genomic_position, value) {
  is_fill <- identical(position, "fill") || inherits(position, "PositionFill")
  is_stack <- identical(position, "stack") ||
    inherits(position, "PositionStack")
  if (!is_fill && !is_stack) return(NULL)
  if (is_fill) {
    return(data.frame(chr = unique(as.character(chr)), value = 1))
  }
  key <- interaction(chr, genomic_position, drop = TRUE, lex.order = TRUE)
  positive <- tapply(pmax(value, 0, na.rm = FALSE), key, sum, na.rm = TRUE)
  negative <- tapply(pmin(value, 0, na.rm = FALSE), key, sum, na.rm = TRUE)
  first <- match(levels(key), key)
  data.frame(
    chr = rep(as.character(chr[first]), each = 2L),
    value = as.vector(rbind(negative, positive)),
    stringsAsFactors = FALSE
  )
}

validate_track_semantics <- function(layout, chr, position, value) {
  chr <- as.character(chr)
  match_layout_chr(layout, chr)
  if (!is.numeric(position) || anyNA(position) ||
      any(!is.finite(position))) {
    stopf("Track `position` must contain finite numeric coordinates.")
  }
  project_positions_checked(layout, chr, position)
  if (!is.numeric(value)) {
    stopf("Track `value` must be numeric.")
  }
  finite <- !is.na(value)
  if (any(!is.finite(value[finite]))) {
    stopf("Track `value` must contain only finite numbers or missing values.")
  }
  invisible(TRUE)
}

require_explicit_track_limits <- function(layout, track_id, component) {
  spec <- track_table_row(layout, track_id)
  if (is.null(spec$limits[[1]]) && !track_id %in% (layout$scope_fixed_tracks %||% character())) {
    stopf(paste0(
      "`%s` requires explicit `limits` in `geom_track(%s = ...)` because its ",
      "standard ggplot2 Stat/Geom reconstructs coordinates before drawing."),
      component, "limits")
  }
  spec
}
