# Rebuild source components after a named track changes.

copy_ideogram_layer <- function(layer) {
  # ggproto retains its parent; force it before a caller reuses a loop variable.
  force(layer)
  ggplot2::ggproto(NULL, layer)
}

record_ideogram_aesthetics <- function(layer) {
  layer$ideogram_aesthetic_names <- list(mapping = names(layer$mapping),
    aes_params = names(layer$aes_params), defaults = names(layer$geom$default_aes))
  layer
}

preserve_ideogram_aesthetics <- function(layer, previous) {
  if (!inherits(previous$geom, class(layer$geom)[1]) &&
      !inherits(layer$geom, class(previous$geom)[1])) return(layer)
  original <- previous$ideogram_aesthetic_names
  current <- list(mapping = names(previous$mapping),
    aes_params = names(previous$aes_params), defaults = names(previous$geom$default_aes))
  # Keep the effective native extension prototype when it has renamed an
  # aesthetic. Source geometry and semantic Stat parameters still come from
  # the newly computed layout.
  pairs <- lapply(names(original), function(field) {
    before <- original[[field]]
    after <- utils::head(current[[field]], length(before))
    if (length(before) != length(after)) return(character())
    stats::setNames(after[before != after], before[before != after])
  })
  renamed <- do.call(c, pairs)
  if (!length(renamed)) return(layer)
  renamed <- renamed[!duplicated(names(renamed))]
  rename <- function(values) {
    index <- match(names(values), names(renamed))
    names(values)[!is.na(index)] <- unname(renamed[index[!is.na(index)]])
    values
  }
  effective <- ggplot2::ggproto(NULL, previous)
  for (field in c("data", "position", "stat_params", "geom_params", "show.legend",
      "inherit.aes", "ideogram_base", "ideogram_scope", "ideogram_scope_slot",
      "ideogram_recipe_id", "ideogram_recipe", "ideogram_aesthetic_names")) {
    effective[[field]] <- layer[[field]]
  }
  effective$mapping <- rename(layer$mapping)
  effective$aes_params <- rename(layer$aes_params)
  effective
}

clear_global_frozen_ranges <- function(specs) {
  for (key in names(specs)) {
    if (identical(specs[[key]]$value_scale, "global")) {
      specs[[key]]$frozen_ranges <- NULL
    }
  }
  specs
}

#' @method ggplot_add ggideogram_track_scope
#' @export
ggplot_add.ggideogram_track_scope <- function(object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  if (is.null(object$track)) {
    stopf("Supply a unique `track` identifier when adding geom_track() with +.")
  }
  specs <- layout$base_spec$tracks %||% track_layout()
  id <- object$track
  if (id %in% names(specs)) {
    previous <- specs[[id]]
    old_scale <- previous$value_scale
    for (arg in object$geometry_args) previous[arg] <- object[arg]
    supplied <- object$supplied_args
    if ("width" %in% supplied) previous$auto_width <- object$auto_width
    if ("gap" %in% supplied) previous$auto_gap <- object$auto_gap
    for (arg in intersect(supplied, c("chr", "clip"))) previous[arg] <- object[arg]
    if ("label" %in% supplied) {
      previous$label <- if (is.null(object$label)) NULL else {
        utils::modifyList(track_label_settings(previous$label),
          track_label_settings(object$label), keep.null = TRUE)
      }
    }
    if ("axis" %in% supplied) {
      previous$axis <- if (is.list(object$axis)) {
        settings <- if (is.list(previous$axis)) previous$axis else list()
        utils::modifyList(settings, object$axis, keep.null = TRUE)
      } else object$axis
    }
    if (isTRUE(object$replace)) {
      previous$parts <- NULL
      previous$layers <- object$layers
      for (arg in intersect(supplied, c("data", "mapping"))) previous[arg] <- object[arg]
    } else if (length(object$layers)) {
      part <- object
      was_empty <- !length(previous$layers) && is.null(previous$parts)
      if (was_empty && is.null(previous$data) && "data" %in% supplied) {
        previous$data <- object$data
        part$inherit_data <- TRUE
      } else {
        part$inherit_data <- !"data" %in% supplied
      }
      if (was_empty && is.null(previous$mapping) && "mapping" %in% supplied) {
        previous$mapping <- object$mapping
        part$inherit_mapping <- TRUE
      } else {
        part$inherit_mapping <- !"mapping" %in% supplied
      }
      previous$parts <- if (was_empty) list(part) else
        c(previous$parts %||% list(previous), list(part))
    } else {
      for (arg in intersect(supplied, c("data", "mapping"))) previous[arg] <- object[arg]
    }
    range_changed <- any(c("limits", "value_scale", "transform", "chr") %in% supplied)
    data_changed <- !length(object$layers) && any(c("data", "mapping") %in% supplied)
    if (isTRUE(object$replace) || range_changed || data_changed) {
      previous$frozen_ranges <- NULL
      if (identical(previous$value_scale, "global") || identical(old_scale, "global")) {
        specs <- clear_global_frozen_ranges(specs)
      }
    }
    specs[[id]] <- validate_track_spec(previous)
  } else {
    if (identical(object$value_scale, "global")) {
      specs <- clear_global_frozen_ranges(specs)
    }
    specs[[id]] <- object
  }
  spec <- layout$base_spec
  spec$tracks <- specs
  rebuild_ideogram(plot, spec, id, isTRUE(object$replace))
}

is_body_base_layer <- function(layer) {
  isTRUE(layer$ideogram_base) && is.null(layer$ideogram_name_gap) &&
    is.null(layer$ideogram_axis_spec) && is.null(layer$ideogram_scope) &&
    !isTRUE(layer$ideogram_serial_row_label)
}

rebuild_ideogram <- function(plot, spec, id = NULL, replace = FALSE) {
  layout <- ideogram_plot_layout(plot)
  omit_base_body <- isTRUE(spec$.omit_base_body)
  call_spec <- spec
  call_spec$.omit_base_body <- NULL
  call_spec$data <- layout$data
  rebuilt <- do.call(ggideogram, call_spec)
  updated <- ideogram_plot_layout(rebuilt)
  updated$base_spec$.omit_base_body <- if (omit_base_body) TRUE else NULL
  updated$attachments <- layout$attachments
  metadata <- intersect(setdiff(names(layout$chrom), names(updated$chrom)),
    names(updated$data$karyotype))
  for (column in metadata) {
    updated$chrom[[column]] <- updated$data$karyotype[[column]][
      match(updated$chrom$.chr, updated$data$karyotype$.chr)]
  }
  rebuilt <- update_plot_ideogram_layout(rebuilt, updated)
  generated <- rebuilt$layers
  if (omit_base_body) {
    generated <- Filter(function(l) !is_body_base_layer(l), generated)
  }
  used <- rep(FALSE, length(generated))
  rebuilt$layers <- list()
  base <- which(vapply(generated, function(layer) isTRUE(layer$ideogram_base), logical(1)))
  base_index <- 0L
  seen <- integer()
  visited_slots <- character()
  last_replace_content <- 0L
  if (replace && !is.null(id)) {
    rc_idx <- which(vapply(plot$layers, function(l) {
      identical(l$ideogram_scope, id) && is.character(l$ideogram_scope_slot) &&
        startsWith(l$ideogram_scope_slot, "content:")
    }, logical(1)))
    if (length(rc_idx)) last_replace_content <- max(rc_idx)
  }
  for (pos in seq_along(plot$layers)) {
    layer <- plot$layers[[pos]]
    if (isTRUE(layer$ideogram_base)) {
      base_index <- base_index + 1L
      index <- base[base_index]
      if (!is.na(index)) {
        rebuilt$layers <- c(rebuilt$layers, list(preserve_ideogram_aesthetics(generated[[index]], layer)))
        used[index] <- TRUE
      }
      next
    }
    if (!is.null(layer$ideogram_scope)) {
      scope_id <- layer$ideogram_scope
      index <- which(!used & vapply(generated, function(candidate) {
        identical(candidate$ideogram_scope, scope_id) &&
          identical(candidate$ideogram_scope_slot, layer$ideogram_scope_slot)
      }, logical(1)))
      if (length(index)) {
        target <- index[1]
        fresh <- generated[[target]]
        replacing <- identical(scope_id, id) && replace &&
          startsWith(layer$ideogram_scope_slot, "content:")
        if (!replacing) fresh <- preserve_ideogram_aesthetics(fresh, layer)
        rebuilt$layers <- c(rebuilt$layers, list(fresh))
        used[target] <- TRUE
        visited_slots <- union(visited_slots, paste0(scope_id, "::", layer$ideogram_scope_slot))
      }
      if (pos == last_replace_content) {
        extra_content <- which(!used & vapply(generated, function(candidate) {
          identical(candidate$ideogram_scope, id) &&
            is.character(candidate$ideogram_scope_slot) &&
            startsWith(candidate$ideogram_scope_slot, "content:")
        }, logical(1)))
        if (length(extra_content)) {
          rebuilt$layers <- c(rebuilt$layers, generated[extra_content])
          used[extra_content] <- TRUE
          for (ec in extra_content) {
            visited_slots <- union(visited_slots,
              paste0(id, "::", generated[[ec]]$ideogram_scope_slot))
          }
        }
      }
      next
    }
    if (!is.null(layer$ideogram_recipe_id)) {
      if (layer$ideogram_recipe_id %in% seen) next
      seen <- c(seen, layer$ideogram_recipe_id)
      prev_tracks <- ideogram_plot_layout(rebuilt)$base_spec$tracks
      prev_markers <- ideogram_plot_layout(rebuilt)$marker_extent
      before <- length(rebuilt$layers)
      rebuilt <- rebuilt + new_ideogram_object(layer$ideogram_recipe)
      cur_layout <- ideogram_plot_layout(rebuilt)
      if (!identical(cur_layout$base_spec$tracks, prev_tracks) ||
          !identical(cur_layout$marker_extent, prev_markers)) {
        keep_layers <- vapply(rebuilt$layers, function(candidate) {
          is.null(candidate$ideogram_scope) ||
            paste0(candidate$ideogram_scope, "::", candidate$ideogram_scope_slot) %in% visited_slots
        }, logical(1))
        rebuilt$layers <- rebuilt$layers[keep_layers]
        refresh_spec <- cur_layout$base_spec
        refresh_spec$.omit_base_body <- NULL
        refresh_spec$data <- cur_layout$data
        generated <- do.call(ggideogram, refresh_spec)$layers
        if (omit_base_body) {
          generated <- Filter(function(l) !is_body_base_layer(l), generated)
        }
        base <- which(vapply(generated, function(g) isTRUE(g$ideogram_base), logical(1)))
        used <- vapply(seq_along(generated), function(gi) {
          cand <- generated[[gi]]
          if (isTRUE(cand$ideogram_base)) return(match(gi, base) <= base_index)
          !is.null(cand$ideogram_scope) &&
            paste0(cand$ideogram_scope, "::", cand$ideogram_scope_slot) %in% visited_slots
        }, logical(1))
      }
      previous <- Filter(function(candidate) identical(candidate$ideogram_recipe_id,
        layer$ideogram_recipe_id), plot$layers)
      count <- min(length(previous), length(rebuilt$layers) - before)
      for (j in seq_len(count)) {
        rebuilt$layers[[before + j]] <-
          preserve_ideogram_aesthetics(rebuilt$layers[[before + j]], previous[[j]])
      }
    } else if (!is.null(layer$ideogram_axis_spec)) {
      axes <- do.call(chromosome_axis_layers,
        c(list(layout = ideogram_plot_layout(rebuilt)), layer$ideogram_axis_spec))
      if (length(axes) >= layer$ideogram_axis_part) {
        fresh <- preserve_ideogram_aesthetics(axes[[layer$ideogram_axis_part]], layer)
        rebuilt <- reserve_chromosome_axis_space(rebuilt + fresh, list(fresh))
      }
    } else {
      rebuilt <- rebuilt + layer
    }
  }
  rebuilt$layers <- c(rebuilt$layers, generated[!used])
  rebuilt$scales <- plot$scales
  rebuilt$labels <- plot$labels
  rebuilt$guides <- plot$guides
  rebuilt$facet <- plot$facet
  rebuilt$theme <- rebuilt$theme + plot$theme
  rebuilt$data <- plot$data
  rebuilt$mapping <- plot$mapping
  if (!is.null(plot$ggnewscale_scales)) rebuilt$ggnewscale_scales <- plot$ggnewscale_scales
  reserve_chromosome_axis_space(rebuilt, list())
}
