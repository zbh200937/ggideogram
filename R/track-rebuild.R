# Rebuild source components after a named track changes.

record_ideogram_aesthetics <- function(layer) {
  layer$ideogram_aesthetic_names <- list(mapping = names(layer$mapping),
    aes_params = names(layer$aes_params), defaults = names(layer$geom$default_aes))
  layer
}

preserve_ideogram_aesthetics <- function(layer, previous) {
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
    for (arg in object$geometry_args) previous[arg] <- object[arg]
    supplied <- object$supplied_args
    if ("side" %in% supplied && !"gap" %in% supplied) {
      previous$gap <- object$gap
    }
    if ("side" %in% supplied && previous$side != "overlay" && !"offset" %in% supplied) {
      previous$offset <- 0
    }
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
      part$inherit_data <- !"data" %in% supplied
      part$inherit_mapping <- !"mapping" %in% supplied
      previous$parts <- c(previous$parts %||% list(previous), list(part))
    } else {
      for (arg in intersect(supplied, c("data", "mapping"))) previous[arg] <- object[arg]
    }
    range_changed <- any(c("limits", "value_scale", "transform", "chr") %in% supplied)
    data_changed <- !length(object$layers) && any(c("data", "mapping") %in% supplied)
    if (isTRUE(object$replace) || range_changed || data_changed) {
      previous$frozen_ranges <- NULL
      if (previous$value_scale == "global" || specs[[id]]$value_scale == "global") {
        for (key in names(specs)) {
          if (specs[[key]]$value_scale == "global") specs[[key]]$frozen_ranges <- NULL
        }
      }
    }
    specs[[id]] <- validate_track_spec(previous)
  } else {
    specs[[id]] <- object
  }
  spec <- layout$base_spec
  spec$tracks <- specs
  spec$data <- layout$data
  rebuilt <- do.call(ggideogram, spec)
  updated <- ideogram_plot_layout(rebuilt)
  updated$attachments <- layout$attachments
  metadata <- intersect(setdiff(names(layout$chrom), names(updated$chrom)),
    names(updated$data$karyotype))
  for (column in metadata) {
    updated$chrom[[column]] <- updated$data$karyotype[[column]][
      match(updated$chrom$.chr, updated$data$karyotype$.chr)]
  }
  rebuilt <- update_plot_ideogram_layout(rebuilt, updated)
  generated <- rebuilt$layers
  used <- rep(FALSE, length(generated))
  rebuilt$layers <- list()
  base <- which(vapply(generated, function(layer) isTRUE(layer$ideogram_base), logical(1)))
  base_index <- 0L
  seen <- integer()
  for (layer in plot$layers) {
    if (isTRUE(layer$ideogram_base)) {
      base_index <- base_index + 1L
      index <- base[base_index]
      rebuilt$layers <- c(rebuilt$layers, list(preserve_ideogram_aesthetics(generated[[index]], layer)))
      used[index] <- TRUE
      next
    }
    if (!is.null(layer$ideogram_scope)) {
      index <- which(!used & vapply(generated, function(candidate) {
        identical(candidate$ideogram_scope, layer$ideogram_scope) &&
          identical(candidate$ideogram_scope_slot, layer$ideogram_scope_slot)
      }, logical(1)))
      if (length(index)) {
        fresh <- generated[[index[1]]]
        replacing <- identical(layer$ideogram_scope, id) && isTRUE(object$replace) &&
          startsWith(layer$ideogram_scope_slot, "content:")
        if (!replacing) fresh <- preserve_ideogram_aesthetics(fresh, layer)
        rebuilt$layers <- c(rebuilt$layers, list(fresh))
        used[index[1]] <- TRUE
      }
      next
    }
    if (!is.null(layer$ideogram_recipe_id)) {
      if (layer$ideogram_recipe_id %in% seen) next
      seen <- c(seen, layer$ideogram_recipe_id)
      before <- length(rebuilt$layers)
      rebuilt <- rebuilt + new_ideogram_object(layer$ideogram_recipe)
      previous <- Filter(function(candidate) identical(candidate$ideogram_recipe_id,
        layer$ideogram_recipe_id), plot$layers)
      count <- min(length(previous), length(rebuilt$layers) - before)
      for (j in seq_len(count)) {
        rebuilt$layers[[before + j]] <-
          preserve_ideogram_aesthetics(rebuilt$layers[[before + j]], previous[[j]])
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
