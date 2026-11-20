# Local horizontal comparisons share ordering and annotation-side decisions.
is_microsynteny_layout <- function(orientation, rows, chromosomes, edges) {
  identical(orientation, "horizontal") && rows == 1L &&
    chromosomes >= 2L && !is.null(edges) && nrow(edges) > 0
}

is_microsynteny <- function(layout) {
  is_microsynteny_layout(layout$orientation, layout$nrow, nrow(layout$chrom), layout$data$pair_edges)
}

microsynteny_order <- function(data, spec) {
  edges <- data$pair_edges
  keys <- data$chr_selection %||% data$karyotype$.chr
  if (spec$orientation != "horizontal" || length(keys) != 3L ||
      spec$link_layout != "auto" || !is.null(spec$chr) || spec$order_by != "input" ||
      (!is.null(spec$ncol) && spec$ncol < 3L) || is.null(edges)) return(data)
  edges <- edges[edges$a %in% keys & edges$b %in% keys, , drop = FALSE]
  degree <- table(factor(c(edges$a, edges$b), levels = keys))
  hub <- keys[degree == 2L]
  if (nrow(edges) == 2L && length(hub) == 1L) {
    leaves <- keys[keys != hub]
    data$chr_selection <- c(leaves[1], hub, leaves[2])
  }
  data
}

microsynteny_side <- function(layout, chr, kind = "locus") {
  if (!isTRUE(layout$microsynteny)) {
    return(rep(if (kind == "axis" && !is_circular_layout(layout)) "left" else "right", length(chr)))
  }
  index <- match(chr, layout$chrom$.chr)
  if (kind == "axis") ifelse(index == nrow(layout$chrom), "right", "left") else
    ifelse(index == 1L, "left", "right")
}

register_microsynteny_pairs <- function(plot, component) {
  layout <- ideogram_plot_layout(plot)
  if (is_serial_layout(layout) || layout$orientation != "horizontal" || layout$nrow != 1L ||
      nrow(layout$chrom) < 2L) return(plot)
  comp <- resolve_chr_pair_nodes(component, layout)
  d <- comp$data
  if (!is.data.frame(d) || !nrow(d)) return(plot)
  # Fast chromosome-pair discovery before full coordinate validation/clipping
  f <- mapped_fields(d, comp$mapping, c("chr1", "chr2"), what = "mapping")
  a <- validate_chr(f$chr1, "chr1")
  b <- validate_chr(f$chr2, "chr2")
  visible_chr <- layout$chrom$.chr
  keep <- (a %in% visible_chr) & (b %in% visible_chr) & (a != b)
  if (!any(keep)) return(plot)
  edges <- unique(data.frame(a = pmin(a[keep], b[keep]), b = pmax(a[keep], b[keep]),
    stringsAsFactors = FALSE))
  edges <- unique(rbind(layout$data$pair_edges, edges))
  rownames(edges) <- NULL
  if (!nrow(edges) || identical(edges, layout$data$pair_edges)) return(plot)
  # Confirm at least one cross-chromosome link survives source & view bounds
  pairs <- prepare_chr_pairs(comp, layout)
  a <- pairs$.pair_chr1; b <- pairs$.pair_chr2
  edges <- unique(data.frame(a = pmin(a, b), b = pmax(a, b), stringsAsFactors = FALSE))
  edges <- edges[edges$a != edges$b, , drop = FALSE]
  edges <- unique(rbind(layout$data$pair_edges, edges))
  rownames(edges) <- NULL
  if (!nrow(edges) || identical(edges, layout$data$pair_edges)) return(plot)
  layout$data$pair_edges <- edges
  plot <- update_plot_ideogram_layout(plot, layout)
  spec <- layout$base_spec
  # Automatically declared annotation bands are recreated from source recipes.
  prev_ext <- attr(spec$tracks, "marker_extent") %||% layout$marker_extent
  prev_ext_chr <- attr(spec$tracks, "marker_extent_by_chr") %||% layout$marker_extent_by_chr
  filtered <- Filter(function(x) is.null(x$.labels_source), spec$tracks)
  spec$tracks <- do.call(track_layout, filtered %||% list())
  attr(spec$tracks, "marker_extent") <- prev_ext
  attr(spec$tracks, "marker_extent_by_chr") <- prev_ext_chr
  rebuild_ideogram(plot, spec)
}

add_auto_locus <- function(component, plot) {
  layout <- ideogram_plot_layout(plot)
  data <- component$data
  if (!isTRUE(layout$microsynteny) || !is.data.frame(data) || !nrow(data)) {
    return(ggplot2::ggplot_add(component, plot))
  }
  if (is.null(component$mapping$position)) return(ggplot2::ggplot_add(component, plot))
  data <- scope_point_data(data, component$mapping, layout)
  if (!nrow(data)) return(plot)
  chr <- mapped_fields(data, component$mapping, "chr", what = "mapping")$chr
  sides <- microsynteny_side(layout, chr)
  text <- inherits(component, "ggideogram_labels_component") ||
    inherits(component, "ggideogram_repel_component") ||
    (inherits(component, "Layer") && inherits(component$geom, "GeomText"))
  middle <- nrow(layout$chrom) >= 3L &
    !chr %in% c(layout$chrom$.chr[1], utils::tail(layout$chrom$.chr, 1))
  if (text && any(middle)) {
    mid_data <- data[middle, , drop = FALSE]
    lbl <- as.character(mapped_fields(mid_data, component$mapping, "label", what = "mapping")$label)
    valid_mid <- !is.na(lbl) & nzchar(lbl)
    if (any(valid_mid)) {
      mapping <- track_mapping_without(component$mapping,
        setdiff(names(component$mapping), c("chr", "position", ggplot2::GeomPoint$aesthetics())))
      style <- component$aes_params %||% component$params
      style <- style[intersect(names(style), c("colour", "fill", "alpha", "stroke"))]
      defaults <- list(shape = 23, size = 2, fill = "white", colour = ideogram_colour("text"))
      defaults <- defaults[!names(defaults) %in% names(mapping)]
      plot <- plot + do.call(geom_chr_marker, c(list(mapping = mapping,
        data = mid_data[valid_mid, , drop = FALSE], gap = 0, show.legend = component$show.legend),
        utils::modifyList(defaults, style)))
    }
    data <- data[!middle, , drop = FALSE]
    sides <- sides[!middle]
    if (!nrow(data)) return(plot)
  }
  for (side in unique(sides)) {
    part <- if (inherits(component, "Layer")) ggplot2::ggproto(NULL, component) else component
    part$data <- data[sides == side, , drop = FALSE]
    if (inherits(part, "Layer")) {
      part$stat_params$side <- side
      if (!is.null(part$ideogram_marker_spec)) part$ideogram_marker_spec$side <- side
    } else part$side <- side
    plot <- ggplot2::ggplot_add(part, plot)
  }
  plot
}

microsynteny_link_gap <- function(layout, chr, sign) {
  sides <- ifelse(sign < 0, "left", "right")
  gap <- numeric(length(chr))
  for (side in unique(sides)) {
    rows <- sides == side
    keys <- unique(chr[rows])
    gap[rows] <- annotation_extent(layout, keys, side)[match(chr[rows], keys)]
  }
  axis <- layout$base_spec$axis
  selected <- if (isTRUE(axis)) layout$chrom$.chr else if (identical(axis, FALSE)) character() else axis
  axis_side <- layout$base_spec$axis_side
  if (shared_microsynteny_axis(layout)) {
    selected <- utils::tail(layout$chrom$.chr, 1)
    axis_side <- "right"
  } else if (axis_side == "auto") axis_side <- microsynteny_side(layout, chr, "axis")
  has_axis <- chr %in% selected & sides == axis_side
  gap[has_axis] <- gap[has_axis] + layout$base_spec$axis_gap + layout$base_spec$axis_width
  gap
}

annotation_extent <- function(layout, chr, side) {
  side <- rep(side, length.out = length(chr))
  if (!isTRUE(layout$microsynteny)) return(unname(layout$track_extent[side]))
  vapply(seq_along(chr), function(i) {
    rows <- layout$tracks[layout$tracks$side == side[i], , drop = FALSE]
    keep <- vapply(rows$id, function(id) chr[i] %in% track_chromosomes(layout, id), logical(1))
    m_map <- (layout$marker_extent_by_chr %||% list())[[side[i]]]
    m_ext <- if (!is.null(m_map)) {
      if (chr[i] %in% names(m_map)) unname(m_map[[chr[i]]]) else 0
    } else {
      (layout$marker_extent %||% c(left = 0, right = 0))[[side[i]]]
    }
    max(c(m_ext, abs(rows$high_offset[keep]) - layout$chromosome_width / 2))
  }, numeric(1))
}

pack_microsynteny <- function(layout) {
  keys <- layout$chrom$.chr
  above <- microsynteny_link_gap(layout, keys, rep(-1, length(keys)))
  below <- microsynteny_link_gap(layout, keys, rep(1, length(keys)))
  center <- c(0, -cumsum(layout$chromosome_width + layout$chromosome_gap +
    utils::head(below, -1) + utils::tail(above, -1)))
  shift <- center - layout$chrom$.axis_start_y
  fields <- grep("(^\\.y_|_y$)", names(layout$chrom), value = TRUE)
  for (field in fields) layout$chrom[[field]] <- layout$chrom[[field]] + shift
  outer_above <- annotation_extent(layout, keys[1], "left")
  outer_below <- annotation_extent(layout, utils::tail(keys, 1), "right")
  layout$bounds$y <- c(utils::tail(center, 1) - outer_below - layout$chromosome_width / 2,
    center[1] + outer_above + layout$chromosome_width / 2)
  layout
}

can_share_microsynteny_axis <- function(layout) {
  isTRUE(layout$microsynteny) && nrow(layout$chrom) >= 3L &&
    identical(layout$scale_length, "global")
}

shared_microsynteny_axis <- function(layout) {
  can_share_microsynteny_axis(layout) &&
    isTRUE(layout$base_spec$axis) &&
    identical(layout$base_spec$axis_side, "auto") && is.null(layout$base_spec$axis_breaks)
}
