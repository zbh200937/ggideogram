#' @noRd
geom_chr_transcript <- function(
    mapping = NULL, data = NULL, track = NULL, stat = "identity", position = "identity",
    ..., lane_order = NULL, block_height = 0.35,
    line_colour = ideogram_colour("annotation"), line_width = ideogram_linewidth("annotation"),
    arrow = grid::arrow(length = grid::unit(0.65, "mm"), type = "open"),
    arrow_min_bp = 0, arrow_spacing_bp = NULL, arrow_margin_bp = 0,
    labels = TRUE, label_size = ideogram_text_size("label"), label_gap = 0.3,
    na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  if (!identical(stat, "identity")) stopf("Gene structures require stat = 'identity'.")
  if (length(block_height) != 1L || !is.finite(block_height) ||
      block_height <= 0 || block_height > 1) {
    stopf("`block_height` must be greater than zero and at most one lane.")
  }
  check_nonnegative_layout(arrow_min_bp, "arrow_min_bp")
  if (!is.null(arrow_spacing_bp)) check_positive_layout(arrow_spacing_bp, "arrow_spacing_bp")
  check_nonnegative_layout(arrow_margin_bp, "arrow_margin_bp")
  check_positive_layout(label_size, "label_size")
  check_nonnegative_layout(label_gap, "label_gap")
  structure(list(mapping = mapping, data = data, track = track,
    id = "transcript", position = position, params = list(...),
    lane_order = lane_order, block_height = block_height,
    line_colour = line_colour, line_width = line_width, arrow = arrow,
    arrow_min_bp = arrow_min_bp, arrow_spacing_bp = arrow_spacing_bp,
    arrow_margin_bp = arrow_margin_bp,
    labels = labels, label_size = label_size, label_gap = label_gap,
    na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes),
    class = "ggideogram_gene_component")
}

#' @noRd
geom_chr_gene <- function(
    mapping = NULL, data = NULL, track = NULL, stat = "identity", position = "identity",
    ..., lane_order = NULL, block_height = 0.35,
    line_colour = ideogram_colour("annotation"), line_width = ideogram_linewidth("annotation"),
    arrow = grid::arrow(length = grid::unit(0.65, "mm"), type = "open"),
    arrow_min_bp = 0, arrow_spacing_bp = NULL, arrow_margin_bp = 0,
    labels = TRUE, label_size = ideogram_text_size("label"), label_gap = 0.3,
    na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  object <- geom_chr_transcript(mapping, data, track, stat, position, ...,
    lane_order = lane_order, block_height = block_height,
    line_colour = line_colour, line_width = line_width, arrow = arrow,
    arrow_min_bp = arrow_min_bp, arrow_spacing_bp = arrow_spacing_bp,
    arrow_margin_bp = arrow_margin_bp,
    labels = labels, label_size = label_size, label_gap = label_gap,
    na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes)
  object$id <- "gene"
  object
}

prepare_gene_rows <- function(object, layout) {
  data <- object$data
  if (!is.data.frame(data)) stopf("Gene structures require data-frame `data`.")
  fields <- mapped_fields(data, object$mapping,
    c("chr", "start", "end", "type", "strand", object$id), what = "mapping")
  data$.model_chr <- validate_chr(fields$chr, "chr")
  data$.model_start <- validate_coordinate(fields$start, "start")
  data$.model_end <- validate_coordinate(fields$end, "end")
  data$.model_id <- validate_chr(fields[[object$id]], object$id)
  data$.model_type <- tolower(as.character(fields$type))
  data$.model_strand <- as.character(fields$strand)
  if (anyNA(data$.model_type) || anyNA(data$.model_strand) ||
      any(!data$.model_strand %in% c("+", "-", ".", "*", "?"))) {
    stopf("Types cannot be missing; strand must be +, -, ., * or ?.")
  }
  if (any(data$.model_start < 1 | data$.model_start > data$.model_end |
      data$.model_start != floor(data$.model_start) |
      data$.model_end != floor(data$.model_end))) {
    stopf("Gene annotations require positive integer start <= end.")
  }
  source <- layout$data$karyotype
  if (is.null(source)) source <- layout$chrom
  index <- match(data$.model_chr, source$.chr)
  if (anyNA(index) || any(data$.model_start - 1 < source$.start[index] |
                         data$.model_end > source$.end[index])) {
    stopf("Gene annotations are outside the source chromosome bounds.")
  }
  data$.model_label <- if ("label" %in% names(object$mapping)) {
    recycle_semantic(rlang::eval_tidy(object$mapping$label, data), nrow(data), "label")
  } else data$.model_id
  data
}

merge_gene_blocks <- function(data) {
  if (!nrow(data)) return(data)
  data <- data[order(data$.model_start, data$.model_end), , drop = FALSE]
  out <- data[1, , drop = FALSE]
  for (i in seq_len(nrow(data))[-1]) {
    j <- nrow(out)
    if (data$.model_start[i] <= out$.model_end[j] + 1) {
      out$.model_end[j] <- max(out$.model_end[j], data$.model_end[i])
    } else out <- rbind(out, data[i, , drop = FALSE])
  }
  out
}

intron_arrow_rows <- function(gaps, object) {
  pieces <- lapply(seq_len(nrow(gaps)), function(i) {
    gap <- gaps[i, , drop = FALSE]
    span <- gap$.model_end - gap$.model_start
    available <- span - 2 * object$arrow_margin_bp
    if (span <= object$arrow_min_bp || available < 0) return(NULL)
    spacing <- object$arrow_spacing_bp
    count <- if (is.null(spacing)) 1L else floor(available / spacing) + 1L
    tips <- if (is.null(spacing)) gap$.midpoint else
      gap$.midpoint + (seq_len(count) - (count + 1) / 2) * spacing
    tips <- tips[tips > gap$.model_start & tips < gap$.model_end]
    out <- gap[rep(1L, length(tips)), , drop = FALSE]
    out$.tip <- tips
    out
  })
  do.call(rbind, pieces)
}

#' @method ggplot_add ggideogram_gene_component
#' @export
ggplot_add.ggideogram_gene_component <- function(object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  spec <- track_table_row(layout, object$track)
  if (spec$side == "overlay") {
    stopf("Gene structures require a beside track; use internal interval/tile layers for body annotation.")
  }
  data <- prepare_gene_rows(object, layout)
  if (!nrow(data)) return(plot)
  ids <- unique(data$.model_id)
  order <- object$lane_order %||% ids
  if (anyNA(order) || anyDuplicated(order) || !setequal(order, ids)) {
    stopf("`lane_order` must list each input identifier exactly once.")
  }
  key <- interaction(match(data$.model_chr, sort(unique(data$.model_chr))),
    match(data$.model_id, sort(ids)), drop = TRUE, lex.order = TRUE)
  models <- split(data, key)
  visible <- list()
  blocks <- list()
  introns <- list()
  for (model in models) {
    strands <- unique(model$.model_strand)
    if (length(strands) != 1L) stopf("Each gene/transcript must have one strand.")
    g <- layout$chrom[match(model$.model_chr[1], layout$chrom$.chr), , drop = FALSE]
    if (is.na(g$.chr) || !nrow(g)) next
    first <- min(model$.model_start) - 1
    last <- max(model$.model_end)
    if (last <= g$.start || first >= g$.end) next
    row <- model[1, , drop = FALSE]
    row$.model_start <- max(first, g$.start)
    row$.model_end <- min(last, g$.end)
    visible[[length(visible) + 1L]] <- row
    exon_union <- merge_gene_blocks(model[model$.model_type == "exon", , drop = FALSE])
    if (nrow(exon_union) > 1L && strands %in% c("+", "-")) {
      gaps <- exon_union[-nrow(exon_union), , drop = FALSE]
      gaps$.model_start <- utils::head(exon_union$.model_end, -1)
      gaps$.model_end <- utils::tail(exon_union$.model_start, -1) - 1
      gaps$.midpoint <- (gaps$.model_start + gaps$.model_end) / 2
      gaps <- intron_arrow_rows(gaps, object)
      if (!is.null(gaps)) {
        gaps <- gaps[gaps$.tip > g$.start & gaps$.tip < g$.end, , drop = FALSE]
        gaps$.model_start <- pmax(gaps$.model_start, g$.start)
        gaps$.model_end <- pmin(gaps$.model_end, g$.end)
        introns[[length(introns) + 1L]] <- gaps
      }
    }
    selected <- model[model$.model_type %in%
      c("exon", "cds", "utr", "five_prime_utr", "three_prime_utr"), , drop = FALSE]
    if (object$id == "gene" && nrow(selected)) {
      selected <- do.call(rbind, lapply(split(selected, selected$.model_type), merge_gene_blocks))
    }
    selected$.model_start <- pmax(selected$.model_start - 1, g$.start)
    selected$.model_end <- pmin(selected$.model_end, g$.end)
    blocks[[length(blocks) + 1L]] <- selected[selected$.model_start < selected$.model_end, , drop = FALSE]
  }
  if (!length(visible)) return(plot)
  backbone <- do.call(rbind, visible)
  blocks <- do.call(rbind, blocks)
  lane_count <- max(table(backbone$.model_chr))
  # Share one lane pitch across chromosomes in this track.
  lane_offsets <- function(rows) {
    vapply(seq_len(nrow(rows)), function(i) {
      chr_ids <- order[order %in% backbone$.model_id[backbone$.model_chr == rows$.model_chr[i]]]
      fraction <- (match(rows$.model_id[i], chr_ids) - 0.5) / lane_count
      spec$low_offset + fraction * (spec$high_offset - spec$low_offset)
    }, numeric(1))
  }
  project <- function(rows, positions, offsets) {
    p <- project_positions_checked(layout, rows$.model_chr, positions)
    offset_chr_points_signed(layout, rows$.model_chr, p, offsets)
  }
  offset <- lane_offsets(backbone)
  lo <- project(backbone, backbone$.model_start, offset)
  hi <- project(backbone, backbone$.model_end, offset)
  negative <- backbone$.model_strand == "-"
  backbone$.gx <- ifelse(negative, hi$x, lo$x)
  backbone$.gy <- ifelse(negative, hi$y, lo$y)
  backbone$.gxend <- ifelse(negative, lo$x, hi$x)
  backbone$.gyend <- ifelse(negative, lo$y, hi$y)
  segment_mapping <- ggplot2::aes(x = .data$.gx, y = .data$.gy,
    xend = .data$.gxend, yend = .data$.gyend)
  plot <- plot + ggplot2::geom_segment(data = backbone, mapping = segment_mapping,
    colour = object$line_colour, linewidth = object$line_width, inherit.aes = FALSE)
  if (length(introns) && !is.null(object$arrow)) {
    arrows <- do.call(rbind, introns)
    if (nrow(arrows)) {
      origin <- ifelse(arrows$.model_strand == "+", arrows$.model_start, arrows$.model_end)
      a <- project(arrows, origin, lane_offsets(arrows))
      b <- project(arrows, arrows$.tip, lane_offsets(arrows))
      arrows$.gx <- a$x; arrows$.gy <- a$y
      arrows$.gxend <- b$x; arrows$.gyend <- b$y
      plot <- plot + ggplot2::geom_segment(data = arrows, mapping = segment_mapping,
        colour = object$line_colour, linewidth = object$line_width,
        arrow = object$arrow, inherit.aes = FALSE)
    }
  }
  if (nrow(blocks)) {
    blocks <- blocks[order(blocks$.model_type == "cds"), , drop = FALSE]
    half <- object$block_height * abs(spec$high_offset - spec$low_offset) / lane_count / 2
    offset <- lane_offsets(blocks)
    a <- project(blocks, blocks$.model_start, offset - half)
    b <- project(blocks, blocks$.model_end, offset + half)
    blocks$.gxmin <- pmin(a$x, b$x); blocks$.gxmax <- pmax(a$x, b$x)
    blocks$.gymin <- pmin(a$y, b$y); blocks$.gymax <- pmax(a$y, b$y)
    mapping <- track_mapping_without(object$mapping,
      c("chr", "start", "end", "type", "strand", "gene", "transcript", "label"))
    mapping <- combine_track_mapping(mapping, ggplot2::aes(
      xmin = .data$.gxmin, xmax = .data$.gxmax, ymin = .data$.gymin, ymax = .data$.gymax))
    for (outline in c(FALSE, TRUE)) {
      rows <- blocks[(blocks$.model_type == "exon") == outline, , drop = FALSE]
      if (!nrow(rows)) next
      rect_mapping <- if (outline) track_mapping_without(mapping, "fill") else mapping
      params <- utils::modifyList(list(fill = ideogram_colour("gene"), colour = NA), object$params)
      if (outline) {
        params$fill <- NA
        if (!"colour" %in% names(object$params)) params$colour <- object$line_colour
        if (!"linewidth" %in% names(object$params)) params$linewidth <- object$line_width
      } else if ("fill" %in% names(rect_mapping) && !"fill" %in% names(object$params)) {
        params$fill <- NULL
      }
      if ("colour" %in% names(rect_mapping) && !"colour" %in% names(object$params)) params$colour <- NULL
      plot <- plot + do.call(ggplot2::geom_rect, c(list(data = rows, mapping = rect_mapping,
        position = object$position, na.rm = object$na.rm,
        show.legend = object$show.legend, inherit.aes = object$inherit.aes), params))
    }
  }
  if (isTRUE(object$labels)) {
    g <- layout$chrom[match(backbone$.model_chr, layout$chrom$.chr), ]
    backbone$.label_position <- ifelse(g$.direction < 0, g$.end, g$.start)
    anchor <- project(backbone, backbone$.label_position, lane_offsets(backbone))
    backbone$.lx <- anchor$x; backbone$.ly <- anchor$y
    vertical <- layout$orientation == "vertical"
    angle <- rep(0, nrow(backbone))
    hjust <- if (vertical) 0.5 else 1 + object$label_gap
    if (is_circular_layout(layout)) {
      angle <- (-circular_theta(layout, anchor$x) * 180 / pi +
        if (layout$circular$clockwise) 0 else 180) %% 360
      angle <- (angle + 180) %% 360 - 180
      flip <- abs(angle) > 90
      angle[flip] <- (angle[flip] + 360) %% 360 - 180
      hjust <- ifelse(flip, -object$label_gap, 1 + object$label_gap)
    }
    label <- ggplot2::geom_text(data = backbone,
      mapping = ggplot2::aes(x = .data$.lx, y = .data$.ly, label = .data$.model_label),
      angle = angle, hjust = hjust,
      vjust = if (vertical) -object$label_gap else 0.5,
      size = object$label_size, family = layout$base_spec$base_family %||% "",
      colour = "black", inherit.aes = FALSE)
    label$ideogram_gene_label <- TRUE
    plot <- reserve_chromosome_axis_space(plot + label, list())
  }
  plot
}
