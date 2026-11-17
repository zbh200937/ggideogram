#' @noRd
geom_chr_transcript <- function(
    mapping = NULL, data = NULL, track = NULL, stat = "identity", position = "identity",
    ..., lane_order = NULL, block_height = 0.35,
    line_colour = ideogram_colour("annotation"), line_width = ideogram_linewidth("annotation"),
    arrow = grid::arrow(length = grid::unit(0.65, "mm"), type = "open"),
    arrow_min_bp = 0, arrow_spacing_bp = NULL, arrow_margin_bp = 0,
    labels = TRUE, label_size = ideogram_text_size("label"), label_gap = 0.3,
    label_colour = ideogram_colour("text"), label_family = NULL, label_fontface = "plain",
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
    label_colour = label_colour, label_family = label_family, label_fontface = label_fontface,
    na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes,
    ideogram_no_scope_inherit = !missing(inherit.aes) && isFALSE(inherit.aes)),
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
    label_colour = ideogram_colour("text"), label_family = NULL, label_fontface = "plain",
    na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  object <- geom_chr_transcript(mapping, data, track, stat, position, ...,
    lane_order = lane_order, block_height = block_height,
    line_colour = line_colour, line_width = line_width, arrow = arrow,
    arrow_min_bp = arrow_min_bp, arrow_spacing_bp = arrow_spacing_bp,
    arrow_margin_bp = arrow_margin_bp,
    labels = labels, label_size = label_size, label_gap = label_gap,
    label_colour = label_colour, label_family = label_family, label_fontface = label_fontface,
    na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes)
  object$id <- "gene"
  object$ideogram_no_scope_inherit <- !missing(inherit.aes) && isFALSE(inherit.aes)
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
  raw_id <- fields[[object$id]]
  data$.model_id <- validate_chr(raw_id, object$id)
  attr(data, "model_id_levels") <- if (is.factor(raw_id)) as.character(levels(raw_id)) else NULL
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

merge_intervals_idx <- function(row_idx, starts, ends) {
  n <- length(row_idx)
  if (!n) return(list(idx = integer(0), start = numeric(0), end = numeric(0)))
  if (n == 1L) return(list(idx = row_idx, start = starts, end = ends))
  ord <- order(starts, ends)
  row_idx <- row_idx[ord]
  starts <- starts[ord]
  ends <- ends[ord]
  cum_max <- cummax(ends)
  new_grp <- c(TRUE, starts[-1L] > cum_max[-n] + 1)
  first_i <- which(new_grp)
  last_i <- c(first_i[-1L] - 1L, n)
  list(idx = row_idx[first_i], start = starts[first_i], end = cum_max[last_i])
}

merge_gene_blocks <- function(data) {
  n <- nrow(data)
  if (n <= 1L) return(data)
  merged <- merge_intervals_idx(seq_len(n), data$.model_start, data$.model_end)
  out <- data[merged$idx, , drop = FALSE]
  out$.model_start <- merged$start
  out$.model_end <- merged$end
  out
}

intron_arrow_spec <- function(row_idx, starts, ends, object, g_start, g_end) {
  spans <- ends - starts
  avail <- spans - 2 * object$arrow_margin_bp
  keep <- which(spans > object$arrow_min_bp & avail >= 0)
  if (!length(keep)) {
    return(list(idx = integer(0), start = numeric(0), end = numeric(0),
                midpoint = numeric(0), tip = numeric(0)))
  }
  mids <- (starts[keep] + ends[keep]) / 2
  spacing <- object$arrow_spacing_bp
  if (is.null(spacing)) {
    tips <- mids
    vis <- tips > starts[keep] & tips < ends[keep] & tips > g_start & tips < g_end
    k <- keep[vis]
    return(list(
      idx = row_idx[k],
      start = pmax(starts[k], g_start),
      end = pmin(ends[k], g_end),
      midpoint = mids[vis],
      tip = tips[vis]
    ))
  }
  counts <- floor(avail[keep] / spacing) + 1L
  rep_k <- rep(seq_along(keep), counts)
  seq_j <- sequence(counts)
  tips <- mids[rep_k] + (seq_j - (counts[rep_k] + 1) / 2) * spacing
  orig <- keep[rep_k]
  vis <- tips > starts[orig] & tips < ends[orig] & tips > g_start & tips < g_end
  orig_v <- orig[vis]
  list(
    idx = row_idx[orig_v],
    start = pmax(starts[orig_v], g_start),
    end = pmin(ends[orig_v], g_end),
    midpoint = mids[rep_k][vis],
    tip = tips[vis]
  )
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
  id_levels <- attr(data, "model_id_levels")
  default_order <- if (!is.null(id_levels)) c(intersect(id_levels, ids), setdiff(ids, id_levels)) else ids
  if (!is.null(object$lane_order)) {
    if (anyNA(object$lane_order) || anyDuplicated(object$lane_order)) {
      stopf("`lane_order` must list each identifier at most once without missing values.")
    }
  }
  order_ids <- object$lane_order %||% default_order
  key <- interaction(match(data$.model_chr, sort(unique(data$.model_chr))),
    match(data$.model_id, sort(ids)), drop = TRUE, lex.order = TRUE)
  row_groups <- split(seq_len(nrow(data)), key)
  n_grp <- length(row_groups)

  all_chr <- data$.model_chr
  all_start <- data$.model_start
  all_end <- data$.model_end
  all_type <- data$.model_type
  all_strand <- data$.model_strand
  block_types <- c("exon", "cds", "utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr")

  vis_idx <- integer(n_grp)
  vis_start <- numeric(n_grp)
  vis_end <- numeric(n_grp)
  n_vis <- 0L

  blk_idx_list <- vector("list", n_grp)
  blk_start_list <- vector("list", n_grp)
  blk_end_list <- vector("list", n_grp)

  arr_idx_list <- vector("list", n_grp)
  arr_start_list <- vector("list", n_grp)
  arr_end_list <- vector("list", n_grp)
  arr_mid_list <- vector("list", n_grp)
  arr_tip_list <- vector("list", n_grp)

  chrom_lookup <- stats::setNames(seq_len(nrow(layout$chrom)), layout$chrom$.chr)

  for (gi in seq_len(n_grp)) {
    r_idx <- row_groups[[gi]]
    st_vec <- all_strand[r_idx]
    st1 <- st_vec[1L]
    if (any(st_vec != st1)) stopf("Each gene/transcript must have one strand.")
    ci <- unname(chrom_lookup[all_chr[r_idx[1L]]])
    if (is.na(ci)) next
    g_start <- layout$chrom$.start[ci]
    g_end <- layout$chrom$.end[ci]
    first <- min(all_start[r_idx]) - 1
    last <- max(all_end[r_idx])
    if (last <= g_start || first >= g_end) next

    n_vis <- n_vis + 1L
    vis_idx[n_vis] <- r_idx[1L]
    vis_start[n_vis] <- max(first, g_start)
    vis_end[n_vis] <- min(last, g_end)

    t_vec <- all_type[r_idx]
    ex_mask <- t_vec == "exon"
    if (any(ex_mask) && st1 %in% c("+", "-")) {
      ex_u <- merge_intervals_idx(r_idx[ex_mask], all_start[r_idx[ex_mask]], all_end[r_idx[ex_mask]])
      if (length(ex_u$idx) > 1L) {
        n_ex <- length(ex_u$idx)
        ar <- intron_arrow_spec(ex_u$idx[-n_ex], ex_u$end[-n_ex], ex_u$start[-1L] - 1,
          object, g_start, g_end)
        if (length(ar$idx)) {
          arr_idx_list[[gi]] <- ar$idx
          arr_start_list[[gi]] <- ar$start
          arr_end_list[[gi]] <- ar$end
          arr_mid_list[[gi]] <- ar$midpoint
          arr_tip_list[[gi]] <- ar$tip
        }
      }
    }

    sel_mask <- t_vec %in% block_types
    if (any(sel_mask)) {
      s_idx <- r_idx[sel_mask]
      s_start <- all_start[s_idx]
      s_end <- all_end[s_idx]
      if (object$id == "gene" && length(s_idx) > 1L) {
        s_types <- t_vec[sel_mask]
        u_types <- unique(s_types)
        m_i <- vector("list", length(u_types))
        m_s <- vector("list", length(u_types))
        m_e <- vector("list", length(u_types))
        for (ti in seq_along(u_types)) {
          tm <- s_types == u_types[ti]
          mg <- merge_intervals_idx(s_idx[tm], s_start[tm], s_end[tm])
          m_i[[ti]] <- mg$idx; m_s[[ti]] <- mg$start; m_e[[ti]] <- mg$end
        }
        s_idx <- unlist(m_i, use.names = FALSE)
        s_start <- unlist(m_s, use.names = FALSE)
        s_end <- unlist(m_e, use.names = FALSE)
      }
      c_start <- pmax(s_start - 1, g_start)
      c_end <- pmin(s_end, g_end)
      ok <- c_start < c_end
      if (any(ok)) {
        blk_idx_list[[gi]] <- s_idx[ok]
        blk_start_list[[gi]] <- c_start[ok]
        blk_end_list[[gi]] <- c_end[ok]
      }
    }
  }

  if (n_vis == 0L) {
    if (!is.null(object$lane_order) && !all(ids %in% object$lane_order)) {
      stopf("`lane_order` must include every visible gene/transcript identifier.")
    }
    return(plot)
  }

  vis_idx <- vis_idx[seq_len(n_vis)]
  backbone <- data[vis_idx, , drop = FALSE]
  backbone$.model_start <- vis_start[seq_len(n_vis)]
  backbone$.model_end <- vis_end[seq_len(n_vis)]
  if (!is.null(object$lane_order) && !all(unique(backbone$.model_id) %in% object$lane_order)) {
    stopf("`lane_order` must list each input identifier exactly once.")
  }
  check_track_chromosomes(layout, object$track, backbone$.model_chr)

  blk_all_idx <- unlist(blk_idx_list, use.names = FALSE) %||% integer(0)
  blocks <- data[blk_all_idx, , drop = FALSE]
  if (length(blk_all_idx)) {
    blocks$.model_start <- unlist(blk_start_list, use.names = FALSE)
    blocks$.model_end <- unlist(blk_end_list, use.names = FALSE)
  }

  lane_count <- max(tabulate(match(backbone$.model_chr, unique(backbone$.model_chr))))
  # Share one lane pitch across chromosomes in this track with O(1) lookup.
  ord_rank <- match(backbone$.model_id, order_ids)
  lane_slot <- stats::ave(ord_rank, backbone$.model_chr,
    FUN = function(r) match(r, sort(unique(r))))
  bb_offsets <- spec$low_offset + ((lane_slot - 0.5) / lane_count) * (spec$high_offset - spec$low_offset)
  offset_lut <- stats::setNames(bb_offsets, paste(backbone$.model_chr, backbone$.model_id, sep = "\r"))
  lane_offsets <- function(rows) {
    unname(offset_lut[paste(rows$.model_chr, rows$.model_id, sep = "\r")])
  }
  project <- function(rows, positions, offsets) {
    p <- project_positions_checked(layout, rows$.model_chr, positions)
    offset_chr_points_signed(layout, rows$.model_chr, p, offsets)
  }
  offset <- bb_offsets
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
  arr_all_idx <- unlist(arr_idx_list, use.names = FALSE) %||% integer(0)
  if (length(arr_all_idx) && !is.null(object$arrow)) {
    arrows <- data[arr_all_idx, , drop = FALSE]
    arrows$.model_start <- unlist(arr_start_list, use.names = FALSE)
    arrows$.model_end <- unlist(arr_end_list, use.names = FALSE)
    arrows$.midpoint <- unlist(arr_mid_list, use.names = FALSE)
    arrows$.tip <- unlist(arr_tip_list, use.names = FALSE)
    if (nrow(arrows)) {
      arr_off <- lane_offsets(arrows)
      origin <- ifelse(arrows$.model_strand == "+", arrows$.model_start, arrows$.model_end)
      a <- project(arrows, origin, arr_off)
      b <- project(arrows, arrows$.tip, arr_off)
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
      size = object$label_size,
      family = object$label_family %||% layout$base_spec$base_family %||% "",
      fontface = object$label_fontface, colour = object$label_colour, inherit.aes = FALSE)
    label$ideogram_gene_label <- TRUE
    label$ideogram_gene_track <- object$track
    plot <- reserve_chromosome_axis_space(plot + label, list())
  }
  plot
}
