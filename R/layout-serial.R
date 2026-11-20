# Serial end-to-end chromosome layout for genome-wide tracks and macrosynteny.

is_serial_layout <- function(layout) identical(layout$chr_layout, "serial")

serial_layout_slots <- function(k, ncol, order_by, genome_order, homolog_order) {
  n <- nrow(k)
  if (order_by == "input") {
    if (!is.null(genome_order) || !is.null(homolog_order)) {
      stopf("Explicit group orders require `order_by = 'genome'` or 'homolog'.")
    }
    if (!is.null(k$.genome) && is.null(ncol) && length(unique(chr_key(k$.genome, "group", k$.assembly))) > 1L) {
      # Multi-genome serial defaults cleanly to one row per genome in input order.
      gkey <- chr_key(k$.genome, "group", k$.assembly)
      glevels <- unique(gkey)
      row <- match(gkey, glevels) - 1L
      ord <- order(row, seq_len(n))
      row <- row[ord]
      col <- stats::ave(seq_len(n), row, FUN = seq_along) - 1L
      glabels <- vapply(glevels, function(gl) {
        idx <- which(gkey == gl)[1L]
        if (!is.null(k$.assembly) && length(unique(k$.assembly[k$.genome == k$.genome[idx]])) > 1L) {
          paste0(k$.genome[idx], " (", k$.assembly[idx], ")")
        } else {
          k$.genome[idx]
        }
      }, character(1), USE.NAMES = FALSE)
      return(list(index = ord, row = row, col = col, ncol = max(tabulate(row + 1L)),
                  row_labels = glabels))
    }
    cols <- validate_layout_ncol(ncol, n)
    row <- (seq_len(n) - 1L) %/% cols
    col <- (seq_len(n) - 1L) %% cols
    return(list(index = seq_len(n), row = row, col = col, ncol = cols, row_labels = NULL))
  }
  if (is.null(k$.genome)) stopf("Grouped layout requires a `genome` mapping.")
  genomes <- layout_level_order(k$.genome, genome_order, "genome_order")
  group_key <- chr_key(k$.genome, "group", k$.assembly)
  group_levels <- unlist(lapply(genomes, function(genome) {
    unique(group_key[k$.genome == genome])
  }), use.names = FALSE)
  group <- match(group_key, group_levels) - 1L
  has_homolog <- !is.null(k$.homolog_group)
  if (order_by == "homolog" && !has_homolog) {
    stopf("`order_by = 'homolog'` requires an explicit `homolog` mapping.")
  }
  if (!has_homolog && !is.null(homolog_order)) {
    stopf("`homolog_order` requires a `homolog` mapping.")
  }
  if (order_by == "homolog") {
    homologs <- layout_level_order(k$.homolog_group, homolog_order, "homolog_order")
    homolog <- match(k$.homolog_group, homologs) - 1L
    ord <- order(homolog, group, seq_len(n))
    row <- homolog[ord]
    col <- stats::ave(seq_len(n), row, FUN = seq_along) - 1L
    return(list(index = ord, row = row, col = col, ncol = max(tabulate(row + 1L)),
                row_labels = homologs))
  }
  # order_by == "genome": one row per genome/assembly
  ord <- if (has_homolog) {
    homologs <- layout_level_order(k$.homolog_group, homolog_order, "homolog_order")
    order(group, match(k$.homolog_group, homologs), seq_len(n))
  } else {
    order(group, seq_len(n))
  }
  row <- group[ord]
  col <- stats::ave(seq_len(n), row, FUN = seq_along) - 1L
  glabels <- vapply(group_levels, function(gl) {
    idx <- which(group_key == gl)[1L]
    if (!is.null(k$.assembly) && length(unique(k$.assembly[k$.genome == k$.genome[idx]])) > 1L) {
      paste0(k$.genome[idx], " (", k$.assembly[idx], ")")
    } else {
      k$.genome[idx]
    }
  }, character(1), USE.NAMES = FALSE)
  list(index = ord, row = row, col = col, ncol = max(tabulate(row + 1L)),
       row_labels = glabels)
}

serial_row_track_extents <- function(track_geometry, row, chr_keys, chromosome_width, row_count) {
  left_ext <- numeric(row_count)
  right_ext <- numeric(row_count)
  tbl <- track_geometry$table
  owners <- track_geometry$owners
  m_ext <- attr(track_geometry, "marker_extent") %||% c(left = 0, right = 0)
  half <- chromosome_width / 2
  for (r in seq_len(row_count) - 1L) {
    r_chrs <- chr_keys[row == r]
    l_val <- m_ext[["left"]]
    r_val <- m_ext[["right"]]
    if (nrow(tbl)) {
      for (ti in seq_len(nrow(tbl))) {
        sd <- tbl$side[ti]
        if (sd == "overlay") next
        if (any(r_chrs %in% owners[[tbl$id[ti]]])) {
          ext <- abs(tbl$high_offset[ti]) - half
          if (sd == "left") l_val <- max(l_val, ext) else r_val <- max(r_val, ext)
        }
      }
    }
    left_ext[r + 1L] <- l_val
    right_ext[r + 1L] <- r_val
  }
  list(left = left_ext, right = right_ext)
}

serial_ideogram_layout <- function(
    karyotype, k, ncol, chromosome_width, chromosome_gap, max_chr_length,
    row_gap, orientation, scale_length, curve_points, tracks,
    order_by, genome_order, homolog_order, reverse_chr, min_gap, row_align) {
  chromosome_count <- nrow(k)
  slots <- serial_layout_slots(k, ncol, order_by, genome_order, homolog_order)
  k <- k[slots$index, , drop = FALSE]
  row <- slots$row
  column <- slots$col
  per_row <- slots$ncol
  row_count <- max(row) + 1L

  track_geometry <- resolve_track_geometry(tracks, chromosome_width,
    karyotype$karyotype$.chr, k$.chr, compact = TRUE, min_gap = min_gap)
  attr(track_geometry, "marker_extent") <- attr(tracks, "marker_extent")
  axis_gaps <- vapply(tracks, function(spec) {
    if (is.list(spec$axis)) return(track_axis_gap(spec$axis$gap))
    if (isTRUE(spec$axis)) track_axis_gap() else 0
  }, numeric(1))
  row_gap <- max(row_gap, 3 * axis_gaps)

  biological_length <- k$.end - k$.start
  row_total_bp <- vapply(seq_len(row_count) - 1L, function(r) {
    sum(biological_length[row == r])
  }, numeric(1))

  if (scale_length == "global") {
    units_per_bp <- rep(max_chr_length / max(row_total_bp), chromosome_count)
    display_length <- biological_length * units_per_bp
  } else if (scale_length == "per_genome") {
    units_per_bp <- max_chr_length / row_total_bp[row + 1L]
    display_length <- biological_length * units_per_bp
  } else {
    max_in_row <- max(tabulate(row + 1L))
    display_length <- rep(max_chr_length / max_in_row, chromosome_count)
    units_per_bp <- display_length / biological_length
  }

  seg_start <- numeric(chromosome_count)
  seg_end <- numeric(chromosome_count)
  row_span <- numeric(row_count)
  for (r in seq_len(row_count) - 1L) {
    idx <- which(row == r)
    idx <- idx[order(column[idx])]
    lens <- display_length[idx]
    starts <- c(0, utils::head(cumsum(lens + chromosome_gap), -1L))
    ends <- starts + lens
    seg_start[idx] <- starts
    seg_end[idx] <- ends
    row_span[r + 1L] <- utils::tail(ends, 1L)
  }
  total_long <- max(row_span)
  row_shift <- if (identical(row_align, "center")) (total_long - row_span) / 2 else numeric(row_count)
  long_start <- seg_start + row_shift[row + 1L]
  long_end <- seg_end + row_shift[row + 1L]

  extents <- serial_row_track_extents(track_geometry, row, k$.chr, chromosome_width, row_count)
  row_left <- extents$left
  row_right <- extents$right

  if (row_count > 1L) {
    step_cross <- chromosome_width + utils::head(row_right, -1L) + utils::tail(row_left, -1L) + row_gap
    offsets <- c(0, cumsum(step_cross))
  } else {
    offsets <- 0
  }
  total_cross <- utils::tail(offsets, 1L) + chromosome_width + row_left[1L] + row_right[row_count]

  if (orientation == "horizontal") {
    row_center <- total_cross - row_left[1L] - chromosome_width / 2 - offsets
    axis_start_x <- long_start
    axis_end_x <- long_end
    axis_start_y <- row_center[row + 1L]
    axis_end_y <- row_center[row + 1L]
    x_min <- axis_start_x
    x_max <- axis_end_x
    y_min <- axis_start_y - chromosome_width / 2
    y_max <- axis_start_y + chromosome_width / 2
    bounds <- list(x = c(0, total_long), y = c(0, total_cross))
  } else {
    row_center <- row_left[1L] + chromosome_width / 2 + offsets
    axis_start_x <- row_center[row + 1L]
    axis_end_x <- row_center[row + 1L]
    axis_start_y <- total_long - long_start
    axis_end_y <- total_long - long_end
    x_min <- axis_start_x - chromosome_width / 2
    x_max <- axis_start_x + chromosome_width / 2
    y_min <- axis_end_y
    y_max <- axis_start_y
    bounds <- list(x = c(0, total_cross), y = c(0, total_long))
  }

  chrom <- data.frame(
    .chr = k$.chr,
    .start = k$.start,
    .end = k$.end,
    .centromere_start = k$.centromere_start,
    .centromere_end = k$.centromere_end,
    .length = biological_length,
    .display_length = display_length,
    .units_per_bp = units_per_bp,
    .row = row,
    .col = column,
    .axis_start_x = axis_start_x,
    .axis_start_y = axis_start_y,
    .axis_end_x = axis_end_x,
    .axis_end_y = axis_end_y,
    .x_min = pmin(x_min, x_max),
    .x_max = pmax(x_min, x_max),
    .y_min = pmin(y_min, y_max),
    .y_max = pmax(y_min, y_max),
    stringsAsFactors = FALSE
  )
  if (!is.null(karyotype$view)) {
    chrom$.source_start <- k$.source_start
    chrom$.source_end <- k$.source_end
    chrom$.cut_start <- chrom$.start > chrom$.source_start
    chrom$.cut_end <- chrom$.end < chrom$.source_end
  }
  semantic_fields <- intersect(c(".chr_name", ".genome", ".assembly",
                                 ".homolog_group", ".display_name"), names(k))
  chrom[semantic_fields] <- k[semantic_fields]
  chrom$.name_x <- (chrom$.axis_start_x + chrom$.axis_end_x) / 2
  chrom$.name_y <- (chrom$.axis_start_y + chrom$.axis_end_y) / 2

  reverse_chr <- validate_chr(reverse_chr, "reverse_chr")
  unknown <- setdiff(reverse_chr, karyotype$karyotype$.chr)
  if (length(unknown)) stopf("Unknown reverse chromosome: %s.", format_chr_rows(unknown))
  reverse <- chrom$.chr %in% reverse_chr
  chrom$.direction <- ifelse(reverse, -1, 1)
  for (axis in c("x", "y")) {
    start_f <- paste0(".axis_start_", axis)
    end_f <- paste0(".axis_end_", axis)
    original <- chrom[[start_f]][reverse]
    chrom[[start_f]][reverse] <- chrom[[end_f]][reverse]
    chrom[[end_f]][reverse] <- original
  }

  centromere_start <- project_positions_raw(chrom, k$.chr, k$.centromere_start, allow_na = TRUE)
  centromere_end <- project_positions_raw(chrom, k$.chr, k$.centromere_end, allow_na = TRUE)
  chrom$.centromere_start_x <- centromere_start$x
  chrom$.centromere_start_y <- centromere_start$y
  chrom$.centromere_end_x <- centromere_end$x
  chrom$.centromere_end_y <- centromere_end$y

  structure(
    list(
      data = karyotype,
      chrom = chrom,
      index = stats::setNames(seq_len(chromosome_count), chrom$.chr),
      ncol = per_row,
      nrow = row_count,
      row_height = rep(total_long, row_count),
      chromosome_width = chromosome_width,
      chromosome_gap = chromosome_gap,
      row_gap = row_gap,
      min_gap = min_gap,
      max_chr_length = max_chr_length,
      curve_points = curve_points,
      orientation = orientation,
      scale_length = scale_length,
      order_by = order_by,
      chr_layout = "serial",
      row_align = row_align,
      length_comparable = identical(scale_length, "global"),
      track_layout = tracks %||% track_layout(),
      tracks = track_geometry$table,
      track_owners = track_geometry$owners,
      track_ranges = initialize_track_ranges(track_geometry$table, chrom$.chr),
      track_extent = c(left = track_geometry$left, right = track_geometry$right),
      attachments = list(chr = list(), locus = list()),
      bounds = bounds,
      units = "chromosome_width",
      serial = list(
        row_labels = slots$row_labels,
        row_centers = row_center,
        row_spans = row_span,
        row_shifts = row_shift,
        row_left = row_left,
        row_right = row_right,
        row_align = row_align
      )
    ),
    class = c("ideogram_layout_v2", "ideogram_layout")
  )
}
