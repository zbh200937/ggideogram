circular_kar <- data.frame(Chr = c('A', 'B', 'C'), Start = 0, End = c(100, 80, 60),
  CE_start = c(40, 30, 20), CE_end = c(50, 40, 30))

expect_circular_gap_boundaries <- function(layout) {
  gap <- layout$circular$label_gap
  radius <- layout$circular$radius
  left <- circular_xy(layout, gap$left_arc, radius)
  right <- circular_xy(layout, gap$right_arc, radius)
  expect_equal(unlist(left), radius * c(x = sin(gap$left_angle*pi/180),
    y = cos(gap$left_angle*pi/180)), tolerance = 1e-10)
  expect_equal(unlist(right), radius * c(x = sin(gap$right_angle*pi/180),
    y = cos(gap$right_angle*pi/180)), tolerance = 1e-10)
  expect_equal(gap$right_angle - gap$left_angle, gap$angle)
}

test_that('automatic circular gaps have a vertical left boundary and open to the right', {
  for (k in list(circular_kar[1, ], circular_kar)) for (clockwise in c(TRUE, FALSE)) {
    l <- ideogram_layout(k, orientation = 'circular', clockwise = clockwise)
    expect_equal(unname(l$circular$gap_angle), c(rep(2, nrow(k)-1), 20))
    expect_equal(l$circular$start_angle, if (clockwise) 20 else 0)
    expect_equal(l$circular$label_gap$angle, 20)
    expect_equal(l$circular$label_gap$left_arc,
      if (clockwise) utils::tail(l$chrom$.arc_end, 1) else 0)
    expect_equal(l$circular$label_gap$right_arc,
      if (clockwise) 0 else utils::tail(l$chrom$.arc_end, 1))
    expect_circular_gap_boundaries(l)
    expect_equal(unlist(circular_xy(l, l$circular$label_gap$left_arc, 20)),
      c(x = 0, y = 20), tolerance = 1e-10)
    expect_equal(unlist(circular_xy(l, l$circular$label_gap$right_arc, 20)),
      20*c(x = sin(20*pi/180), y = cos(20*pi/180)), tolerance = 1e-10)
    expect_equal(l$data$karyotype$.start, k$Start)
    expect_equal(l$data$karyotype$.end, k$End)
  }
})

test_that('opening angles change only the closing gap and preserve genomic projection', {
  d <- data.frame(Chr = rep(c('A', 'B', 'C'), each = 3),
    Pos = c(0, 50, 100, 0, 40, 80, 0, 30, 60), ID = letters[1:9])
  for (clockwise in c(TRUE, FALSE)) for (opening in c(0, 8, 20, 55)) {
    p <- ggideogram(circular_kar, orientation = 'circular', clockwise = clockwise,
      reverse_chr = 'B', show_names = FALSE,
      gap_angle = c(C = 12, A = 4, B = 6), opening_angle = opening) +
      geom_locus(ggplot2::aes(chr = Chr, position = Pos), d)
    q <- p + geom_track(track = 'signal', side = 'inner', width = 2,
      limits = c(0, 1))
    for (plot in list(p, q)) {
      layout <- plot$coordinates$layout
      expect_equal(unname(layout$circular$gap_angle), c(4, 6, opening))
      expect_equal(layout$circular$opening_angle, opening)
      expect_equal(layout$base_spec$opening_angle, opening)
      expect_equal(layout$circular$start_angle, if (clockwise) opening else 0)
      expect_circular_gap_boundaries(layout)
      expect_equal(unlist(circular_xy(layout, layout$circular$label_gap$left_arc, 20)),
        c(x = 0, y = 20), tolerance = 1e-10)
      projected <- project_chr_point(plot, d, 'Chr', 'Pos')
      expect_equal(projected$ID, d$ID)
      expect_equal(projected$Pos, d$Pos)
      expect_equal(unproject_chr_point(plot, projected, 'Chr', '.x', '.y')$.position,
        d$Pos, tolerance = 1e-10)
    }
    expect_equal(q$coordinates$layout$chrom$.units_per_bp,
      p$coordinates$layout$chrom$.units_per_bp)
  }
  for (opening in list(-1, Inf, NA_real_, c(10, 20), '20')) {
    expect_error(ideogram_layout(circular_kar, orientation = 'circular',
      opening_angle = opening), 'opening_angle')
  }
  expect_error(ideogram_layout(circular_kar, orientation = 'circular',
    gap_angle = 4, opening_angle = 352), 'less than 360')
})

test_that('opening overrides the final visible chromosome after ordering and selection', {
  k <- data.frame(Genome = c('G1', 'G2'), Chr = '1', Start = 0, End = c(100, 80))
  keys <- with(k, chr_key(Genome, Chr))
  semantic <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start,
    end = End, genome = Genome))
  p <- ggideogram(semantic, orientation = 'circular', order_by = 'genome',
    genome_order = c('G2', 'G1'), gap_angle = stats::setNames(c(12, 27), keys),
    opening_angle = 15)
  expect_equal(p$coordinates$layout$chrom$.chr, rev(keys))
  expect_equal(unname(p$coordinates$layout$circular$gap_angle), c(27, 15))
  selected <- ggideogram(circular_kar, chr = c('C', 'A'), orientation = 'circular',
    gap_angle = c(A = 4, B = 6, C = 12), opening_angle = 10)
  expect_equal(selected$coordinates$layout$chrom$.chr, c('C', 'A'))
  expect_equal(unname(selected$coordinates$layout$circular$gap_angle), c(12, 10))
  view <- ggideogram(chr_view(circular_kar, 'B', 20, 60), orientation = 'circular',
    opening_angle = 8) + geom_track(track = 'signal', side = 'outer', width = 2)
  expect_equal(view$coordinates$layout$circular$opening_angle, 8)
  expect_equal(view$coordinates$layout$chrom$.start, 20)
  expect_equal(view$coordinates$layout$chrom$.end, 60)
})

test_that('explicit gaps and angles preserve closing boundaries after chromosome ordering', {
  k <- data.frame(Genome = c('G1', 'G2'), Chr = '1', Start = 0, End = c(100, 80))
  keys <- with(k, chr_key(Genome, Chr))
  semantic <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome))
  gaps <- stats::setNames(c(12, 27), keys)
  for (clockwise in c(TRUE, FALSE)) {
    l <- ideogram_layout(semantic, orientation = 'circular', clockwise = clockwise,
      order_by = 'genome', genome_order = c('G2', 'G1'), gap_angle = gaps, reverse_chr = keys[1])
    expect_equal(l$chrom$.chr, rev(keys))
    expect_equal(unname(l$circular$gap_angle), c(27, 12))
    expect_equal(l$circular$start_angle, if (clockwise) 12 else 0)
    expect_equal(l$circular$label_gap$angle, 12)
    expect_circular_gap_boundaries(l)
    expect_equal(unlist(circular_xy(l, l$circular$label_gap$left_arc, 20)),
      c(x = 0, y = 20), tolerance = 1e-10)
    d <- data.frame(Key = rep(keys, each = 3), Pos = c(0, 50, 100, 0, 40, 80), ID = 1:6)
    q <- project_chr_point(l, d, 'Key', 'Pos')
    expect_equal(q$ID, d$ID)
    expect_equal(q$Pos, d$Pos)
    expect_equal(unproject_chr_point(l, q, 'Key', '.x', '.y')$.position, d$Pos)
    explicit <- ideogram_layout(semantic, orientation = 'circular', clockwise = clockwise,
      start_angle = 73, gap_angle = 7)
    expect_equal(explicit$circular$start_angle, 73)
    expect_equal(unname(explicit$circular$gap_angle), c(7, 7))
    expect_equal(explicit$circular$label_gap$angle, 7)
    expect_circular_gap_boundaries(explicit)
  }
  expect_error(ideogram_layout(circular_kar, orientation = 'circular', start_angle = Inf), 'finite angle')
  expect_error(ideogram_layout(circular_kar, orientation = 'circular', start_angle = numeric()), 'finite angle')
})

test_that('circular default gaps survive local views and whole-layout track replay', {
  for (clockwise in c(TRUE, FALSE)) {
    view <- chr_view(circular_kar, 'B', 20, 60)
    d <- data.frame(Chr = 'B', Pos = c(20, 40, 60), ID = letters[1:3])
    p <- ggideogram(view, orientation = 'circular', clockwise = clockwise,
      reverse_chr = 'B', show_names = FALSE) +
      geom_locus(ggplot2::aes(chr = Chr, position = Pos), d, size = 1.7)
    layout <- p$coordinates$layout
    expect_null(layout$base_spec$start_angle)
    expect_null(layout$base_spec$gap_angle)
    expect_equal(unname(layout$circular$gap_angle), 20)
    expect_equal(layout$data$karyotype$.end, circular_kar$End)
    q <- p + geom_track(track = 'signal', side = 'outer', width = 2, limits = c(0, 1))
    expect_equal(q$coordinates$layout$circular$label_gap, layout$circular$label_gap)
    expect_equal(q$coordinates$layout$chrom$.units_per_bp, layout$chrom$.units_per_bp)
    point <- Filter(function(layer) inherits(layer$geom, 'GeomPoint'), q$layers)[[1]]
    expect_equal(point$data$Pos, d$Pos)
    expect_equal(point$data$ID, d$ID)
    projected <- project_chr_point(q, d, 'Chr', 'Pos')
    expect_equal(unproject_chr_point(q, projected, 'Chr', '.x', '.y')$.position, d$Pos)
    expect_circular_gap_boundaries(q$coordinates$layout)
  }
})

test_that('circular sectors preserve shared bp scale, gaps and reversible projection', {
  l <- ideogram_layout(circular_kar, orientation = 'circular', radius = 20,
    gap_angle = c(C = 30, A = 10, B = 20), reverse_chr = 'B')
  expect_equal(l$chrom$.display_length / sum(l$chrom$.display_length), c(100, 80, 60) / 240)
  expect_equal(unique(l$chrom$.units_per_bp), 20 * 300 * pi / 180 / 240)
  expect_equal(unname(l$circular$gap_angle), c(10, 20, 30))
  expect_equal(l$chrom$.arc_start[2] - l$chrom$.arc_end[1], 20 * 10 * pi / 180)
  expect_equal(l$chrom$.direction, c(1, -1, 1))
  d <- data.frame(Chr = rep(c('A', 'B', 'C'), each = 3), Pos = c(0, 50, 100, 0, 40, 80, 0, 30, 60))
  q <- project_chr_point(l, d, 'Chr', 'Pos')
  q$.y <- q$.y + 2
  expect_equal(unproject_chr_point(l, q, 'Chr', '.x', '.y')$.position, d$Pos)
  expect_equal(q$Pos, d$Pos)
  expect_true(all(diff(q$.x[q$Chr == 'B']) < 0))
  ccw <- ideogram_layout(circular_kar, orientation = 'circular', clockwise = FALSE, start_angle = 0)
  cw <- ideogram_layout(circular_kar, orientation = 'circular', start_angle = 0)
  a <- circular_xy(cw, cw$chrom$.arc_start, 20)
  b <- circular_xy(ccw, ccw$chrom$.arc_start, 20)
  expect_equal(a$x, -b$x); expect_equal(a$y, b$y)
  rotated <- ideogram_layout(circular_kar, orientation = 'circular', start_angle = 90)
  expect_equal(unlist(circular_xy(rotated, 0, 20)), c(x = 20, y = 0), tolerance = 1e-12)
  equal <- ideogram_layout(circular_kar, orientation = 'circular', scale_length = 'per_chr')
  expect_equal(length(unique(equal$chrom$.display_length)), 1L)
  expect_false(equal$length_comparable)
  expect_error(ideogram_layout(circular_kar, orientation = 'circular', gap_angle = 120), 'less than 360')
  expect_error(ideogram_layout(circular_kar, orientation = 'circular', gap_angle = c(A = 1)), 'Name every')
  expect_error(ideogram_layout(circular_kar, orientation = 'circular', ncol = 2), 'one sector ring')
  expect_error(ideogram_layout(circular_kar, orientation = 'circular', radius = 1,
    tracks = track_layout(x = track('inner', width = 2))), 'accommodate')
  g <- l$chrom[1, ]
  body <- chromosome_section_polygon(g, 1, 12)
  expect_equal(range(body$y[abs(body$x - g$.axis_start_x) < 1e-12]), c(19.5, 20.5))
  waist <- project_positions_raw(l$chrom, 'A', 45)
  expect_equal(unique(body$y[abs(body$x - waist$x) < 1e-12]), 20)
})

test_that('circular tracks keep radial direction and native non-linear geom drawing', {
  d <- data.frame(Chr = rep(c('A', 'B', 'C'), each = 4),
    Pos = rep(c(10, 20, 30, 40), 3), Value = rep(1:4, 3), Low = 0, High = 5)
  tracks <- track_layout(inner = track('inner', width = 2, limits = c(0, 5)),
    outer = track('outer', width = 2, limits = c(0, 5)),
    overlay = track('overlay', width = 1, limits = c(0, 5)))
  l <- ideogram_layout(circular_kar, orientation = 'circular', tracks = tracks, reverse_chr = 'A')
  q <- data.frame(Chr = 'A', Pos = 50, Value = c(0, 5))
  inner <- project_chr_track(l, q, 'Chr', 'Pos', 'Value', 'inner')
  outer <- project_chr_track(l, q, 'Chr', 'Pos', 'Value', 'outer')
  expect_lt(diff(inner$.y), 0); expect_gt(diff(outer$.y), 0)
  expect_equal(inner$.x, outer$.x)
  expect_equal(inner$Value, q$Value)
  base <- ggideogram(circular_kar, orientation = 'circular', tracks = tracks,
    axis = TRUE, axis_n = 2, name_size = 2.8)
  layers <- list(
    geom_track_point(ggplot2::aes(chr = Chr, position = Pos, value = Value), d, track = 'inner'),
    geom_track_line(ggplot2::aes(chr = Chr, position = Pos, value = Value), d, track = 'inner'),
    geom_track_col(ggplot2::aes(chr = Chr, position = Pos, value = Value), d, track = 'outer', width = 4),
    geom_track_area(ggplot2::aes(chr = Chr, position = Pos, value = Value), d, track = 'outer'),
    geom_track_ribbon(ggplot2::aes(chr = Chr, position = Pos, ymin = Low, ymax = High), d, track = 'inner'),
    geom_track_tile(ggplot2::aes(chr = Chr, position = Pos, value = 2.5, fill = Value), d,
      track = 'overlay', width = 10, height = 5),
    geom_track_text(ggplot2::aes(chr = Chr, position = Pos, value = Value, label = Value), d, track = 'outer'))
  for (layer in layers) expect_no_error(ggplot2::ggplotGrob(base + layer))
  expect_s3_class(utils::tail((base + layers[[1]])$layers, 1)[[1]]$geom, 'GeomPoint')
  distribution <- data.frame(Chr = rep(c('A', 'B'), each = 20), Pos = 25,
    Value = rep(seq(1, 4, length.out = 20), 2))
  expect_no_error(ggplot2::ggplotGrob(base +
    geom_track_boxplot(ggplot2::aes(chr = Chr, position = Pos, value = Value), distribution,
      track = 'inner', width = 4)))
  expect_no_error(ggplot2::ggplotGrob(base +
    geom_track_violin(ggplot2::aes(chr = Chr, position = Pos, value = Value), distribution,
      track = 'outer', width = 4)))
  expect_no_error(ggplot2::ggplotGrob(base + geom_track_axis('inner', chr = 'A', breaks = c(0, 5))))
  expect_no_error(ggplot2::ggplotGrob(base + geom_chr_marker(
    ggplot2::aes(chr = Chr, position = Pos), d, side = 'outer', size = 1.7)))
  ends <- data.frame(Chr = 'A', Pos = c(0, 100))
  expect_no_error(ggplot2::ggplotGrob(base + geom_track_tile(
    ggplot2::aes(chr = Chr, position = Pos, value = 2.5), ends,
    track = 'overlay', width = 20, height = 5)))
})

test_that('circular links attach inside tracks and retain interval arcs and correspondence', {
  tracks <- track_layout(x = track('inner', width = 3, gap = 0.5))
  base <- ggideogram(circular_kar, orientation = 'circular', tracks = tracks, show_names = FALSE)
  d <- data.frame(C1 = c('A', 'B'), P1 = c(20, 30), C2 = c('C', 'B'), P2 = c(40, 60))
  p <- base + geom_chr_connection(ggplot2::aes(chr1 = C1, position1 = P1,
    chr2 = C2, position2 = P2), d, curve_points = 12, linewidth = 0.3)
  layer <- utils::tail(p$layers, 1)[[1]]
  expect_s3_class(layer$geom, 'GeomSegment')
  expect_equal(nrow(layer$data), 22)
  first <- layer$data[c(1, 12), ]
  expect_equal(sqrt(first$.pair_x^2 + first$.pair_y^2), rep(16, 2))
  expect_equal(first$P1, d$P1); expect_equal(first$P2, d$P2)
  expect_true(all(utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]$ideogram_cartesian))
  expect_no_error(ggplot2::ggplotGrob(p))
  straight <- base + geom_chr_connection(ggplot2::aes(chr1 = C1, position1 = P1,
    chr2 = C2, position2 = P2), d, curvature = 0)
  expect_equal(nrow(utils::tail(straight$layers, 1)[[1]]$data), 2)
  intervals <- data.frame(C1 = 'A', S1 = 11, E1 = 21, C2 = 'C', S2 = 31, E2 = 41)
  ribbon <- base + geom_chr_synteny(ggplot2::aes(chr1 = C1, start1 = S1, end1 = E1,
    chr2 = C2, start2 = S2, end2 = E2), intervals, orientation = '-', curve_points = 12)
  layer <- utils::tail(ribbon$layers, 1)[[1]]
  expect_s3_class(layer$geom, 'GeomPolygon')
  expect_equal(layer$data$.pair_from2[1], 41)
  expect_equal(layer$data$.pair_to2[1], 30)
  expect_equal(nrow(layer$data), 48)
  expect_equal(sqrt(layer$data$.pair_x[13:24]^2 + layer$data$.pair_y[13:24]^2), rep(16, 12))
  expect_no_error(ggplot2::ggplotGrob(ribbon))
})

test_that('circular semantic segments draw both interval and leader endpoints', {
  base <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100),
    orientation = 'circular', show_names = FALSE)
  interval <- base + geom_chr_interval(ggplot2::aes(chr = Chr, start = Start, end = End),
    data.frame(Chr = 'A', Start = 10, End = 90), placement = 'beside', side = 'inner', gap = 2)
  build <- ggplot2::ggplot_build(interval)
  data <- utils::tail(build$data, 1)[[1]]
  expect_equal(data$start, 10); expect_equal(data$end, 90)
  expect_equal(data$y, 17.5); expect_equal(data$yend, 17.5)
  expect_gt(abs(data$xend - data$x), 0)
  grob <- ggplot2::GeomSegment$draw_panel(data,
    build$layout$panel_params[[1]], build$plot$coordinates)
  expect_gt(length(unique(as.numeric(grob$x))), 2)
  expected <- build$plot$coordinates$transform(data, build$layout$panel_params[[1]])
  ends <- c(1, length(grob$x))
  expect_equal(as.numeric(grob$x)[ends], c(expected$x, expected$xend))
  expect_equal(as.numeric(grob$y)[ends], c(expected$y, expected$yend))

  markers <- data.frame(Chr = 'A', Pos = c(49, 50, 51))
  leader <- base + geom_chr_link(ggplot2::aes(chr = Chr, position = Pos), markers,
    side = 'outer', gap = 4, position = position_chr_repel(8, units = 'bp'))
  build <- ggplot2::ggplot_build(leader)
  data <- utils::tail(build$data, 1)[[1]]
  expect_equal(data$ideogram_anchor_position, markers$Pos)
  expect_equal(data$ideogram_position, c(42, 50, 58))
  expect_equal(data$y, rep(20.5, 3)); expect_equal(data$yend, rep(24.5, 3))
  expect_s3_class(utils::tail(leader$layers, 1)[[1]]$geom, 'GeomSegment')
  grob <- ggplot2::GeomSegment$draw_panel(data,
    build$layout$panel_params[[1]], build$plot$coordinates)
  expected <- build$plot$coordinates$transform(data, build$layout$panel_params[[1]])
  for (i in seq_len(nrow(data))) {
    index <- which(grob$id == i)
    ends <- index[c(1, length(index))]
    expect_equal(as.numeric(grob$x)[ends], c(expected$x[i], expected$xend[i]))
    expect_equal(as.numeric(grob$y)[ends], c(expected$y[i], expected$yend[i]))
  }
})

test_that('circular views and composite keys preserve original coordinates and gene geometry', {
  k <- data.frame(Genome = c('G1', 'G2'), Chr = 'A', Start = 0, End = c(1000, 800))
  semantic <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome))
  l <- ideogram_layout(semantic, orientation = 'circular', order_by = 'genome', genome_order = c('G2', 'G1'))
  expect_equal(l$chrom$.genome, c('G2', 'G1'))
  expect_equal(nrow(l$chrom), 2)
  view <- chr_view(circular_kar, 'A', 20, 70)
  p <- ggideogram(view, orientation = 'circular', axis = TRUE, axis_breaks = c(20, 40, 70),
    tracks = track_layout(models = track('outer', width = 2)))
  expect_equal(p$coordinates$layout$chrom$.start, 20)
  expect_equal(p$coordinates$layout$data$karyotype$.end, circular_kar$End)
  labels <- utils::tail(p$layers, 1)[[1]]$data$label
  expect_equal(labels, c('20 bp', '40 bp', '70 bp'))
  gene <- data.frame(Chr = 'A', Start = c(21, 51), End = c(35, 69), Transcript = 't', Type = 'exon', Strand = '+')
  p <- p + geom_chr_transcript(ggplot2::aes(chr = Chr, start = Start, end = End,
    transcript = Transcript, type = Type, strand = Strand), gene, track = 'models', labels = FALSE)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, 'GeomRect'), logical(1))))
  expect_no_error(ggplot2::ggplotGrob(p))
  q <- project_chr_point(p$coordinates$layout, data.frame(Chr = 'A', Pos = c(20, 40, 70)), 'Chr', 'Pos')
  expect_equal(unproject_chr_point(p$coordinates$layout, q, 'Chr', '.x', '.y')$.position, c(20, 40, 70))
})

test_that('circular insets use radial edges and true chromosome collision geometry', {
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  d <- data.frame(Chr = 'A', Pos = 12.5); d$Plot <- list(grid::rectGrob())
  base <- ggideogram(k, orientation = 'circular', show_names = FALSE,
    tracks = track_layout(inset = track('outer', width = 8)))
  layer <- geom_chr_inset(ggplot2::aes(chr = Chr, position = Pos, plot = Plot), d,
    track = 'inset', width = grid::unit(0.08, 'npc'), height = grid::unit(0.08, 'npc'))
  built <- ggplot2::ggplot_build(base + layer)
  coords <- built$plot$coordinates$transform(utils::tail(built$data, 1)[[1]], built$layout$panel_params[[1]])
  just <- inset_justification(built$plot$coordinates, coords, 'beside', NULL, NULL)
  expect_equal(just$hjust, 0); expect_equal(just$vjust, 0)
  expect_no_error(ggplot2::ggplotGrob(base + layer))
  bad <- geom_chr_inset(ggplot2::aes(chr = Chr, position = Pos, plot = Plot), d,
    width = grid::unit(0.4, 'npc'), height = grid::unit(0.4, 'npc'), hjust = 0.5, vjust = 0.5)
  expect_error(ggplot2::ggplotGrob(base + bad), 'overlaps chromosome')
  skip_if_not_installed('patchwork')
  child <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_s3_class(patchwork::patchworkGrob(patchwork::wrap_plots(base, child)), 'gtable')
  expect_s3_class(patchwork::patchworkGrob(child + patchwork::inset_element(base, 0.4, 0.4, 0.9, 0.9)), 'gtable')
})

test_that('circular interactive points and columns retain SVG and full source metadata', {
  skip_if_not_installed('ggiraph'); skip_if_not_installed('htmlwidgets')
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = c(100, 80))
  records <- data.frame(ID = c('a', 'b'), Chr = c('A', 'B'), Start = c(20, 40),
    End = c(30, 50), Pos = c(25, 45), Value = c(2, 4))
  p <- ggideogram(k, orientation = 'circular', name_size = 2.8,
    tracks = track_layout(signal = track('inner', width = 2, limits = c(0, 5)))) +
    geom_chr_track(ggplot2::aes(chr = Chr, position = Pos, value = Value, tooltip = ID, data_id = ID),
      records, track = 'signal', geom = ggiraph::geom_point_interactive, size = 2.2)
  expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, 'GeomPoint')
  p <- p + geom_chr_connection(data = data.frame(C1 = 'A', P1 = 25, C2 = 'B', P2 = 45),
    ggplot2::aes(chr1 = C1, position1 = P1, chr2 = C2, position2 = P2),
    colour = '#123456', linewidth = 0.4)
  widgets <- lapply(c(90, 180), function(mm) as_ideogram_widget(p, records, width = mm, height = mm))
  for (w in widgets) {
    expect_match(w$x$html, 'data-id=.a.')
    expect_match(w$x$html, 'data-id=.b.')
    expect_equal(w$jsHooks$render[[1]]$data$records[[2]]$end, 50)
  }
  circles <- lapply(widgets, function(w) regmatches(w$x$html, gregexpr('<circle[^>]+>', w$x$html))[[1]])
  expect_equal(length(circles[[1]]), length(circles[[2]]))
  radii <- lapply(circles, function(tags) sub('.* r=.([0-9.]+)..*', '\\1', tags))
  expect_gt(length(radii[[1]]), 0)
  expect_equal(radii[[1]], radii[[2]])
  fonts <- lapply(widgets, function(w) unique(regmatches(w$x$html, gregexpr('font-size=.[0-9.]+(pt|px).', w$x$html))[[1]]))
  expect_gt(length(fonts[[1]]), 0)
  expect_equal(fonts[[1]], fonts[[2]])
  strokes <- lapply(circles, function(tags) sub('.*stroke-width=.([0-9.]+)..*', '\\1', tags))
  expect_equal(strokes[[1]], strokes[[2]])
  links <- lapply(widgets, function(w) regmatches(w$x$html,
    gregexpr('<polyline[^>]*stroke=.#123456.[^>]*>', w$x$html))[[1]])
  expect_gt(length(links[[1]]), 0)
  widths <- lapply(links, function(tags) unique(sub('.*stroke-width=.([0-9.]+)..*', '\\1', tags)))
  expect_equal(widths[[1]], widths[[2]])
  expect_equal(p$coordinates$ratio, 1)
  expect_equal(diff(p$coordinates$limits$x), diff(p$coordinates$limits$y))
  col <- ggideogram(k, orientation = 'circular', tracks = p$coordinates$layout$track_layout) +
    geom_chr_track(ggplot2::aes(chr = Chr, position = Pos, value = Value, data_id = ID), records,
      track = 'signal', geom = ggiraph::geom_col_interactive, width = 4)
  expect_match(as_ideogram_widget(col, width = 100, height = 100)$x$html, 'data-id=.a.')
})

test_that('circular physical labels keep annular columns and true Cartesian anchors', {
  skip_if_not_installed('ggrepel')
  d <- data.frame(Chr = 'A', Pos = c(10, 30, 50, 70), Label = c('one', 'two', 'three', 'four'))
  p <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100), orientation = 'circular',
    show_names = FALSE, tracks = track_layout(labels = track('outer', width = 12, gap = 1))) +
    geom_chr_text_repel(ggplot2::aes(chr = Chr, position = Pos, label = Label), d,
      track = 'labels', size = 2.9, seed = 12)
  draw <- function(width) {
    file <- tempfile(fileext = '.pdf')
    grDevices::cairo_pdf(file, width = width / 25.4, height = width / 25.4)
    on.exit({ grDevices::dev.off(); unlink(file) })
    b <- ggplot2::ggplot_build(p)
    grid::grid.draw(ggplot2::ggplot_gtable(b))
    state <- as.list(b$layout$panel_params[[1]]$ideogram_labels)[[1]]
    list(positions = state$cache$positions, state = state)
  }
  expect_no_warning(a <- draw(150))
  expect_no_warning(b <- draw(110))
  expect_equal(a$positions$position, d$Pos)
  expect_equal(a$positions$fontsize, b$positions$fontsize)
  s <- a$state
  expected <- s$coord$transform(data.frame(x = s$data$x_orig, y = s$data$y_orig), s$panel)
  expect_equal(a$positions$anchor_x, expected$x)
  expect_equal(a$positions$anchor_y, expected$y)
})
