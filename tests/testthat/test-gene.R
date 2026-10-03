gene_fixture <- function() data.frame(
  Chr = '1', Start = c(11, 11, 41, 16, 46, 21, 21, 61),
  End = c(80, 30, 80, 25, 70, 90, 35, 90),
  Type = c('mRNA', 'exon', 'exon', 'CDS', 'CDS', 'mRNA', 'exon', 'exon'),
  Strand = c(rep('+', 5), rep('-', 3)),
  Tx = c(rep('a', 5), rep('b', 3)), Gene = c(rep('g1', 5), rep('g2', 3)))

gene_plot <- function(data = gene_fixture(), orientation = 'horizontal', view = NULL, ...) {
  k <- data.frame(Chr = '1', Start = 0, End = 100)
  ggideogram(if (is.null(view)) k else chr_view(k, '1', view[1], view[2]),
    orientation = orientation, tracks = track_layout(genes = track(width = 4))) +
    geom_chr_transcript(data = data, track = 'genes',
      ggplot2::aes(chr = Chr, start = Start, end = End, type = Type,
                   strand = Strand, transcript = Tx, fill = Tx), ...)
}

test_that('structures retain native geoms, closed interval widths and strand direction', {
  p <- gene_plot()
  rect <- Filter(function(l) inherits(l$geom, 'GeomRect'), p$layers)
  expect_equal(nrow(rect[[1]]$data), 2)
  expect_equal(nrow(rect[[2]]$data), 4)
  blocks <- do.call(rbind, lapply(rect, function(l) l$data))
  expect_equal(blocks$.gxmax - blocks$.gxmin,
    (blocks$.model_end - blocks$.model_start) * 0.4)
  expect_equal(blocks$.gymax - blocks$.gymin, rep(0.35 * 4 / 2, nrow(blocks)))
  index <- which(vapply(p$layers, function(l) inherits(l$geom, 'GeomRect'), logical(1)))
  built <- ggplot2::ggplot_build(p)$data
  expect_setequal(built[[index[1]]]$fill, scales::hue_pal()(1))
  expect_true(all(is.na(built[[index[2]]]$fill)))
  seg <- Filter(function(l) inherits(l$geom, 'GeomSegment'), p$layers)[[1]]$data
  expect_true(seg$.gxend[1] > seg$.gx[1])
  expect_true(seg$.gxend[2] < seg$.gx[2])
  expect_s3_class(ggplot2::ggplotGrob(p), 'gtable')
  expect_s3_class(ggplot2::ggplotGrob(gene_plot(orientation = 'vertical')), 'gtable')
})

test_that('window clips structures without inventing transcription ends', {
  p <- gene_plot(view = c(20, 60))
  rect <- do.call(rbind, lapply(Filter(function(l) inherits(l$geom, 'GeomRect'), p$layers), function(l) l$data))
  expect_true(all(rect$.model_start >= 20 & rect$.model_end <= 60))
  seg <- Filter(function(l) inherits(l$geom, 'GeomSegment'), p$layers)
  expect_equal(nrow(seg[[1]]$data), 2)
  expect_null(seg[[1]]$geom_params$arrow)
  expect_setequal(seg[[2]]$data$.model_id, c('a', 'b'))
  expect_equal(seg[[2]]$data$.midpoint, c(35, 47.5))
  expect_s3_class(ggplot2::ggplotGrob(p), 'gtable')
})

test_that('gene mode merges overlapping exons and validates source bounds', {
  d <- gene_fixture(); d$Strand <- '+'; d$Gene <- 'g'
  p <- ggideogram(data.frame(Chr = '1', Start = 0, End = 100),
    tracks = track_layout(genes = track(width = 4))) +
    geom_chr_gene(data = d, track = 'genes', ggplot2::aes(chr = Chr,
      start = Start, end = End, type = Type, strand = Strand, gene = Gene))
  rect <- do.call(rbind, lapply(Filter(function(l) inherits(l$geom, 'GeomRect'), p$layers), function(l) l$data))
  expect_equal(sum(rect$.model_type == 'exon'), 2)
  expect_equal(rect$.model_end[rect$.model_type == 'exon'], c(35, 90))
  d$End[1] <- 110
  expect_error(gene_plot(d), 'source chromosome bounds')
  one <- gene_fixture()[2, ]; one$End <- one$Start
  p <- gene_plot(one)
  rect <- Filter(function(l) inherits(l$geom, 'GeomRect'), p$layers)[[1]]$data
  expect_equal(rect$.gxmax - rect$.gxmin, 0.4)
})

test_that('UTR fills retain their annotation and exon outlines retain full extents', {
  d <- gene_fixture()
  utr <- d[2, ]; utr$End <- 15; utr$Type <- 'five_prime_UTR'
  p <- gene_plot(rbind(d, utr))
  rect <- Filter(function(l) inherits(l$geom, 'GeomRect'), p$layers)
  expect_setequal(rect[[1]]$data$.model_type, c('cds', 'five_prime_utr'))
  expect_true(all(rect[[2]]$data$.model_type == 'exon'))
  expect_equal(rect[[2]]$data$.model_start, c(10, 40, 20, 60))
  expect_equal(rect[[2]]$data$.model_end, c(30, 80, 35, 90))
  expect_false('fill' %in% names(rect[[2]]$mapping))
})

test_that('repeated intron arrows retain spacing, strand and original window positions', {
  p <- gene_plot(arrow_spacing_bp = 3, arrow_margin_bp = 1)
  arrows <- Filter(function(l) inherits(l$geom, 'GeomSegment'), p$layers)[[2]]$data
  expect_equal(arrows$.tip[arrows$.model_id == 'a'], c(32, 35, 38))
  expect_equal(diff(arrows$.tip[arrows$.model_id == 'b']), rep(3, 7))
  expect_true(all(arrows$.gxend[arrows$.model_strand == '+'] > arrows$.gx[arrows$.model_strand == '+']))
  expect_true(all(arrows$.gxend[arrows$.model_strand == '-'] < arrows$.gx[arrows$.model_strand == '-']))
  window <- gene_plot(view = c(40, 55), arrow_spacing_bp = 3, arrow_margin_bp = 1)
  clipped <- Filter(function(l) inherits(l$geom, 'GeomSegment'), window$layers)[[2]]$data
  expect_equal(clipped$.tip, arrows$.tip[arrows$.tip > 40 & arrows$.tip < 55])
})

test_that('structure thickness is shared across chromosomes with different model counts', {
  d <- rbind(gene_fixture(), transform(gene_fixture()[1:5, ], Chr = '2'))
  k <- data.frame(Chr = c('1', '2'), Start = 0, End = 100)
  for (direction in c('horizontal', 'vertical')) {
    p <- ggideogram(k, orientation = direction,
      tracks = track_layout(genes = track(width = 4))) +
      geom_chr_transcript(data = d, track = 'genes', ggplot2::aes(chr = Chr,
        start = Start, end = End, type = Type, strand = Strand, transcript = Tx))
    blocks <- do.call(rbind, lapply(Filter(function(l) inherits(l$geom, 'GeomRect'), p$layers), function(l) l$data))
    thickness <- if (direction == 'horizontal') blocks$.gymax - blocks$.gymin else blocks$.gxmax - blocks$.gxmin
    expect_setequal(blocks$.model_chr, c('1', '2'))
    expect_equal(thickness, rep(0.7, nrow(blocks)))
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that('lane labels retain their display-start clearance on reversed chromosomes', {
  d <- rbind(gene_fixture(), transform(gene_fixture()[1:5, ], Chr = '2'))
  k <- data.frame(Chr = c('1', '2'), Start = 0, End = 100)
  for (orientation in c('horizontal', 'vertical', 'circular')) {
    p <- ggideogram(k, orientation = orientation, show_names = FALSE,
      reverse_chr = '1', tracks = track_layout(genes = track(width = 4))) +
      geom_chr_transcript(data = d, track = 'genes', ggplot2::aes(chr = Chr,
        start = Start, end = End, type = Type, strand = Strand, transcript = Tx))
    index <- which(vapply(p$layers, function(l) inherits(l$geom, 'GeomText'), logical(1)))
    expect_length(index, 1)
    label <- p$layers[[index]]$data
    g <- p$coordinates$layout$chrom[match(label$.model_chr, p$coordinates$layout$chrom$.chr), ]
    expect_equal(label$.label_position, ifelse(label$.model_chr == '1', 100, 0))
    if (orientation == 'horizontal') expect_equal(label$.lx, g$.x_min)
    if (orientation == 'vertical') expect_equal(label$.ly, g$.y_max)
    if (orientation == 'circular') expect_equal(label$.lx, g$.arc_start)
    expect_equal(label$.model_strand, c('+', '-', '+'))
    built <- ggplot2::ggplot_build(p)$data[[index]]
    expect_equal(built$size, rep(3, nrow(label)))
    if (orientation == 'horizontal') expect_equal(built$hjust, rep(1.3, nrow(label)))
    if (orientation == 'vertical') expect_equal(built$vjust, rep(-0.3, nrow(label)))
    expect_identical(class(p$layers[[index]]$geom)[1], 'GeomText')
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that('circular lane labels precede the displayed sector in both directions', {
  k <- data.frame(Chr = '1', Start = 0, End = 100)
  for (clockwise in c(TRUE, FALSE)) for (start_angle in c(0, 90, 180, 270)) {
    p <- ggideogram(chr_view(k, '1', 20, 60), orientation = 'circular',
      clockwise = clockwise, start_angle = start_angle, reverse_chr = '1',
      show_names = FALSE, tracks = track_layout(genes = track(width = 4))) +
      geom_chr_transcript(data = gene_fixture(), track = 'genes', ggplot2::aes(chr = Chr,
        start = Start, end = End, type = Type, strand = Strand, transcript = Tx))
    index <- which(vapply(p$layers, function(l) inherits(l$geom, 'GeomText'), logical(1)))
    label <- p$layers[[index]]$data
    expect_equal(label$.label_position, rep(60, nrow(label)))
    expect_equal(label$.lx, rep(p$coordinates$layout$chrom$.arc_start, nrow(label)))
    built <- ggplot2::ggplot_build(p)$data[[index]]
    theta <- start_angle * pi / 180
    angle <- built$angle * pi / 180
    extend <- ifelse(built$hjust > 1, -1, 1)
    direction <- if (clockwise) 1 else -1
    along <- extend * direction * (cos(angle) * cos(theta) - sin(angle) * sin(theta))
    expect_equal(along, rep(-1, nrow(label)))
    expect_true(all(abs(built$angle) <= 90))
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})
