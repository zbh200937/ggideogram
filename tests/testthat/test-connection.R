pair_kar <- data.frame(Chr = c('A', 'B', 'C'), Start = 0, End = c(100, 80, 60))

test_that('circular links do not wrap across the centre after endpoint offsets', {
  p <- ggideogram(pair_kar, orientation = 'circular', radius = 10,
    tracks = list(inner = geom_track(side = 'inner', width = 3)))
  d <- data.frame(a = 'A', b = 'B', p = 20, q = 30)
  m <- ggplot2::aes(chr1 = a, chr2 = b, position1 = p, position2 = q)
  expect_error(p + geom_chrlink(m, d, gap = c(1, 7)), 'reaches or crosses the circle center')
  expect_no_warning(ggplot2::ggplotGrob(p + geom_chrlink(m, d, gap = c(1, 2))))
})

test_that("public links use translucent interval bands by default", {
  base <- ggideogram(pair_kar[1:2, ], orientation = "horizontal")
  d <- data.frame(a = "A", b = "B", lo = 10, hi = 20, left = 30, right = 40)
  interval <- ggplot2::aes(chr1 = a, start1 = lo, end1 = hi,
    chr2 = b, start2 = left, end2 = right)
  layer <- function(plot) utils::tail(plot$layers, 1)[[1]]
  band <- layer(base + geom_chrlink(interval, d))
  expect_s3_class(band$geom, "GeomPolygon")
  expect_equal(band$aes_params$alpha, 0.35)
  expect_equal(band$aes_params$fill, ideogram_colour("synteny"))
  expect_equal(unique(band$data$.pair_start1), 10)
  expect_equal(unique(band$data$.pair_end2), 40)
  styled <- layer(base + geom_chrlink(interval, d, alpha = 1, fill = "red"))
  expect_equal(styled$aes_params$alpha, 1)
  expect_equal(styled$aes_params$fill, "red")
  mapped <- layer(base + geom_chrlink(ggplot2::aes(chr1 = a, start1 = lo, end1 = hi,
    chr2 = b, start2 = left, end2 = right, fill = a, alpha = lo), d))
  expect_null(mapped$aes_params$fill)
  expect_null(mapped$aes_params$alpha)

  nodes <- chr_nodes(data.frame(id = c("x", "y"), chr = c("A", "B"),
    start = c(10, 30), end = c(20, 40)),
    ggplot2::aes(id = id, chr = chr, start = start, end = end))
  edges <- data.frame(a = "x", b = "y")
  ids <- ggplot2::aes(from = a, to = b)
  expect_s3_class(layer(base + geom_chrlink(ids, edges, nodes = nodes))$geom, "GeomPolygon")
  points <- layer(base + geom_chrlink(ids, edges, nodes = nodes, type = "point"))
  expect_s3_class(points$geom, "GeomSegment")
  expect_equal(points$data$.pair_start1, 14.5)
  expect_equal(points$data$.pair_start2, 34.5)
})

test_that('paired points use native segments and preserve both original positions', {
  d <- data.frame(Chr1 = c('A', 'A'), Pos1 = c(10, 25), Chr2 = c('B', 'A'), Pos2 = c(30, 75))
  for (orientation in c('vertical', 'horizontal')) {
    p <- ggideogram(pair_kar, orientation = orientation, ncol = 2) +
      geom_chr_connection(data = d,
        ggplot2::aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2))
    layer <- utils::tail(p$layers, 1)[[1]]
    expect_identical(class(layer$geom)[1], 'GeomSegment')
    ends1 <- layer$data[!duplicated(layer$data$.pair_id), , drop = FALSE]
    ends2 <- layer$data[!duplicated(layer$data$.pair_id, fromLast = TRUE), , drop = FALSE]
    expect_equal(ends1$.pair_start1, d$Pos1)
    expect_equal(ends1$.pair_start2, d$Pos2)
    projected <- unproject_chr_point(p$coordinates$layout, ends1,
      '.pair_chr1', '.pair_x', '.pair_y')
    expect_equal(projected$.position, d$Pos1)
    projected <- unproject_chr_point(p$coordinates$layout, ends2,
      '.pair_chr2', '.pair_xend', '.pair_yend')
    expect_equal(projected$.position, d$Pos2)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  d$Chr2[1] <- 'C'
  p <- ggideogram(pair_kar, ncol = 2) + geom_chr_connection(data = d,
    ggplot2::aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that('synteny preserves closed widths and explicit reverse correspondence', {
  d <- data.frame(Chr1 = 'A', Start1 = c(11, 40), End1 = c(20, 40),
    Chr2 = 'B', Start2 = c(31, 50), End2 = c(40, 50), Direction = c('+', '-'))
  for (orientation in c('vertical', 'horizontal')) {
    base <- ggideogram(pair_kar, orientation = orientation)
    mapping <- ggplot2::aes(chr1 = Chr1, start1 = Start1, end1 = End1,
      chr2 = Chr2, start2 = Start2, end2 = End2, orientation = Direction, fill = Direction)
    p <- base + geom_chrlink(mapping, d)
    layer <- utils::tail(p$layers, 1)[[1]]
    expect_identical(class(layer$geom)[1], 'GeomPolygon')
    pairs <- layer$data[!duplicated(layer$data$.pair_id), , drop = FALSE]
    expect_equal(as.integer(table(layer$data$.pair_id)), c(4L, 4L))
    expect_equal(pairs$.pair_from2, c(30, 50))
    expect_equal(pairs$.pair_to2, c(40, 49))
    expect_equal(pairs$.pair_to1 - pairs$.pair_from1, c(10, 1))
    expect_equal(pairs$.pair_curvature, c(0, 0))
    expect_equal(unique(layer$data$Start1), d$Start1)
    expect_no_error(ggplot2::ggplotGrob(p))
    q <- base + geom_chrlink(mapping, d, curvature = 0.25, curve_points = 9)
    curved <- utils::tail(q$layers, 1)[[1]]$data
    expect_equal(as.integer(table(curved$.pair_id)), c(18L, 18L))
    expect_equal(unique(curved$.pair_curvature), 0.25)
    fields <- c(".pair_from1", ".pair_to1", ".pair_from2", ".pair_to2")
    expect_equal(curved[!duplicated(curved$.pair_id), fields], pairs[fields],
      ignore_attr = TRUE)
    expect_no_error(ggplot2::ggplotGrob(q))
  }
})

test_that('linear local curves keep endpoints and obey explicit straight or bend choices', {
  points <- data.frame(Chr1 = 'A', Pos1 = 25, Chr2 = 'A', Pos2 = 75)
  mapping <- ggplot2::aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2)
  for (orientation in c('vertical', 'horizontal')) for (reverse in c(FALSE, TRUE)) {
    base <- ggideogram(chr_view(pair_kar, 'A', 20, 80), orientation = orientation,
      reverse_chr = if (reverse) 'A' else character())
    p <- base + geom_chr_connection(mapping, points, curve_points = 9)
    d <- utils::tail(p$layers, 1)[[1]]$data
    expect_equal(nrow(d), 8)
    expect_equal(unique(d$.pair_start1), 25)
    expect_equal(unique(d$.pair_start2), 75)
    expect_equal(unproject_chr_point(p$coordinates$layout, d[1, ],
      '.pair_chr1', '.pair_x', '.pair_y')$.position, 25)
    expect_equal(unproject_chr_point(p$coordinates$layout, d[8, ],
      '.pair_chr2', '.pair_xend', '.pair_yend')$.position, 75)
    a <- c(d$.pair_x[1], d$.pair_y[1]); b <- c(d$.pair_xend[8], d$.pair_yend[8])
    mid <- c(d$.pair_xend[4], d$.pair_yend[4])
    expect_gt(abs((b[1] - a[1]) * (mid[2] - a[2]) -
      (b[2] - a[2]) * (mid[1] - a[1])), 0)
    expect_no_error(ggplot2::ggplotGrob(p))
    straight <- base + geom_chr_connection(mapping, points, curvature = 0)
    explicit_zero <- base + geom_chr_connection(mapping, points, bend = 0)
    expect_equal(nrow(utils::tail(straight$layers, 1)[[1]]$data), 1)
    expect_equal(nrow(utils::tail(explicit_zero$layers, 1)[[1]]$data), 1)
    curved <- base + geom_chr_connection(mapping, points, curvature = 0.4, curve_points = 9)
    expect_equal(nrow(utils::tail(curved$layers, 1)[[1]]$data), 8)
  }
})

test_that('paired window clipping uses the same correspondence fraction', {
  d <- data.frame(Chr1 = 'A', Start1 = 11, End1 = 50, Chr2 = 'A', Start2 = 51, End2 = 90)
  p <- ggideogram(chr_view(pair_kar, 'A', 20, 80)) + geom_chr_synteny(data = d,
    ggplot2::aes(chr1 = Chr1, start1 = Start1, end1 = End1, chr2 = Chr2, start2 = Start2, end2 = End2))
  x <- utils::tail(p$layers, 1)[[1]]$data
  expect_equal(unique(x$.pair_from1), 20)
  expect_equal(unique(x$.pair_to1), 40)
  expect_equal(unique(x$.pair_from2), 60)
  expect_equal(unique(x$.pair_to2), 80)
  expect_equal(unique(x$.pair_start1), 11)
  expect_equal(unique(x$.pair_end2), 90)
  p <- ggideogram(chr_view(pair_kar, 'A', 20, 80)) + geom_chr_synteny(data = d,
    ggplot2::aes(chr1 = Chr1, start1 = Start1, end1 = End1, chr2 = Chr2, start2 = Start2, end2 = End2),
    orientation = '-')
  x <- utils::tail(p$layers, 1)[[1]]$data
  expect_equal(unique(x$.pair_from2), 80)
  expect_equal(unique(x$.pair_to2), 50)
  points <- data.frame(Chr1 = c('A', 'A'), Pos1 = c(5, 30), Chr2 = c('A', 'B'), Pos2 = c(40, 40))
  base <- ggideogram(chr_view(pair_kar, 'A', 20, 80))
  p <- base + geom_chr_connection(data = points,
    ggplot2::aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2))
  expect_equal(length(p$layers), length(base$layers))
  points$Pos2[1] <- 110
  expect_error(base + geom_chr_connection(data = points,
    ggplot2::aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2)), 'source bounds')
})

test_that('same-chromosome bends retain the true outer endpoints', {
  d <- data.frame(Chr1 = 'A', Pos1 = 20, Chr2 = 'A', Pos2 = 80)
  p <- ggideogram(pair_kar[1, ], orientation = 'horizontal') + geom_chr_connection(data = d,
    ggplot2::aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2), bend = 4)
  layer <- utils::tail(p$layers, 1)[[1]]
  expect_identical(class(layer$geom)[1], 'GeomSegment')
  expect_equal(nrow(layer$data), 2)
  expect_equal(layer$data$.pair_xend[1], layer$data$.pair_x[2])
  expect_equal(layer$data$.pair_yend[1], layer$data$.pair_y[2])
  a <- unproject_chr_point(p$coordinates$layout, layer$data[1, ], '.pair_chr1', '.pair_x', '.pair_y')
  b <- unproject_chr_point(p$coordinates$layout, layer$data[2, ], '.pair_chr2', '.pair_xend', '.pair_yend')
  expect_equal(a$.position, 20)
  expect_equal(b$.position, 80)
  expect_lte(min(layer$data$.pair_y), min(p$coordinates$layout$bounds$y))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that('sampled curves draw native arrowheads only at whole-pair endpoints', {
  d <- data.frame(Chr1 = c('A', 'A'), Pos1 = c(20, 50), Chr2 = c('A', 'B'), Pos2 = c(80, 50),
    Score = c(0.4, 0.8))
  for (direction in c('horizontal', 'vertical', 'circular')) for (ends in c('first', 'last', 'both')) {
    base <- ggideogram(pair_kar, orientation = direction, reverse_chr = 'A')
    p <- base + geom_chr_connection(ggplot2::aes(chr1 = Chr1, position1 = Pos1,
      chr2 = Chr2, position2 = Pos2, alpha = Score), d, curve_points = 8,
      arrow = grid::arrow(length = grid::unit(1.5, 'mm'), ends = ends, type = 'closed'),
      colour = '#123456', linewidth = 0.3)
    layers <- p$layers[-seq_along(base$layers)]
    heads <- c(first = 0L, last = 0L)
    for (layer in layers) {
      expect_identical(class(layer$geom)[1], 'GeomSegment')
      expect_equal(layer$aes_params$colour, '#123456')
      expect_equal(layer$aes_params$linewidth, 0.3)
      expect_equal(layer$data$Score, d$Score[layer$data$.pair_id])
      arrow <- layer$geom_params$arrow
      if (is.null(arrow)) next
      if (arrow$ends %in% c(1, 3)) heads['first'] <- heads['first'] + nrow(layer$data)
      if (arrow$ends %in% c(2, 3)) heads['last'] <- heads['last'] + nrow(layer$data)
      expect_equal(arrow$length, grid::unit(1.5, 'mm'))
      expect_equal(arrow$type, 2L)
    }
    expect_equal(unname(heads), c(if (ends %in% c('first', 'both')) 2L else 0L,
      if (ends %in% c('last', 'both')) 2L else 0L))
    expected_rows <- if (direction == 'circular') 14 else 8
    expect_equal(sum(vapply(layers, function(layer) nrow(layer$data), integer(1))), expected_rows)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})
