test_that('overlay offsets place independent lanes inside the body', {
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  specs <- track_layout(
    a = track('overlay', width = 0.4, offset = -0.3, limits = c(0, 1)),
    b = track('overlay', width = 0.4, offset = 0.3, limits = c(0, 1)))
  layout <- ideogram_layout(k, tracks = specs)
  expect_equal(layout$tracks$low_offset, c(-0.5, 0.1))
  expect_equal(layout$tracks$high_offset, c(-0.1, 0.5))
  d <- data.frame(Chr = 'A', Position = 50, Value = 0.5)
  a <- project_chr_track(layout, d, 'Chr', 'Position', 'Value', 'a')
  b <- project_chr_track(layout, d, 'Chr', 'Position', 'Value', 'b')
  expect_equal(b$.x - a$.x, 0.6)
  expect_error(track('right', offset = 0.1), 'Only an overlay')
  expect_error(ideogram_layout(k,
    tracks = track_layout(a = track('overlay', width = 0.4, offset = 0.4))), 'beyond')
})

test_that('region fills preserve annotation coordinates and use native masked rectangles', {
  k <- data.frame(Chr = 'A', Start = 0, End = 100, CE_start = 45, CE_end = 55)
  d <- data.frame(Chr = 'A', Start = c(1, 46, 100), End = c(10, 55, 100), Value = 1:3)
  specs <- track_layout(body = track('overlay', width = 1))
  for (orientation in c('vertical', 'horizontal')) {
    p <- ggideogram(k, orientation = orientation, tracks = specs) +
      geom_chr_fill(data = d, track = 'body',
        ggplot2::aes(chr = Chr, start = Start, end = End, fill = Value))
    rect <- p$layers[[length(p$layers)]]
    expect_s3_class(rect$geom, 'GeomRect')
    expect_equal(rect$data$.source_start, d$Start)
    expect_equal(rect$data$.source_end, d$End)
    long_width <- if (orientation == 'vertical') rect$data$.ymax - rect$data$.ymin else
      rect$data$.xmax - rect$data$.xmin
    expect_equal(long_width, (d$End - d$Start + 1) * 0.4)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  p <- ggideogram(chr_view(k, 'A', 5, 70), orientation = 'horizontal', tracks = specs) +
    geom_chr_fill(data = d, track = 'body',
      ggplot2::aes(chr = Chr, start = Start, end = End, fill = Value))
  rect <- p$layers[[length(p$layers)]]$data
  expect_equal(nrow(rect), 2)
  expect_equal(rect$.xmin[1], 0)
  expect_equal(rect$.source_start, c(1, 46))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that('internal points retain GeomPoint and share the offset projection', {
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  specs <- track_layout(a = track('overlay', width = 0.4, offset = 0.3, limits = c(0, 1)))
  d <- data.frame(Chr = 'A', Position = 50, Value = 0.5)
  p <- ggideogram(k, tracks = specs) + geom_track_point(data = d, track = 'a',
    ggplot2::aes(chr = Chr, position = Position, value = Value), size = 1.5)
  expect_identical(class(p$layers[[length(p$layers)]]$geom)[1], 'GeomPoint')
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that('native track labels accept exact ends after large chromosome separation', {
  k <- data.frame(Chr = as.character(1:4), Start = 0,
    End = c(248956422, 242193529, 198295559, 190214555))
  d <- transform(k, Pos = End, Value = 0.5)
  specs <- track_layout(a = track('overlay', width = 1, limits = c(0, 1)))
  p <- ggideogram(k, tracks = specs) + geom_track_text(data = d, track = 'a',
    ggplot2::aes(chr = Chr, position = Pos, value = Value), label = 'A')
  expect_no_error(ggplot2::ggplotGrob(p))
  d$Pos[2] <- d$Pos[2] + 1
  expect_error(ggideogram(k, tracks = specs) + geom_track_text(data = d, track = 'a',
    ggplot2::aes(chr = Chr, position = Pos, value = Value), label = 'A'), 'outside')
})
