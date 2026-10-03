test_that('interactive points retain native geoms and complete source intervals', {
  skip_if_not_installed('ggiraph')
  skip_if_not_installed('htmlwidgets')
  k <- data.frame(Chr = 'A', Start = 0, End = 1000)
  records <- data.frame(ID = c('a', 'b'), Chr = 'A', Start = c(100, 750),
    End = c(250, 900), Pos = c(230, 770), Category = c('+', '-'),
    URL = c('https://example.org/a', 'https://example.org/b'))
  p <- ggideogram(chr_view(k, 'A', 200, 800), orientation = 'horizontal',
    tracks = track_layout(body = track('overlay', width = 1, limits = c(0, 1)))) +
    geom_chr_track(data = records, track = 'body', geom = ggiraph::geom_point_interactive,
      ggplot2::aes(chr = Chr, position = Pos, value = 0.5, tooltip = ID, data_id = ID), size = 2)
  expect_true(inherits(utils::tail(p$layers, 1)[[1]]$geom, 'GeomPoint'))
  w <- as_ideogram_widget(p, records, category = 'Category', url = 'URL', width = 120, height = 50)
  expect_s3_class(w, 'girafe')
  expect_match(w$x$html, 'data-id=.a.')
  expect_match(w$x$html, 'data-id=.b.')
  expect_false(w$x$settings$sizing$rescale)
  metadata <- w$jsHooks$render[[1]]$data
  expect_equal(metadata$width, 120)
  expect_equal(metadata$records[[1]]$start, 100)
  expect_equal(metadata$records[[2]]$end, 900)
  expect_equal(metadata$fields, names(records))
  expect_equal(metadata$records[[1]]$values[[3]], records$Start[1])
  expect_equal(metadata$records[[2]]$url, records$URL[2])
  expect_error(as_ideogram_widget(p, transform(records, ID = 'a')), 'unique')
  expect_error(as_ideogram_widget(p, transform(records, Chr = 'B')), 'unknown source')
  expect_error(as_ideogram_widget(p, transform(records, End = 1100)), 'source chromosome bounds')
  expect_error(as_ideogram_widget(p, records, url = 'ID'), 'HTTP')
})

test_that('native interactive columns remain compatible with both orientations', {
  skip_if_not_installed('ggiraph')
  skip_if_not_installed('htmlwidgets')
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  d <- data.frame(ID = c('a', 'b'), Chr = 'A', Pos = c(20, 80), Value = c(10, 30))
  for (orientation in c('vertical', 'horizontal')) {
    p <- ggideogram(k, orientation = orientation,
      tracks = track_layout(signal = track('right', width = 3, limits = c(0, 50)))) +
      geom_chr_track(data = d, track = 'signal', geom = ggiraph::geom_col_interactive,
        ggplot2::aes(chr = Chr, position = Pos, value = Value, data_id = ID), width = 8)
    expect_true(inherits(utils::tail(p$layers, 1)[[1]]$geom, 'GeomCol'))
    expect_no_error(ggplot2::ggplotGrob(p))
    w <- as_ideogram_widget(p, width = 120, height = 60)
    expect_match(w$x$html, 'data-id=.a.')
    expect_match(w$x$html, 'data-id=.b.')
  }
})
