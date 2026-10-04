test_that('native path geoms retain missing observations as breaks', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = 100)
  observations <- data.frame(Chr = rep(c('A', 'B'), each = 3),
    Pos = rep(c(20, 50, 80), 2), Value = c(2, NA_real_, 8, 6, NA_real_, 3))
  source <- observations
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    for (geom in list(ggplot2::geom_line(), ggplot2::geom_path(), ggplot2::geom_step())) {
      p <- ggideogram(k, orientation = orientation, reverse_chr = 'B') +
        geom_track(track = 'signal', data = observations,
          mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = geom)
      built <- ggplot2::ggplot_build(p)
      actual <- utils::tail(built$data, 1)[[1]]
      expect_equal(actual$position, observations$Pos)
      expect_equal(actual$value, observations$Value)
      expect_equal(length(unique(actual$group)), 2L)
      layer <- utils::tail(p$layers, 1)[[1]]
      retained <- layer$geom$handle_na(actual, layer$geom_params)
      expect_equal(nrow(retained), 6L)
      projected <- built$layout$coord$transform(retained, built$layout$panel_params[[1]])
      expect_equal(which(is.na(projected$x) | is.na(projected$y)), c(2L, 5L))
      expect_no_warning(ggplot2::ggplotGrob(p))
    }
  }
  expect_identical(observations, source)
})

test_that('native rectangle geometry can cross local boundaries without changing source observations', {
  features <- read_chr_features(system.file('extdata', 'arabidopsis-first-genes.gff3',
    package = 'ggideogram'))
  genes <- features[features$Type == 'gene', ]
  k <- data.frame(Chr = '1', Start = 0, End = 30427671)
  semantic <- as_ideogram_data(k)
  counts <- bin_genome(genes, semantic, window = 2000, method = 'count')
  counts$Pos <- (counts$Start - 1 + counts$End) / 2
  view <- chr_view(semantic, '1', 4500, 15500)
  visible <- view_chr_data(counts, view, position = 'Pos')
  source <- visible
  expect_true(any(visible$Pos - visible$Width / 2 < view$view$start))
  expect_true(any(visible$Pos + visible$Width / 2 > view$view$end))
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    for (reverse in c(FALSE, TRUE)) {
      base <- ggideogram(view, orientation = orientation,
        reverse_chr = if (reverse) '1' else character(),
        tracks = track_layout(counts = track('right', width = 2, limits = c(0, 2))))
      p <- base + geom_track_col(
        ggplot2::aes(chr = Chr, position = Pos, value = Value, width = Width),
        visible, track = 'counts')
      layer <- utils::tail(p$layers, 1)[[1]]
      actual <- utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]
      expect_s3_class(layer$geom, 'GeomCol')
      expect_equal(actual$position, visible$Pos)
      expect_equal(actual$value, visible$Value)
      expect_equal(actual$xmax - actual$xmin, visible$Width)
      expect_no_warning(ggplot2::ggplotGrob(p))
      expect_error(base + geom_track_point(
        ggplot2::aes(chr = Chr, position = Pos, value = Value),
        transform(visible[1, ], Pos = 4499), track = 'counts'), 'outside')
    }
  }
  expect_identical(visible, source)
})

test_that('short final bins retain native widths in circular rectangle tracks', {
  data('human_karyotype', package = 'ggideogram')
  data('LTR_density', package = 'ggideogram')
  k <- human_karyotype[2, ]
  last <- utils::tail(LTR_density[LTR_density$Chr == k$Chr, ], 1)
  last$Pos <- (last$Start + last$End) / 2
  expect_lt(last$End - last$Start, 1000000)
  expect_gt(last$Pos + 500000, k$End)
  source <- last
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(k, orientation = orientation,
      tracks = track_layout(counts = track('right', limits = c(0, last$Value + 1))))
    p <- base + geom_track_col(
      ggplot2::aes(chr = Chr, position = Pos, value = Value),
      last, track = 'counts', width = 1000000)
    actual <- utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]
    expect_equal(actual$position, last$Pos)
    expect_equal(actual$xmax - actual$xmin, 1000000)
    expect_no_warning(ggplot2::ggplotGrob(p))
    expect_error(base + geom_chr_track(
      ggplot2::aes(chr = Chr, position = Pos, value = Value,
        xmin = Start, xmax = End + 1, ymin = 0, ymax = Value),
      last, geom = ggplot2::geom_rect, track = 'counts'), 'outside')
  }
  expect_identical(last, source)
})

test_that('identity geom constructors without a stat argument retain native columns', {
  native_col <- function(mapping = NULL, data = NULL, position = 'stack', ...,
                         na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
    ggplot2::geom_col(mapping = mapping, data = data, position = position,
      ..., na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes)
  }
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  observations <- data.frame(Chr = 'A', Pos = c(10, 30), Value = c(2, 4),
    Width = c(4, 6))
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    p <- ggideogram(k, orientation = orientation,
      tracks = track_layout(counts = track('right', limits = c(0, 5)))) +
      geom_chr_track(ggplot2::aes(chr = Chr, position = Pos, value = Value,
        width = Width), observations, geom = native_col, track = 'counts',
        position = 'stack', fill = 'steelblue')
    built <- ggplot2::ggplot_build(p)
    layer <- utils::tail(p$layers, 1)[[1]]
    actual <- utils::tail(built$data, 1)[[1]]
    projected <- built$layout$coord$transform(actual, built$layout$panel_params[[1]])
    expected <- project_chr_track(built$layout$coord$layout, observations,
      'Chr', 'Pos', 'Value', track = 'counts')
    expected <- built$layout$coord$transform(
      data.frame(x = expected$.x, y = expected$.y), built$layout$panel_params[[1]])
    expect_s3_class(layer$geom, 'GeomCol')
    expect_equal(actual$position, observations$Pos)
    expect_equal(actual$xmax - actual$xmin, observations$Width)
    expect_equal(projected$x, expected$x)
    expect_equal(projected$y, expected$y)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
})
