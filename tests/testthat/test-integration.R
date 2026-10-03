test_that('chromosome attachments retain source identity after ordering and local views', {
  k <- data.frame(Genome = c('B', 'A'), Chr = '1', Start = 0, End = c(90, 100))
  semantic <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome))
  keys <- semantic$karyotype$.chr
  metadata <- data.frame(Key = keys, Species = c('species-B', 'species-A'))
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(semantic, orientation = orientation,
      order_by = 'genome', genome_order = c('A', 'B'), reverse_chr = keys[1])
    p <- base + chr_data(metadata, by = c(chr = 'Key'), name = 'species')
    layout <- p$coordinates$layout
    expect_equal(layout$data$karyotype$Species, metadata$Species)
    expect_equal(layout$chrom$Species, rev(metadata$Species))
    expect_equal(get_ideogram_attachment(p, 'species')$.chr, keys)
    expect_null(base$coordinates$layout$data$karyotype$Species)
    expect_no_error(ggplot2::ggplotGrob(p))

    view <- chr_view(semantic, keys[1], 10, 80)
    p <- ggideogram(view, orientation = orientation, reverse_chr = keys[1]) +
      chr_data(metadata[1, ], by = c(chr = 'Key'), name = 'species')
    expect_equal(p$coordinates$layout$data$karyotype$Species, c('species-B', NA_character_))
    expect_equal(p$coordinates$layout$chrom$Species, 'species-B')
    expect_equal(p$coordinates$layout$data$karyotype$End, k$End)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that('distribution positions stay independent across chromosomes and retain within-chromosome dodge', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = 100)
  d <- data.frame(Chr = rep(c('A', 'B'), each = 6), Pos = 50,
    Value = c(1, 2, 3, 2, 3, 4, 6, 7, 8, 7, 8, 9))
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(k, orientation = orientation, reverse_chr = 'B', max_chr_length = 40,
      tracks = track_layout(dist = track('right', width = 2, limits = c(0, 10))))
    long <- if (orientation == 'vertical') 'y' else 'x'
    for (wrapper in list(geom_track_boxplot, geom_track_violin)) {
      p <- base + wrapper(ggplot2::aes(chr = Chr, position = Pos, value = Value), d,
        track = 'dist', width = 1.8)
      actual <- utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]
      expected <- project_chr_track(base$coordinates$layout, d, 'Chr', 'Pos', 'Value', 'dist')
      expected_long <- expected[[paste0('.', long)]][!duplicated(expected$Chr)]
      expect_equal(as.numeric(tapply(actual[[long]], actual$group, unique)),
        expected_long)
      expect_no_warning(ggplot2::ggplotGrob(p))
      grouped <- rbind(transform(d, Category = 'one'), transform(d, Value = Value + 0.2, Category = 'two'))
      p <- base + wrapper(ggplot2::aes(chr = Chr, position = Pos, value = Value, fill = Category),
        grouped, track = 'dist', width = 1.8)
      adjusted <- utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]
      positions <- as.numeric(tapply(adjusted[[long]], adjusted$group, unique))
      expect_length(positions, 4)
      expect_equal(sort(c(mean(positions[1:2]), mean(positions[3:4]))),
        sort(expected_long))
      expect_true(all(abs(diff(positions[c(1, 2)])) > 0))
    }
  }
})

test_that('generic native segment and rectangle tracks preserve all chromosome coordinates', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = 100)
  d <- data.frame(Chr = c('A', 'B'), Pos = 20, Value = 2,
    End = 80, EndValue = 5, Left = 20, Right = 80, Low = 1, High = 5)
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(k, orientation = orientation, reverse_chr = 'B',
      tracks = track_layout(signal = track('right', width = 2)))
    for (constant in c(FALSE, TRUE)) {
      segment <- if (constant) geom_chr_track(
        ggplot2::aes(chr = Chr, position = Pos, value = Value,
          xend = End / 2, yend = EndValue / 2), d,
        geom = ggplot2::geom_segment, track = 'signal', xend = 80, yend = 5)
      else geom_chr_track(ggplot2::aes(chr = Chr, position = Pos, value = Value,
        xend = End, yend = EndValue), d, geom = ggplot2::geom_segment, track = 'signal')
      p <- base + segment
      build <- ggplot2::ggplot_build(p)
      data <- utils::tail(build$data, 1)[[1]]
      expect_equal(p$coordinates$layout$track_ranges[['track:signal']]$limits, c(2, 5))
      expect_equal(data$xend - data$ideogram_position_shift, d$End)
      expected <- project_chr_track(p$coordinates$layout,
        transform(d, Pos = End, Value = EndValue), 'Chr', 'Pos', 'Value', 'signal')
      expected <- p$coordinates$transform(data.frame(x = expected$.x, y = expected$.y),
        build$layout$panel_params[[1]])
      actual <- p$coordinates$transform(data, build$layout$panel_params[[1]])
      expect_equal(actual$xend, expected$x); expect_equal(actual$yend, expected$y)
      expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, 'GeomSegment')
      expect_no_error(ggplot2::ggplotGrob(p))
    }
    p <- base + geom_chr_track(ggplot2::aes(chr = Chr, position = Pos, value = Value,
      xmin = Left, xmax = Right, ymin = Low, ymax = High), d,
      geom = ggplot2::geom_rect, track = 'signal')
    build <- ggplot2::ggplot_build(p)
    data <- utils::tail(build$data, 1)[[1]]
    expect_equal(data$xmin - data$ideogram_position_shift, d$Left)
    expect_equal(data$xmax - data$ideogram_position_shift, d$Right)
    expect_equal(p$coordinates$layout$track_ranges[['track:signal']]$limits, c(1, 5))
    expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, 'GeomRect')
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})
