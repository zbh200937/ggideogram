multi_karyotype <- data.frame(
  Genome = rep(c('A', 'B'), each = 3), Assembly = 'v1',
  Chr = rep(c('1', '2', '3'), 2), Start = 0,
  End = c(100, 80, 60, 95, 75, 65), Homolog = rep(c('1', '2', '3'), 2),
  Label = paste0(rep(c('1', '2', '3'), 2), rep(c('A', 'B'), each = 3)))
multi_semantic <- function(k = multi_karyotype) as_ideogram_data(k,
  ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome,
               assembly = Assembly, homolog = Homolog, label = Label))

test_that('composite identities preserve genome, assembly and source fields', {
  keys <- chr_key(c('ab', 'a', 'ab'), c('c', 'bc', 'c'), c('x', 'x', 'y'))
  expect_equal(length(unique(keys)), 3)
  expect_equal(chr_key('水稻', c('01', '02')), chr_key(rep('水稻', 2), c('01', '02')))
  expect_error(chr_key(letters[1:3], letters[1:2]), 'length')
  d <- multi_semantic()
  expect_equal(d$karyotype$.chr, with(multi_karyotype, chr_key(Genome, Chr, Assembly)))
  expect_equal(d$karyotype$.chr_name, multi_karyotype$Chr)
  expect_equal(d$karyotype$.display_name, multi_karyotype$Label)
  expect_equal(d$karyotype$End, multi_karyotype$End)
  k <- multi_karyotype[c(1, 1), ]; k$Assembly <- c('v1', 'v2')
  expect_equal(length(unique(multi_semantic(k)$karyotype$.chr)), 2)
  expect_error(multi_semantic(multi_karyotype[c(1, 1), ]), 'unique')
})

test_that('grouped order preserves explicit homolog slots and bp comparisons', {
  d <- multi_semantic(multi_karyotype[-5, ])
  g <- ideogram_layout(d, order_by = 'genome', genome_order = c('B', 'A'))
  expect_equal(g$chrom$.genome, c('B', 'B', 'A', 'A', 'A'))
  expect_equal(g$chrom$.col, c(0, 2, 0, 1, 2))
  expect_equal(length(unique(g$chrom$.units_per_bp)), 1)
  expect_true(g$length_comparable)
  h <- ideogram_layout(d, order_by = 'homolog', homolog_order = c('3', '2', '1'))
  expect_equal(h$chrom$.homolog_group, c('3', '3', '2', '1', '1'))
  expect_equal(h$chrom$.col, c(0, 1, 0, 0, 1))
  expect_equal(ideogram_layout(d, order_by = 'homolog', ncol = 1)$ncol, 1)
  expect_false(ideogram_layout(d, scale_length = 'per_chr')$length_comparable)
  plain <- as_ideogram_data(multi_karyotype,
    ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome))
  expect_error(ideogram_layout(plain, order_by = 'homolog'), 'homolog')
  expect_error(ideogram_layout(d, order_by = 'genome', genome_order = 'A'), 'every observed')
})

test_that('reverse display retains names, track side and original bp ticks', {
  d <- multi_semantic()
  keys <- d$karyotype$.chr
  for (orientation in c('vertical', 'horizontal')) {
    base <- ideogram_layout(d, orientation = orientation,
      tracks = track_layout(signal = track('right', width = 2, limits = c(0, 1))))
    rev <- ideogram_layout(d, orientation = orientation, reverse_chr = keys[4],
      tracks = base$track_layout)
    expect_equal(rev$chrom$.start, base$chrom$.start)
    expect_equal(rev$chrom$.name_x, base$chrom$.name_x)
    expect_equal(rev$chrom$.name_y, base$chrom$.name_y)
    points <- data.frame(Key = keys[c(1, 4)], Pos = c(20, 20), Value = 1)
    projected <- project_chr_track(rev, points, 'Key', 'Pos', 'Value', 'signal')
    expect_equal(unproject_chr_point(rev, projected, 'Key', '.x', '.y')$.position, points$Pos)
    normal <- chromosome_right_normal(rev, rev$chrom[c(1, 4), ])
    expect_equal(normal$nx[1], normal$nx[2])
    expect_equal(normal$ny[1], normal$ny[2])
    p <- ggideogram(d, orientation = orientation, reverse_chr = keys[4],
      tracks = rev$track_layout, axis = keys[4], axis_breaks = c(0, 20, 80)) +
      geom_track_point(data = points, track = 'signal',
        ggplot2::aes(chr = Key, position = Pos, value = Value)) +
      geom_chr_marker(data = points, ggplot2::aes(chr = Key, position = Pos))
    expect_true(inherits(utils::tail(p$layers, 1)[[1]]$geom, 'GeomPoint'))
    expect_no_error(ggplot2::ggplotGrob(p))
    labels <- Filter(function(layer) inherits(layer$geom, 'GeomIdeogramAxisText'), p$layers)
    expect_equal(as.character(labels[[1]]$data$label), c('0 bp', '20 bp', '80 bp'))
  }
})

test_that('same-name annotation and paired endpoints match the declared genome', {
  d <- multi_semantic()
  key <- d$karyotype$.chr
  pts <- data.frame(Genome = c('A', 'B'), Chr = '1', Assembly = 'v1', Pos = 25)
  p <- ggideogram(d, order_by = 'genome') + geom_chr_marker(data = pts,
    ggplot2::aes(chr = chr_key(Genome, Chr, Assembly), position = Pos))
  b <- ggplot2::ggplot_build(p)
  expect_equal(utils::tail(b$data, 1)[[1]]$ideogram_chr, key[c(1, 4)])
  pairs <- data.frame(Key1 = key[1], Key2 = key[4], Pos1 = 20, Pos2 = 30)
  p <- p + geom_chr_connection(data = pairs,
    ggplot2::aes(chr1 = Key1, position1 = Pos1, chr2 = Key2, position2 = Pos2))
  layer <- utils::tail(p$layers, 1)[[1]]
  expect_equal(unproject_chr_point(p$coordinates$layout, layer$data,
    '.pair_chr2', '.pair_xend', '.pair_yend')$.position, 30)
  expect_no_error(ggplot2::ggplotGrob(p))
  stats <- data.frame(Key = key, Value = multi_karyotype$End)
  q <- ggplot2::ggplot(stats, ggplot2::aes(Key, Value)) +
    ggplot2::geom_col() + scale_x_chromosome(d)
  built <- ggplot2::ggplot_build(q)
  expect_equal(as.character(built$layout$panel_scales_x[[1]]$get_labels()), multi_karyotype$Label)
  expect_true(inherits(q$layers[[1]]$geom, 'GeomCol'))
  expect_no_error(ggplot2::ggplotGrob(q))
})
