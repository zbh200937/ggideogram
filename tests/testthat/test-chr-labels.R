draw_chr_labels_test <- function(plot, mm = 120) {
  file <- tempfile(fileext = '.pdf')
  grDevices::cairo_pdf(file, width = mm / 25.4, height = mm / 25.4)
  on.exit({ grDevices::dev.off(); unlink(file) })
  built <- ggplot2::ggplot_build(plot)
  grid::grid.draw(ggplot2::ggplot_gtable(built))
  states <- Filter(function(state) !is.null(state$options$padding),
    as.list(built$layout$panel_params[[1]]$ideogram_labels))
  list(built = built, states = states,
    positions = do.call(rbind, lapply(states, function(state) state$cache$positions)))
}

chr_labels_test_anchors <- function(state) {
  layout <- state$coord$layout
  d <- state$data
  p <- project_positions_checked(layout, as.character(d$chr), d$position)
  side <- track_table_row(layout, state$options$track)$side
  p <- offset_chr_points(layout, as.character(d$chr), p, side, layout$chromosome_width / 2)
  p <- state$coord$transform(as.data.frame(p), state$panel)
  mm <- utils::tail(state$device, 2)
  cbind(x = p$x * mm[1], y = p$y * mm[2])
}

expect_chr_labels_native_baseline <- function(state) {
  z <- state$cache$positions
  mm <- utils::tail(state$device, 2)
  text <- state$cache$text
  expect_equal(as.numeric(text$x), z$baseline_x / mm[1])
  expect_equal(as.numeric(text$y), z$baseline_y / mm[2])
  expect_equal(as.numeric(text$hjust), z$hjust)
  expect_equal(as.numeric(text$vjust), z$vjust)
  expect_equal(text$gp$fontsize, z$size * ggplot2::.pt)
  angle <- as.numeric(text$rot) * pi / 180
  distance <- (0.5 - as.numeric(text$hjust)) * z$width
  expect_equal(z$x, as.numeric(text$x) * mm[1] + distance * cos(angle))
  expect_equal(z$y, as.numeric(text$y) * mm[2] + distance * sin(angle))
  last <- !duplicated(state$cache$leaders$group, fromLast = TRUE)
  expect_equal(state$cache$leaders$xend[last], z$baseline_x / mm[1])
  expect_equal(state$cache$leaders$yend[last], z$baseline_y / mm[2])
}

chr_labels_test_overlap <- function(a, b) {
  direction <- function(z) cbind(cos(z$angle * pi / 180), sin(z$angle * pi / 180))
  ua <- direction(a); ub <- direction(b)
  va <- cbind(-ua[, 2], ua[, 1]); vb <- cbind(-ub[, 2], ub[, 1])
  axes <- rbind(ua, va, ub, vb)
  delta <- c(b$x - a$x, b$y - a$y)
  all(vapply(seq_len(nrow(axes)), function(i) {
    axis <- axes[i, ]
    half_a <- abs(sum(ua * axis)) * a$width / 2 + abs(sum(va * axis)) * a$height / 2
    half_b <- abs(sum(ub * axis)) * b$width / 2 + abs(sum(vb * axis)) * b$height / 2
    abs(sum(delta * axis)) < half_a + half_b - 1e-8
  }, logical(1)))
}

test_that('ordered labels are native text and segments and retain real positions', {
  skip_if_not(capabilities('cairo'))
  d <- data.frame(Chr = 'A', Pos = c(48, 49, 50, 51, 52, 53),
    Label = c('G1', 'Gene2', 'NAC001', 'LongGene4', 'G5', 'ARV1'),
    OriginalStart = 1001:1006)
  for (orientation in c('vertical', 'horizontal')) for (side in c('left', 'right'))
    for (reverse in c(FALSE, TRUE)) {
    base <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100),
      orientation = orientation, show_names = FALSE, reverse_chr = if (reverse) 'A' else character(),
      tracks = track_layout(labels = track(side, width = 40, gap = 1)))
    p <- base + geom_chr_labels(ggplot2::aes(chr = Chr, position = Pos, label = Label),
      d, track = 'labels', size = 3)
    expect_s3_class(utils::tail(p$layers, 2)[[1]]$geom, 'GeomSegment')
    expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, 'GeomText')
    expect_equal(utils::tail(p$layers, 1)[[1]]$data$OriginalStart, d$OriginalStart)
    result <- draw_chr_labels_test(p)
    expect_equal(result$positions$position, d$Pos)
    expect_equal(nrow(result$positions), nrow(d))
    anchors <- chr_labels_test_anchors(result$states[[1]])
    expect_equal(result$positions$anchor_x, anchors[, 'x'])
    expect_equal(result$positions$anchor_y, anchors[, 'y'])
    expect_chr_labels_native_baseline(result$states[[1]])
    axis <- if (orientation == 'vertical') 'baseline_y' else 'baseline_x'
    normal <- if (orientation == 'vertical') 'baseline_x' else 'baseline_y'
    near <- if (orientation == 'vertical') 'near_x' else 'near_y'
    expect_equal(result$positions[[normal]], rep(result$positions[[normal]][1], nrow(d)))
    expect_equal(abs(result$positions[[normal]] - result$positions[[near]]), rep(5, nrow(d)))
    expect_equal(result$positions$angle, rep(if (orientation == 'vertical') 0 else 90, nrow(d)))
    extent <- result$positions$height
    i <- order(result$positions[[axis]])
    expect_true(all(diff(result$positions[[axis]][i]) >=
      (utils::head(extent[i], -1) + utils::tail(extent[i], -1)) / 2 + 0.6 - 1e-8))
    expect_true(result$states[[1]]$cache$fits)
  }
})

test_that('circle label near edges share a ring while varied words extend radially', {
  skip_if_not(capabilities('cairo'))
  k <- data.frame(Genome = c('G1', 'G2'), Chr = '1', Start = 0, End = 100)
  keys <- with(k, chr_key(Genome, Chr))
  source <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome))
  d <- data.frame(Key = rep(keys, each = 4), Pos = rep(c(49, 49.2, 50, 51), 2),
    Label = rep(c('G1', 'ARV1', 'NAC001', 'GeneABC3'), 2), Size = rep(c(2.6, 3), 4))
  for (side in c('inner', 'outer')) for (clockwise in c(FALSE, TRUE)) {
    base <- ggideogram(source, orientation = 'circular', clockwise = clockwise,
      radius = 30, start_angle = 31, reverse_chr = keys[2], axis = TRUE, axis_n = 2,
      tracks = track_layout(labels = track(side, width = if (side == 'inner') 20 else 34)))
    p <- base + geom_chr_labels(ggplot2::aes(chr = Key, position = Pos, label = Label, size = Size),
      d, track = 'labels') + ggplot2::scale_size_identity()
    result <- lapply(c(120, 200), function(mm) draw_chr_labels_test(p, mm))
    for (drawn in result) {
    positions <- drawn$positions
    mm <- utils::tail(drawn$states[[1]]$device, 2)
    radius <- sqrt((positions$baseline_x - mm[1] / 2)^2 + (positions$baseline_y - mm[2] / 2)^2)
    near_radius <- sqrt((positions$near_x - mm[1] / 2)^2 + (positions$near_y - mm[2] / 2)^2)
    centre_radius <- sqrt((positions$x - mm[1] / 2)^2 + (positions$y - mm[2] / 2)^2)
    expect_equal(radius, rep(radius[1], nrow(d)), tolerance = 1e-8)
    expect_equal(abs(radius - near_radius), rep(5, nrow(d)), tolerance = 1e-8)
    expect_equal(abs(centre_radius - radius), positions$width / 2, tolerance = 1e-8)
    expect_equal(positions$position, d$Pos)
    expect_equal(positions$chr, d$Key)
    expect_equal(positions$size, d$Size)
    expect_true(all(positions$angle >= -90 & positions$angle <= 90))
    anchors <- chr_labels_test_anchors(drawn$states[[1]])
    expect_equal(positions$anchor_x, anchors[, 'x'])
    expect_equal(positions$anchor_y, anchors[, 'y'])
    state <- drawn$states[[1]]
    expect_chr_labels_native_baseline(state)
    first <- seq(1, nrow(state$cache$leaders), by = state$coord$layout$curve_points - 1)
    expect_equal(as.numeric(state$cache$segments$x0)[first], anchors[, 'x'] / mm[1])
    expect_equal(as.numeric(state$cache$segments$y0)[first], anchors[, 'y'] / mm[2])
    pairs <- utils::combn(seq_len(nrow(positions)), 2)
    expect_false(any(apply(pairs, 2, function(i)
      chr_labels_test_overlap(positions[i[1], ], positions[i[2], ]))))
    expect_true(drawn$states[[1]]$cache$fits)
  }
    expect_equal(result[[1]]$positions$width, result[[2]]$positions$width)
    expect_equal(result[[1]]$positions$height, result[[2]]$positions$height)
  }
})

test_that('own label occupancy preserves tracks and local source records', {
  skip_if_not(capabilities('cairo'))
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  view <- chr_view(as_ideogram_data(k), 'A', 40, 60)
  d <- data.frame(Chr = 'A', Pos = c(10, 49, 50, 51, 90), Label = LETTERS[1:5], ID = 1:5)
  base <- ggideogram(view, orientation = 'circular', clockwise = FALSE, reverse_chr = 'A',
    tracks = track_layout(signal = track('inner', width = 2, limits = c(0, 1))))
  p <- base + geom_chr_labels(ggplot2::aes(chr = Chr, x = Pos, label = Label), d)
  expect_equal(p$coordinates$layout$tracks[1, ], base$coordinates$layout$tracks[1, ])
  expect_equal(p$coordinates$layout$data$karyotype, base$coordinates$layout$data$karyotype)
  expect_true('chr_labels' %in% p$coordinates$layout$tracks$id)
  expect_equal(utils::tail(p$layers, 1)[[1]]$data$ID, 2:4)
  result <- draw_chr_labels_test(p)
  expect_equal(result$positions$position, c(49, 50, 51))
  expect_error(base + geom_chr_labels(ggplot2::aes(chr = Chr, position = Pos, label = Label),
    transform(d, Pos = c(10, 49, 50, 51, 110))), 'source')
  expect_equal(d$Pos, c(10, 49, 50, 51, 90))
})

test_that('insufficient physical capacity warns and keeps all records', {
  skip_if_not(capabilities('cairo'))
  d <- data.frame(Chr = 'A', Pos = 50, Label = paste('Long gene name', 1:20))
  p <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100), orientation = 'horizontal',
    show_names = FALSE, tracks = track_layout(labels = track('right', width = 2))) +
    geom_chr_labels(ggplot2::aes(chr = Chr, position = Pos, label = Label), d, track = 'labels')
  expect_warning(result <- draw_chr_labels_test(p, 70), 'available display span')
  expect_equal(nrow(result$positions), nrow(d))
  expect_equal(result$positions$position, d$Pos)
  expect_equal(p$coordinates$layout$tracks$width, 2)
  expect_false(result$states[[1]]$cache$fits)
})

test_that('automatic outer and inner bands hold native short and long gene IDs', {
  skip_if_not(capabilities('cairo'))
  d <- data.frame(Chr = rep(c('A', 'B'), each = 3), Pos = rep(c(49, 50, 51), 2),
    Label = c('NAC001', 'Os01g0100100', 'ORUFI01G00010', 'ARV1', 'NAC001', 'Os01g0100100'))
  for (side in c('inner', 'outer')) {
    p <- ggideogram(data.frame(Chr = c('A', 'B'), Start = 0, End = 100),
      orientation = 'circular', reverse_chr = 'B', tracks = track_layout(
        value = track('outer', width = 2), genes = track('inner', width = 2))) +
      geom_locus(geom = "text", position = "spread", ggplot2::aes(chr = Chr, position = Pos, label = Label), d, side = side)
    expect_equal(track_table_row(p$coordinates$layout, 'chr_labels')$width,
      if (side == 'inner') 14 else 32)
    for (mm in c(120, 200)) {
      expect_no_warning(result <- draw_chr_labels_test(p, mm))
      expect_true(result$states[[1]]$cache$fits)
      expect_equal(result$positions$position, d$Pos)
      expect_chr_labels_native_baseline(result$states[[1]])
    }
  }
})

test_that('connection height is physical and explicit track widths remain unchanged', {
  skip_if_not(capabilities('cairo'))
  d <- data.frame(Chr = 'A', Pos = c(49, 51), Label = c('G1', 'ARV1'))
  p <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100), orientation = 'vertical',
    tracks = track_layout(labels = track('left', width = 24))) +
    geom_locus(geom = "text", position = "spread", ggplot2::aes(chr = Chr, position = Pos, label = Label), d,
      track = 'labels', connection_height = 7)
  expect_equal(p$coordinates$layout$tracks$width, 24)
  result <- draw_chr_labels_test(p)
  expect_equal(abs(result$positions$baseline_x - result$positions$near_x), rep(7, nrow(d)))
  expect_error(geom_chr_labels(connection_height = -1), 'connection_height')
})

test_that('own label tracks relayout chromosome collections and survive replay', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = 100)
  d <- data.frame(Chr = c('A', 'B'), Pos = c(49, 51), Label = c('alpha', 'beta'))
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    p <- ggideogram(k, orientation = orientation,
      tracks = list(signal = geom_track(side = 'left', width = 2))) +
      geom_locus(ggplot2::aes(chr = Chr, position = Pos), d, size = 1.7) +
      geom_locus(geom = "text", position = "spread", ggplot2::aes(chr = Chr, position = Pos, label = Label), d)
    actual <- p$coordinates$layout
    expected <- ideogram_layout(actual$data, orientation = orientation,
      tracks = actual$track_layout)
    expect_equal(actual$chrom$.axis_start_x, expected$chrom$.axis_start_x)
    expect_equal(actual$chrom$.axis_start_y, expected$chrom$.axis_start_y)
    q <- p + geom_track(track = 'other', side = 'right', width = 2)
    expect_equal(sum(q$coordinates$layout$tracks$id == 'chr_labels'), 1)
    expect_false('chr_labels_' %in% q$coordinates$layout$tracks$id)
    expect_equal(sum(vapply(q$layers, function(layer) inherits(layer$geom, 'GeomChrLabelsText'), logical(1))), 1)
    marker <- Filter(function(layer) inherits(layer$geom, 'GeomPoint'), ggplot2::ggplot_build(q)$plot$layers)
    expect_length(marker, 1)
  }
})

test_that('text spreading accepts native constructors and owned layer data and styles', {
  skip_if_not(capabilities('cairo'))
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  d <- data.frame(Chr = 'A', Pos = c(49, 51), Label = c('NAC001', 'ARV1'), ID = c('g1', 'g2'))
  native <- ggplot2::geom_text(data = d,
    mapping = ggplot2::aes(x = Pos, label = Label),
    size = 3.1, colour = '#4477AA', family = 'sans', fontface = 'bold', check_overlap = TRUE)
  for (orientation in c('horizontal', 'vertical', 'circular')) {
    p <- ggideogram(k, orientation = orientation, show_names = FALSE) +
      geom_locus(ggplot2::aes(chr = Chr), geom = native, position = 'spread',
        label_width = 40, connection_height = 6)
    expect_no_warning(result <- draw_chr_labels_test(p))
    state <- result$states[[1]]
    expect_equal(state$data$position, d$Pos)
    expect_equal(state$options$source$ID, d$ID)
    expect_equal(state$cache$text$gp$fontsize, rep(3.1 * ggplot2::.pt, 2))
    expect_equal(state$cache$text$gp$col, rep('#4477AA', 2))
    expect_true(state$cache$text$check.overlap)
    expect_chr_labels_native_baseline(state)
    q <- ggideogram(k, orientation = orientation, show_names = FALSE) +
      geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos, label = Label),
        geom = ggplot2::geom_text, position = 'spread', label_width = 40,
        size = 8, size.unit = 'pt', colour = '#CC6677')
    expect_no_warning(points <- draw_chr_labels_test(q))
    expect_equal(points$states[[1]]$cache$text$gp$fontsize, rep(8, 2))
    expect_equal(points$positions$size, rep(8 / ggplot2::.pt, 2))
  }
  expect_identical(native$data, d)
  expect_equal(native$aes_params$size, 3.1)
  expect_error(geom_locus(geom = 'point', position = 'spread'), 'text geom')
})
