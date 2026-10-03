repel_test_plot <- function(data, orientation = 'vertical', columns = 1,
                            overflow = 'warn', max_labels = Inf) {
  ggideogram(data.frame(Chr = 'A', Start = 0, End = 100),
    orientation = orientation, show_names = FALSE,
    tracks = track_layout(labels = track('right', width = 12, gap = 1))) +
    geom_chr_text_repel(data = data, track = 'labels', columns = columns,
      ggplot2::aes(chr = Chr, position = Pos, label = Label, priority = Priority),
      size = 2.9, seed = 12, max.iter = 5000, force_pull = 0.2,
      overflow = overflow, max_labels = max_labels)
}

draw_repel_test <- function(plot, width = 120, height = 120) {
  file <- tempfile(fileext = '.pdf')
  grDevices::cairo_pdf(file, width = width / 25.4, height = height / 25.4)
  on.exit({ grDevices::dev.off(); unlink(file) })
  built <- ggplot2::ggplot_build(plot)
  grid::grid.draw(ggplot2::ggplot_gtable(built))
  states <- as.list(built$layout$panel_params[[1]]$ideogram_labels)
  list(states = states, panel = built$layout$panel_params[[1]],
    positions = do.call(rbind, lapply(states, function(s) s$cache$positions)))
}

test_that('device labels keep true anchors and native text/segment geoms', {
  skip_if_not_installed('ggrepel')
  d <- data.frame(Chr = 'A', Pos = c(48, 49, 50, 51, 52, 53),
    Label = paste('Gene', 1:6), Priority = 1)
  for (orientation in c('vertical', 'horizontal')) {
    p <- repel_test_plot(d, orientation)
    leader <- Filter(function(layer) inherits(layer$geom, 'GeomChrRepelSegment'), p$layers)[[1]]
    text <- Filter(function(layer) inherits(layer$geom, 'GeomChrRepelText'), p$layers)[[1]]
    expect_s3_class(leader$geom, 'GeomSegment')
    expect_s3_class(text$geom, 'GeomTextRepel')
    drawn <- draw_repel_test(p)
    labels <- drawn$positions
    expect_equal(nrow(labels), nrow(d))
    expect_equal(sort(labels$position), d$Pos)
    expected <- p$coordinates$transform(data.frame(
      x = leader$data$.anchor_x,
      y = leader$data$.anchor_y), drawn$panel)
    expect_equal(labels$anchor_x, expected$x)
    expect_equal(labels$anchor_y, expected$y)
    expect_equal(ggideogram:::chr_repel_conflicts(labels,
      list(xlim = c(0, 1), ylim = c(0, 1))), integer())
    expect_true(any(abs((labels$top + labels$bottom) / 2 - labels$anchor_y) > 0.001) ||
      any(abs((labels$left + labels$right) / 2 - labels$anchor_x) > 0.001))
  }
})

test_that('repeated rendering is deterministic and font size survives two devices', {
  skip_if_not_installed('ggrepel')
  d <- data.frame(Chr = 'A', Pos = c(30, 31, 32, 70),
    Label = c('one', 'two', 'three', 'four'), Priority = 1)
  p <- repel_test_plot(d)
  a <- draw_repel_test(p)
  b <- draw_repel_test(p)
  c <- draw_repel_test(p, width = 85, height = 90)
  expect_equal(a$positions, b$positions)
  expect_equal(a$positions$fontsize, c$positions$fontsize)
  expect_equal(a$positions$position, c$positions$position)
  expect_false(identical(a$positions$top, c$positions$top))
})

test_that('left/right tracks, columns and explicit limits retain priority', {
  skip_if_not_installed('ggrepel')
  d <- data.frame(Chr = 'A', Pos = seq(10, 80, 10),
    Label = letters[1:8], Priority = c(1, 1, 1, 1, 1, 1, 3, 3))
  p <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100), show_names = FALSE,
    tracks = track_layout(l = track('left', width = 12, gap = 1),
      r = track('right', width = 12, gap = 1))) +
    geom_chr_text_repel(data = d, track = c(left = 'l', right = 'r'), columns = 2,
      ggplot2::aes(chr = Chr, position = Pos, label = Label))
  expect_equal(sum(vapply(p$layers, function(l) inherits(l$geom, 'GeomChrRepelText'), logical(1))), 4)
  expect_no_warning(draw_repel_test(p))
  expect_warning(limited <- repel_test_plot(d, max_labels = 2), 'Label limit omitted')
  leader <- Filter(function(layer) inherits(layer$geom, 'GeomChrRepelSegment'), limited$layers)[[1]]
  expect_equal(sort(leader$data$Label), c('g', 'h'))
  d$Side <- rep(c('left', 'right'), each = 4)
  assigned <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100), show_names = FALSE,
    tracks = track_layout(l = track('left', width = 12), r = track('right', width = 12))) +
    geom_chr_text_repel(data = d, track = c(left = 'l', right = 'r'),
      ggplot2::aes(chr = Chr, position = Pos, label = Label, side = Side))
  leaders <- Filter(function(layer) inherits(layer$geom, 'GeomChrRepelSegment'), assigned$layers)
  expect_equal(leaders[[2]]$data$Pos, seq(50, 80, 10))
})

test_that('overflow is explicit and omission keeps the highest priority', {
  skip_if_not_installed('ggrepel')
  d <- data.frame(Chr = 'A', Pos = rep(50, 12),
    Label = paste('gene', seq_len(12)), Priority = c(10, rep(0, 11)))
  p <- repel_test_plot(d)
  expect_warning(draw_repel_test(p, width = 70, height = 40), 'do not fit')
  p <- repel_test_plot(d, overflow = 'omit')
  expect_warning(result <- draw_repel_test(p, width = 70, height = 40), 'omitted')
  expect_lt(nrow(result$positions), nrow(d))
  expect_true('L1' %in% result$positions$id)
  expect_equal(ggideogram:::chr_repel_conflicts(result$positions,
    list(xlim = c(0, 1), ylim = c(0, 1))), integer())
})

test_that('an annular column detects a text box crossing its inner boundary between corners', {
  skip_if_not_installed('ggrepel')
  d <- data.frame(Chr = 'A', Pos = 0.025, Label = 'AT1G01020 (ARV1)')
  base <- ggideogram(data.frame(Chr = 'A', Start = 0, End = 100),
    orientation = 'circular', radius = 10, start_angle = 0, gap_angle = 2,
    show_names = FALSE,
    tracks = track_layout(labels = track('inner', width = 2, gap = 1))) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  label <- function(overflow) geom_chr_text_repel(
    ggplot2::aes(chr = Chr, position = Pos, label = Label), d,
    track = 'labels', size = 2.9, seed = 9, force = 0, force_pull = 0,
    max.iter = 1, background = NA, position = ggplot2::position_nudge(y = -1.3),
    overflow = overflow)
  expect_warning(warned <- draw_repel_test(base + label('warn'), 80, 80), 'do not fit')
  state <- warned$states[[1]]
  box <- warned$positions
  polygon <- state$coord$transform(state$options$scope_polygon, state$panel)
  expect_true(all(ggideogram:::points_inside_polygon(
    c(box$left, box$right, box$right, box$left),
    c(box$bottom, box$bottom, box$top, box$top), polygon)))
  expect_false(ggideogram:::points_inside_polygon(
    (box$left + box$right) / 2, (box$bottom + box$top) / 2, polygon))
  expect_equal(box$position, d$Pos)
  expect_warning(omitted <- draw_repel_test(base + label('omit'), 80, 80), 'omitted')
  expect_equal(nrow(omitted$positions), 0L)
  expect_equal(omitted$states[[1]]$cache$omitted, 'L1')
})
test_that('locus text repulsion reserves annotation space and retains native styles', {
  skip_if_not_installed('ggrepel')
  k <- data.frame(Chr = 'A', Start = 0, End = 100)
  d <- data.frame(Chr = 'A', Pos = c(25, 75), Label = c('NAC001', 'ARV1'))
  p <- ggideogram(k, orientation = 'vertical', show_names = FALSE) +
    geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos, label = Label),
      geom = ggplot2::geom_text(size = 3.1, colour = '#4477AA'),
      position = 'repel', label_width = 30, seed = 2)
  expect_equal(track_table_row(p$coordinates$layout, 'chr_labels')$width, 30)
  expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, 'GeomTextRepel')
  expect_s3_class(utils::tail(p$layers, 2)[[1]]$geom, 'GeomSegment')
  b <- ggplot2::ggplot_build(p)
  text <- utils::tail(b$data, 1)[[1]]
  expect_equal(text$label, d$Label)
  expect_equal(text$size, rep(3.1, 2))
  expect_equal(text$colour, rep('#4477AA', 2))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that('circular locus text stays in its sector with nodes, an inset and late track changes', {
  skip_if_not_installed('ggrepel')
  k <- read.delim(system.file('extdata', 'rice-comparison-karyotype.tsv',
    package = 'ggideogram'), check.names = FALSE)
  source <- read.delim(system.file('extdata', 'rice-comparison-genes.tsv',
    package = 'ggideogram'), check.names = FALSE)
  d <- source[c(head(which(source$Genome == 'rice'), 3),
    head(which(source$Genome == 'wild'), 3)), ]
  d$Key <- chr_key(d$Genome, d$Chr, d$Assembly)
  d$Mid <- (d$Start - 1 + d$End) / 2
  d$Value <- 2:7
  model <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End,
    genome = Genome, assembly = Assembly, homolog = Homolog, label = Label))
  nodes <- chr_nodes(d, ggplot2::aes(id = Gene, chr = Key, start = Start, end = End))
  edges <- data.frame(From = d$Gene[1:3], To = d$Gene[4:6])
  child <- ggplot2::ggplot(data.frame(x = c(10, 30), y = c(1000, 2000),
    Group = c('child-low', 'child-high')), ggplot2::aes(x, y, colour = Group)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_colour_manual(values = c('child-low' = 'gold', 'child-high' = 'purple')) +
    ggplot2::coord_cartesian(xlim = c(0, 50), ylim = c(0, 3000)) +
    ggplot2::theme_minimal()
  child_data <- ggplot2::ggplot_build(child)$data
  inset <- geom_locus_inset(child, data = d[1, ],
    ggplot2::aes(chr = Key, start = Start, end = End), track = 'child_space',
    width = grid::unit(24, 'mm'), height = grid::unit(22, 'mm'))
  make <- function(mode, width = 2, gap = .2, limits = NULL) {
    text <- do.call(geom_locus, c(list(data = d,
      mapping = ggplot2::aes(chr = Key, position = Mid, label = Gene, colour = Genome),
      geom = 'text', position = mode, label_width = 48, size = 3),
      if (mode == 'repel') list(seed = 11) else list()))
    ggideogram(model, orientation = 'circular', radius = 35, show_names = FALSE,
      reverse_chr = model$karyotype$.chr[2],
      tracks = list(signal = geom_track(side = 'right', width = width, gap = gap,
        limits = limits, data = d,
        ggplot2::aes(chr = Key, x = Mid, y = Value, colour = Genome),
        layers = list(ggplot2::geom_line(linewidth = .35), ggplot2::geom_point(size = 1.8))),
        child_space = geom_track(side = 'right', width = 20))) +
      geom_chrlink(data = edges, nodes = nodes, ggplot2::aes(from = From, to = To),
        side1 = 'right', side2 = 'right', gap = 0, curvature = 0, linewidth = .25) +
      geom_locus(data = d, ggplot2::aes(chr = Key, position = Mid),
        side = 'right', gap = 0, size = 1.4) +
      geom_locus(data = d[1, ], ggplot2::aes(chr = Key, position = Mid),
        track = 'child_space', track_position = 0, size = 1.1) +
      text +
      ggplot2::scale_colour_manual(values = c(rice = 'blue', wild = 'red'))
  }
  select <- function(p, pred) which(vapply(p$layers, pred, logical(1)))
  xy <- function(p, b, i) p$coordinates$transform(b$data[[i]], b$layout$panel_params[[1]])
  for (mode in c('spread', 'repel')) {
    base <- make(mode)
    with_child <- base + inset
    expect_equal(ggplot2::ggplot_build(with_child)$data[seq_along(base$layers)],
      ggplot2::ggplot_build(base)$data)
    expect_equal(with_child$coordinates$layout$track_ranges, base$coordinates$layout$track_ranges)
    after <- with_child + geom_track(track = 'signal', width = 3.5, gap = .8, limits = c(0, 10))
    reference <- make(mode, width = 3.5, gap = .8, limits = c(0, 10)) + inset
    b <- ggplot2::ggplot_build(after)
    ref <- ggplot2::ggplot_build(reference)
    for (i in seq_along(b$data)) {
      fields <- setdiff(intersect(names(b$data[[i]]), names(ref$data[[i]])),
        c('ideogram_inset_layer', 'plot'))
      expect_equal(b$data[[i]][fields], ref$data[[i]][fields])
    }
    text <- select(after, function(l) inherits(l$geom, 'GeomChrLabelsText') ||
      inherits(l$geom, 'GeomChrRepelText'))
    rows <- do.call(rbind, lapply(text, function(i) after$layers[[i]]$data))
    at <- match(d$Gene, rows$Gene)
    expect_equal(rows$.label_chr[at], d$Key)
    expect_equal(rows$.label_position[at], nodes$.node_position)
    expect_equal(rows[c('Gene', 'Start', 'End')][at, ], d[c('Gene', 'Start', 'End')],
      ignore_attr = TRUE)
    for (i in text) {
      expect_s3_class(after$layers[[i]]$geom, if (mode == 'spread') 'GeomText' else 'GeomTextRepel')
      expect_equal(b$data[[i]]$size, rep(3, nrow(b$data[[i]])))
    }
    marker <- select(after, function(l) inherits(l$stat, 'StatChrMarker') &&
      is.null(l$stat_params$track) && nrow(l$data) == 6)
    link <- select(after, function(l) inherits(l$geom, 'GeomSegment') &&
      is.data.frame(l$data) && '.pair_id' %in% names(l$data))
    expect_s3_class(after$layers[[marker]]$geom, 'GeomPoint')
    expect_s3_class(after$layers[[link]]$geom, 'GeomSegment')
    m <- xy(after, b, marker); e <- xy(after, b, link)
    expect_equal(e$x, m$x[1:3]); expect_equal(e$y, m$y[1:3])
    expect_equal(e$xend, m$x[4:6]); expect_equal(e$yend, m$y[4:6])
    child_anchor <- select(after, function(l) inherits(l$stat, 'StatChrMarker') &&
      identical(l$stat_params$track, 'child_space'))
    component <- select(after, function(l) inherits(l$geom, 'GeomChrInset'))
    expect_equal(xy(after, b, child_anchor)[c('x', 'y')],
      xy(after, b, component)[c('x', 'y')])
    expect_equal(b$data[[component]]$ideogram_position, nodes$.node_position[1])
    expect_equal(ggplot2::ggplot_build(child)$data, child_data)
    expect_equal(child$coordinates$limits, list(x = c(0, 50), y = c(0, 3000)))
    expect_no_warning(drawn <- draw_repel_test(after, width = 300, height = 260))
    expect_equal(nrow(drawn$positions), 6L)
    if (mode == 'repel') for (state in drawn$states) {
      box <- state$cache$positions
      polygon <- state$coord$transform(state$options$scope_polygon, state$panel)
      expect_true(all(vapply(seq_len(nrow(box)), function(i)
        ggideogram:::rectangle_inside_polygon(box[i, ], polygon), logical(1))))
      expect_equal(ggideogram:::chr_repel_conflicts(box,
        list(xlim = c(0, 1), ylim = c(0, 1))), integer())
      anchor <- state$coord$transform(data.frame(x = state$data$x_orig,
        y = state$data$y_orig), state$panel)
      expect_equal(box$anchor_x, anchor$x[box$index])
      expect_equal(box$anchor_y, anchor$y[box$index])
      expect_equal(box$fontsize, rep(3 * ggplot2::.pt, nrow(box)))
      expect_equal(box$position, state$options$source$.label_position[box$index])
    }
  }
})
