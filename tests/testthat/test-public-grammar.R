test_that("a chromosome collection creates the same plot through either constructor", {
  k <- data.frame(Chr = c("A", "B", "C"), Start = 0, End = c(100, 90, 80))
  a <- ggplot2::ggplot(k) + geom_chr(chr = c("C", "A"),
    orientation = "horizontal", names = TRUE, axis = TRUE)
  b <- ggideogram(k, chr = c("C", "A"), orientation = "horizontal", axis = TRUE)
  expect_s3_class(a, "ggplot")
  expect_identical(a$coordinates$layout$chrom$.chr, c("C", "A"))
  expect_identical(a$coordinates$layout$data$karyotype$.chr, k$Chr)
  expect_equal(ggplot2::ggplot_build(a)$data, ggplot2::ggplot_build(b)$data)
  expect_no_warning(ggplot2::ggplotGrob(a))
  base <- ggideogram(k, orientation = "horizontal", show_names = FALSE)
  selected <- base + geom_chr(chr = c("C", "A"), component = c("name", "axis"))
  expect_equal(length(selected$layers), length(base$layers) + 4)
  name <- Filter(function(layer) inherits(layer$stat, "StatChromosomeName"), selected$layers)[[1]]
  expect_identical(name$data$name_chr, c("C", "A"))
  expect_no_warning(ggplot2::ggplotGrob(selected))
})

test_that("selected bp axes retain their chromosomes after loci and tracks are added", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  d <- data.frame(Chr = c("A", "B"), Pos = c(20, 40), Value = c(2, 4))
  for (orientation in c("horizontal", "vertical", "circular")) {
    base <- ggideogram(k, orientation = orientation, reverse_chr = "B", show_names = FALSE)
    selected <- base + geom_chr(chr = "B", component = "axis",
      side = "left", breaks = c(0, 40, 80), units = "bp")
    replayed <- selected + geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos)) +
      geom_track(track = "signal", data = d, side = "right",
        mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_point())
    for (plot in list(selected, replayed)) {
      ticks <- Filter(function(layer) inherits(layer$geom, "GeomIdeogramTick"), plot$layers)
      expect_length(ticks, 1)
      expect_equal(nrow(ticks[[1]]$data), 3)
      source <- unproject_chr_point(plot, transform(ticks[[1]]$data, Chr = "B"), "Chr", "x", "y")
      expect_equal(source$.position, c(0, 40, 80))
      expect_identical(plot$coordinates$layout$chrom$.chr, k$Chr)
      expect_no_warning(ggplot2::ggplotGrob(plot))
    }
  }
})

test_that("locus adapters retain native prototypes and circular interval endpoints", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 80), Start = c(11, 61),
    End = c(30, 90), Name = c("gene1", "gene2"))
  native <- ggplot2::geom_point(size = 1.7, stroke = .26,
    position = ggplot2::position_nudge(x = 4))
  for (orientation in c("horizontal", "vertical", "circular")) {
    p <- ggideogram(k, orientation = orientation, show_names = FALSE) +
      geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos), geom = native) +
      geom_locus(data = d, ggplot2::aes(chr = Chr, start = Start, end = End),
        geom = ggplot2::geom_segment(linewidth = .23)) +
      geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos, label = Name),
        geom = "text")
    point <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))
    segment <- which(vapply(p$layers, function(layer) inherits(layer$stat, "StatChrInterval"), logical(1)))
    text <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1)))
    expect_s3_class(p$layers[[point]]$geom, "GeomPoint")
    expect_equal(p$layers[[point]]$aes_params$size, 1.7)
    expect_equal(p$layers[[point]]$position$x, 4)
    expect_s3_class(p$layers[[segment]]$geom, "GeomSegment")
    expect_equal(p$layers[[text]]$aes_params$size, 3)
    b <- ggplot2::ggplot_build(p)
    expect_equal(b$data[[segment]]$start, d$Start)
    expect_equal(b$data[[segment]]$end, d$End)
    if (orientation == "circular") expect_true(all(b$data[[segment]]$ideogram_projected))
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(native$data, ggplot2::waiver())
})

test_that("default track space and font roles are coherent", {
  expect_equal(geom_track()$width, 3)
  expect_equal(geom_track(side = "overlay")$width, 1)
  expect_equal(geom_track(side = "inner")$side, "left")
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  p <- ggideogram(k, orientation = "horizontal", axis = TRUE)
  expect_equal(p$coordinates$layout$chromosome_gap, 4.3)
  sizes <- vapply(Filter(function(layer) inherits(layer$geom, "GeomText"), p$layers),
    function(layer) layer$aes_params$size, numeric(1))
  expect_equal(sort(unique(sizes)), c(2.8, 3.2))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("track edge markers remain attached to circular node ports after track edits", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 75), ID = c("first", "last"))
  nodes <- chr_nodes(d, ggplot2::aes(id = ID, chr = Chr, position = Pos))
  p <- ggideogram(k, orientation = "circular", radius = 20, show_names = FALSE,
    tracks = list(signal = geom_track(side = "inner", width = 2),
      ports = geom_track(side = "inner", width = 1, gap = .3))) +
    geom_chrlink(nodes = nodes, data = data.frame(From = "first", To = "last"),
      ggplot2::aes(from = From, to = To)) +
    geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos),
      track = "ports", track_position = 1, size = 1.7)
  for (plot in list(p, p + geom_track(track = "ports", width = 3, gap = .8))) {
    build <- ggplot2::ggplot_build(plot)
    i <- which(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))
    j <- which(vapply(plot$layers, function(layer) inherits(layer$stat, "StatCartesianPair"), logical(1)))
    xy <- plot$coordinates$transform(build$data[[i]], build$layout$panel_params[[1]])
    links <- plot$coordinates$transform(build$data[[j]], build$layout$panel_params[[1]])
    expect_equal(xy$x, c(links$x[1], utils::tail(links$xend, 1)))
    expect_equal(xy$y, c(links$y[1], utils::tail(links$yend, 1)))
    expect_equal(build$data[[i]]$position, d$Pos)
    expect_equal(build$data[[i]]$size, rep(1.7, 2))
    expect_no_warning(ggplot2::ggplotGrob(plot))
  }
})

test_that("chromosome fills honor selection without changing source data or prior layers", {
  k <- data.frame(Chr = c("A", "B", "C"), Start = 0, End = c(100, 80, 60),
    Source = c("assembly-A", "assembly-B", "assembly-C"))
  d <- data.frame(Chr = c("A", "B", "B"), Start = c(11, 21, 51),
    End = c(30, 40, 70), ID = c("a", "b1", "b2"), Class = c("x", "y", "z"))
  original <- d
  for (orientation in c("horizontal", "vertical", "circular")) {
    base <- ggideogram(k, orientation = orientation, reverse_chr = "B",
      show_names = FALSE, tracks = list(body = geom_track(side = "overlay", width = .8)))
    p <- base + geom_chr(data = d, chr = "B", component = "fill", track = "body",
      ggplot2::aes(chr = Chr, start = Start, end = End, fill = Class))
    fills <- Filter(function(layer) inherits(layer$geom, "GeomIdeogramRect"), p$layers)
    expect_length(fills, 1)
    fill <- fills[[1]]$data
    selected <- d[d$Chr == "B", ]
    expect_equal(fill[names(d)], selected)
    expect_equal(fill$.source_start, selected$Start)
    expect_equal(fill$.source_end, selected$End)
    anchor <- project_chr_point(p, transform(selected, Pos = (Start - 1 + End) / 2), "Chr", "Pos")
    expect_equal((fill$.xmin + fill$.xmax) / 2, anchor$.x)
    expect_equal((fill$.ymin + fill$.ymax) / 2, anchor$.y)
    expect_identical(p$coordinates$layout$data, base$coordinates$layout$data)
    expect_identical(p$coordinates$layout$chrom, base$coordinates$layout$chrom)
    expect_equal(ggplot2::ggplot_build(p)$data[seq_along(base$layers)],
      ggplot2::ggplot_build(base)$data)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(d, original)
  expect_error(base + geom_chr(data = d, chr = "missing", component = "fill", track = "body",
    ggplot2::aes(chr = Chr, start = Start, end = End)), "Unknown selected chromosome")
  invalid <- d; invalid$End[1] <- 101
  expect_error(base + geom_chr(data = invalid, chr = "B", component = "fill", track = "body",
    ggplot2::aes(chr = Chr, start = Start, end = End)), "source chromosome bounds")
})

test_that("body and fill components share selection and replay in the final layout", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  d <- data.frame(Chr = c("A", "B"), Start = c(11, 21), End = c(30, 40),
    ID = c("a", "b"), Class = c("x", "y"))
  mapping <- ggplot2::aes(chr = Chr, start = Start, end = End, fill = Class)
  for (orientation in c("horizontal", "vertical", "circular")) {
    base <- ggideogram(k, orientation = orientation, reverse_chr = "B", show_names = FALSE,
      tracks = list(body = geom_track(side = "overlay", width = .8)))
    combined <- base + geom_chr(data = d, mapping = mapping, chr = "B", track = "body",
      component = c("body", "fill"), linewidth = .2)
    separate <- base + geom_chr(chr = "B", component = "body", linewidth = .2) +
      geom_chr(data = d, mapping = mapping, chr = "B", component = "fill", track = "body", linewidth = .2)
    expect_equal(ggplot2::ggplot_build(combined)$data, ggplot2::ggplot_build(separate)$data)
    expect_identical(combined$coordinates$layout$data, base$coordinates$layout$data)
    added_body <- combined$layers[[length(base$layers) + 1L]]
    expect_identical(unique(added_body$data$.chr), "B")
    points <- transform(d, Pos = (Start - 1 + End) / 2, Value = c(2, 4))
    replayed <- combined + geom_track(track = "signal", data = points, side = "right",
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_point())
    fill <- Filter(function(layer) inherits(layer$geom, "GeomIdeogramRect"), replayed$layers)[[1]]$data
    expect_identical(fill$ID, "b")
    anchor <- project_chr_point(replayed, points[2, ], "Chr", "Pos")
    expect_equal((fill$.xmin + fill$.xmax) / 2, anchor$.x)
    expect_equal((fill$.ymin + fill$.ymax) / 2, anchor$.y)
    expect_identical(replayed$coordinates$layout$data, base$coordinates$layout$data)
    expect_no_warning(ggplot2::ggplotGrob(replayed))
    scoped <- ggideogram(k, orientation = orientation, show_names = FALSE,
      tracks = list(body = geom_track(side = "overlay", width = .8, data = d,
        mapping = mapping, layers = list(geom_chr(chr = "B", component = "fill")))))
    expect_identical(Filter(function(layer) inherits(layer$geom, "GeomIdeogramRect"),
      scoped$layers)[[1]]$data$ID, "b")
    expect_no_warning(ggplot2::ggplotGrob(scoped))
  }
  plain <- ggplot2::ggplot(k) + geom_chr(data = d, mapping = mapping, chr = "B",
    component = c("body", "fill"), track = "body", names = FALSE,
    orientation = "horizontal", tracks = list(body = geom_track(side = "overlay", width = .8)))
  expect_identical(plain$coordinates$layout$chrom$.chr, "B")
  expect_identical(plain$coordinates$layout$data$karyotype$Chr, k$Chr)
  expect_identical(Filter(function(layer) inherits(layer$geom, "GeomIdeogramRect"),
    plain$layers)[[1]]$data$ID, "b")
  expect_no_warning(ggplot2::ggplotGrob(plain))
})

test_that("additive geom_chr axis matches constructor gap, family and plain ggplot dispatch", {
  k <- data.frame(Chr = c("A", "B", "C", "D"), Start = 0, End = c(100e3, 90e3, 80e3, 70e3))
  ctor <- ggideogram(k, ncol = 4, axis = TRUE, base_family = "serif")
  added <- ggideogram(k, ncol = 4, base_family = "serif") + geom_chr(component = "axis")
  expect_equal(added$coordinates$layout$chromosome_gap, ctor$coordinates$layout$chromosome_gap)
  expect_equal(added$coordinates$layout$chromosome_gap, 8)
  axis_text <- Filter(function(l) inherits(l$geom, "GeomIdeogramAxisText"), added$layers)[[1]]
  expect_identical(axis_text$aes_params$family, "serif")
  expect_true(isTRUE(axis_text$geom_params$check_overlap))

  axis_only <- ggplot2::ggplot(k, ggplot2::aes(chr = Chr, start = Start, end = End)) +
    geom_chr(component = "axis", units = "kb", breaks = c(0, 50e3, 100e3))
  expect_false(any(vapply(axis_only$layers, function(l) inherits(l$geom, "GeomChromosome"), logical(1))))
  only_text <- Filter(function(l) inherits(l$geom, "GeomIdeogramAxisText"), axis_only$layers)[[1]]
  expect_true(all(grepl("kb$", only_text$data$label)))
  expect_no_warning(ggplot2::ggplotGrob(axis_only))

  body_and_axis <- ggplot2::ggplot(k) +
    geom_chr(component = c("body", "axis"), ncol = 4, units = "kb")
  expect_equal(body_and_axis$coordinates$layout$chromosome_gap, 8)
  expect_error(ggideogram(k) + geom_chr(chr = c("A", "A")), "unique")
  expect_error(ggideogram(k) + geom_chr(chr = character(0)), "unique|non-empty")
})

test_that("disjoint implicit spread labels share one chr_labels track and beside annotations reserve marker_extent", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 100))
  d_a <- data.frame(Chr = "A", Pos = c(20, 40), Label = c("gA1", "gA2"))
  d_b <- data.frame(Chr = "B", Pos = c(30, 60), Label = c("gB1", "gB2"))
  p <- ggideogram(k) +
    geom_locus(data = d_a, ggplot2::aes(chr = Chr, position = Pos, label = Label),
      geom = "text", position = "spread") +
    geom_locus(data = d_b, ggplot2::aes(chr = Chr, position = Pos, label = Label),
      geom = "text", position = "spread")
  expect_equal(p$coordinates$layout$tracks$id, "chr_labels")
  expect_false("chr_labels_" %in% p$coordinates$layout$tracks$id)
  expect_equal(unname(p$coordinates$layout$track_extent["right"]), 16.2)

  # Same-chromosome implicit spread calls still stack safely on separate lanes
  p_same <- ggideogram(k) +
    geom_locus(data = d_a, ggplot2::aes(chr = Chr, position = Pos, label = Label),
      geom = "text", position = "spread") +
    geom_locus(data = transform(d_a, Pos = Pos + 5),
      ggplot2::aes(chr = Chr, position = Pos, label = Label),
      geom = "text", position = "spread")
  expect_true("chr_labels_" %in% p_same$coordinates$layout$tracks$id)

  # Untracked beside interval and identity text reserve marker_extent
  p_iv <- ggideogram(k, show_names = FALSE) +
    geom_chr_interval(data = data.frame(Chr = "A", Start = 10, End = 40),
      ggplot2::aes(chr = Chr, start = Start, end = End),
      placement = "beside", side = "right", gap = 0.6) +
    geom_track(track = "sig", side = "right", width = 2, gap = 0.2)
  expect_equal(unname(p_iv$coordinates$layout$marker_extent["right"]), 1.2)
  expect_equal(track_table_row(p_iv$coordinates$layout, "sig")$low_offset, 0.6 + 1.2 + 0.2)

  p_txt <- ggideogram(k, show_names = FALSE) +
    geom_locus(data = d_a, ggplot2::aes(chr = Chr, position = Pos, label = Label),
      geom = "text", side = "left", gap = 0.7)
  expect_equal(unname(p_txt$coordinates$layout$marker_extent["left"]), 1.4)
})

test_that("geom_locus forwards third-party interactive aesthetics, na.rm and explicit inherit.aes without warnings", {
  skip_if_not_installed("ggiraph")
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(25, 75), Tip = c("t1", "t2"), ID = c("id1", "id2"))
  expect_no_warning({
    p <- ggideogram(k, show_names = FALSE) +
      geom_locus(
        data = d,
        ggplot2::aes(chr = Chr, position = Pos, tooltip = Tip, data_id = ID),
        geom = ggiraph::geom_point_interactive,
        na.rm = TRUE
      )
  })
  layer <- utils::tail(p$layers, 1)[[1]]
  expect_true(isTRUE(layer$geom_params$na.rm))
  expect_true(isTRUE(layer$stat_params$na.rm))
  expect_false(isTRUE(layer$inherit.aes))
  built <- ggplot2::ggplot_build(p)
  expect_identical(utils::tail(built$data, 1)[[1]]$tooltip, d$Tip)
  expect_identical(utils::tail(built$data, 1)[[1]]$data_id, d$ID)
})
