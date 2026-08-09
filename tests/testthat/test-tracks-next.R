track_kar <- data.frame(
  Chr = c("A", "B"), Start = 0, End = c(100, 80),
  CE_start = c(40, 30), CE_end = c(55, 45)
)

track_signal <- data.frame(
  Chr = rep(c("A", "B"), each = 5),
  Pos = rep(seq(10, 70, length.out = 5), 2),
  Value = c(1:5, 5:1),
  Sample = rep(c("s1", "s2"), each = 5),
  stringsAsFactors = FALSE
)

test_that("declared tracks own stable named layout space", {
  specs <- track_layout(
    first = track(side = "right", width = 1, gap = 0.2),
    second = track(side = "right", width = 0.5, gap = 0.3),
    left = track(side = "left", width = 0.8, gap = 0.1),
    overlay = track(side = "overlay", width = 0.6)
  )
  layout <- ideogram_layout(
    as_ideogram_data(track_kar), chromosome_width = 1,
    chromosome_gap = 0.7, tracks = specs)

  expect_s3_class(specs, "ideogram_track_layout")
  expect_equal(layout$tracks$id, c("first", "second", "left", "overlay"))
  expect_equal(layout$tracks$low_offset, c(0.7, 2, -0.6, -0.3))
  expect_equal(layout$tracks$high_offset, c(1.7, 2.5, -1.4, 0.3))
  expect_equal(layout$track_extent, c(left = 0.9, right = 2))

  centres <- layout$chrom$.axis_start_x
  first_outer <- centres[1] + 0.5 + layout$track_extent[["right"]]
  second_outer <- centres[2] - 0.5 - layout$track_extent[["left"]]
  expect_equal(second_outer - first_outer, 0.7)
  expect_equal(diff(layout$bounds$x), 2 * (0.9 + 1 + 2) + 0.7)
})

test_that("track declarations reject ambiguous or impossible geometry", {
  expect_error(track_layout(track()), "non-empty name")
  expect_error(track(side = "overlay", gap = 0.1), "gap = 0")
  expect_error(track(width = 0), "positive")
  expect_error(track(limits = c(1, 1)), "increasing")
  expect_error(track(transform = "not_a_transform"), "Invalid")
  expect_error(
    ideogram_layout(
      as_ideogram_data(track_kar),
      tracks = track_layout(wide = track("overlay", width = 2))
    ),
    "wider than"
  )
})

test_that("dynamic layers keep their ordinary geom classes", {
  specs <- track_layout(
    signal = track("right", limits = c(0, 5)),
    bars = track("left", limits = c(0, 5)),
    labels = track("overlay", width = 0.8, limits = c(0, 5.5))
  )
  third_party_point <- function(...) ggplot2::geom_point(...)
  p <- ggideogram(track_kar, tracks = specs, show_names = FALSE) +
    geom_chr_track(
      data = track_signal, geom = third_party_point, track = "signal",
      mapping = ggplot2::aes(
        chr = Chr, position = Pos, value = Value, colour = Sample),
      size = 1
    ) +
    geom_track_line(
      data = track_signal, track = "signal",
      ggplot2::aes(
        chr = Chr, position = Pos, value = Value, colour = Sample),
      linewidth = 0.4
    ) +
    geom_track_col(
      data = track_signal, track = "bars",
      ggplot2::aes(chr = Chr, position = Pos, value = Value, fill = Sample),
      width = 5, position = "identity"
    ) +
    geom_track_tile(
      data = track_signal, track = "labels",
      ggplot2::aes(chr = Chr, position = Pos, value = Value, fill = Sample),
      width = 5, height = 0.2
    )

  classes <- utils::tail(
    vapply(p$layers, function(layer) class(layer$geom)[1], character(1)), 4)
  expect_equal(
    unname(classes),
    c("GeomPoint", "GeomLine", "GeomCol", "GeomIdeogramTile"))
  expect_s3_class(p$layers[[length(p$layers)]]$geom, "GeomTile")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("overlay heatmap tiles are clipped to the chromosome silhouette", {
  heat <- data.frame(
    Chr = rep(c("A", "B"), each = 2),
    Pos = c(2.5, 97.5, 2.5, 77.5),
    Value = 0.5,
    Width = 5,
    Fill = rep(c("low", "high"), 2)
  )
  specs <- track_layout(
    inside = track("overlay", width = 1, limits = c(0, 1)),
    outside = track("right", width = 1, limits = c(0, 1))
  )
  base <- ggideogram(
    track_kar, tracks = specs, curve_points = 18, show_names = FALSE)
  inside <- base + geom_track_tile(
    data = heat, track = "inside",
    ggplot2::aes(
      chr = Chr, position = Pos, value = Value,
      width = Width, fill = Fill),
    height = 1
  )
  outside <- base + geom_track_tile(
    data = heat, track = "outside",
    ggplot2::aes(
      chr = Chr, position = Pos, value = Value,
      width = Width, fill = Fill),
    height = 1
  )
  unrestricted <- base + geom_track_tile(
    data = heat, track = "inside", clip = "off",
    ggplot2::aes(
      chr = Chr, position = Pos, value = Value,
      width = Width, fill = Fill),
    height = 1
  )

  expect_s3_class(inside$layers[[length(inside$layers)]]$geom,
                  "GeomIdeogramTile")
  expect_s3_class(inside$layers[[length(inside$layers)]]$geom, "GeomTile")
  expect_identical(
    class(outside$layers[[length(outside$layers)]]$geom)[1], "GeomTile")
  expect_identical(
    class(unrestricted$layers[[length(unrestricted$layers)]]$geom)[1],
    "GeomTile")
  expect_equal(inside$coordinates$layout$curve_points, 18L)
  expect_no_error(ggplot2::ggplotGrob(inside))
  expect_error(
    base + geom_track_tile(
      data = heat, track = "outside", clip = "on",
      ggplot2::aes(
        chr = Chr, position = Pos, value = Value,
        width = Width, fill = Fill),
      height = 1),
    "only for an overlay track"
  )
})

test_that("native markers can occupy a declared track", {
  specs <- track_layout(markers = track("right", width = 0.6, gap = 0.2))
  marker <- data.frame(Chr = "A", Pos = 50, Type = "locus")
  p <- ggideogram(track_kar[1, ], tracks = specs, show_names = FALSE) +
    geom_chr_link(
      data = marker,
      ggplot2::aes(chr = Chr, position = Pos, colour = Type),
      track = "markers"
    ) +
    geom_chr_marker(
      data = marker,
      ggplot2::aes(chr = Chr, position = Pos, colour = Type),
      track = "markers", size = 1
    )
  built <- ggplot2::ggplot_build(p)
  link <- built$data[[length(built$data) - 1L]]
  point <- built$data[[length(built$data)]]
  link_xy <- p$coordinates$transform(
    link, built$layout$panel_params[[1]])
  point_xy <- p$coordinates$transform(
    point, built$layout$panel_params[[1]])

  expect_s3_class(p$layers[[length(p$layers)]]$geom, "GeomPoint")
  expect_equal(link_xy$xend, point_xy$x)
  expect_equal(link_xy$yend, point_xy$y)
  expect_gt(point_xy$x, link_xy$x)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("automatic ranges update all layers without mutating the base plot", {
  specs <- track_layout(signal = track("right", value_scale = "per_track"))
  base <- ggideogram(track_kar, tracks = specs, show_names = FALSE)
  narrow <- track_signal[track_signal$Value >= 2 & track_signal$Value <= 4, ]
  wide <- data.frame(Chr = c("A", "B"), Pos = c(50, 40),
                     Value = c(0, 10))
  first <- base + geom_track_line(
    data = narrow, track = "signal",
    ggplot2::aes(chr = Chr, position = Pos, value = Value)
  )
  combined <- first + geom_track_point(
    data = wide, track = "signal",
    ggplot2::aes(chr = Chr, position = Pos, value = Value)
  )

  expect_null(base$coordinates$layout$track_ranges[["track:signal"]]$limits)
  expect_equal(first$coordinates$layout$track_ranges[["track:signal"]]$limits,
               c(2, 4))
  expect_equal(combined$coordinates$layout$track_ranges[["track:signal"]]$limits,
               c(0, 10))
  expect_no_error(ggplot2::ggplotGrob(combined))
})

test_that("per chromosome and global value scopes are explicit", {
  per_chr <- track_layout(signal = track("right", value_scale = "per_chr"))
  values <- data.frame(
    Chr = rep(c("A", "B"), each = 2), Pos = c(10, 80, 10, 70),
    Value = c(0, 10, 100, 200)
  )
  p_chr <- ggideogram(track_kar, tracks = per_chr, show_names = FALSE) +
    geom_track_point(
      data = values, track = "signal",
      ggplot2::aes(chr = Chr, position = Pos, value = Value)
    )
  expect_equal(
    p_chr$coordinates$layout$track_ranges[["track:signal:chr:A"]]$limits,
    c(0, 10))
  expect_equal(
    p_chr$coordinates$layout$track_ranges[["track:signal:chr:B"]]$limits,
    c(100, 200))

  global <- track_layout(
    one = track("right", value_scale = "global"),
    two = track("left", value_scale = "global")
  )
  p_global <- ggideogram(track_kar, tracks = global, show_names = FALSE) +
    geom_track_point(
      data = values[1:2, ], track = "one",
      ggplot2::aes(chr = Chr, position = Pos, value = Value)
    ) +
    geom_track_point(
      data = values[3:4, ], track = "two",
      ggplot2::aes(chr = Chr, position = Pos, value = Value)
    )
  expect_equal(p_global$coordinates$layout$track_ranges$global$limits,
               c(0, 200))
})

test_that("track projection supports transforms, reversal and orientation", {
  specs <- track_layout(
    signal = track(
      "right", width = 1, limits = c(1, 100),
      transform = "log10", reverse = TRUE)
  )
  semantic <- as_ideogram_data(track_kar[1, ])
  vertical <- ideogram_layout(semantic, tracks = specs)
  horizontal <- ideogram_layout(semantic, tracks = specs,
                                orientation = "horizontal")
  data <- data.frame(Chr = "A", Pos = 50, Value = c(1, 10, 100))
  v <- project_chr_track(vertical, data, "Chr", "Pos", "Value", "signal")
  h <- project_chr_track(horizontal, data, "Chr", "Pos", "Value", "signal")

  expect_equal(v$.value_fraction, c(1, 0.5, 0))
  expect_true(all(diff(v$.x) < 0))
  expect_equal(length(unique(v$.y)), 1L)
  expect_equal(length(unique(h$.x)), 1L)
  expect_true(all(diff(h$.y) > 0))
  expect_equal(h$.x, diff(vertical$bounds$y) - v$.y)
  expect_equal(h$.y, diff(vertical$bounds$x) - v$.x)
  expect_true(all(h$.y >= horizontal$bounds$y[1]))
  expect_true(all(h$.y <= horizontal$bounds$y[2]))
})

test_that("stacked columns are isolated by chromosome and register stack range", {
  stacked <- data.frame(
    Chr = rep(c("A", "B"), each = 4),
    Pos = rep(c(30, 30, 60, 60), 2),
    Value = rep(c(2, 3, 1, 4), 2),
    Class = rep(c("x", "y"), 4)
  )
  specs <- track_layout(bars = track("right"))
  p <- ggideogram(track_kar, tracks = specs, show_names = FALSE) +
    geom_track_col(
      data = stacked, track = "bars",
      ggplot2::aes(
        chr = Chr, position = Pos, value = Value, fill = Class),
      width = 6
    )
  built <- ggplot2::ggplot_build(p)
  bars <- utils::tail(built$data, 1)[[1]]

  expect_equal(p$coordinates$layout$track_ranges[["track:bars"]]$limits,
               c(0, 5))
  expect_lte(max(bars$ymax), 5)
  expect_true(length(unique(bars$ideogram_position_shift)) == 2L)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("area, ribbon, boxplot and violin use standard ggplot2 geoms", {
  line <- transform(track_signal, Low = pmax(0, Value - 0.5),
                    High = Value + 0.5)
  set.seed(42)
  distributions <- data.frame(
    Chr = rep(c("A", "B"), each = 40),
    Pos = rep(rep(c(25, 60), each = 20), 2),
    Value = c(stats::rnorm(40, 2), stats::rnorm(40, 3)),
    Class = rep(rep(c("x", "y"), each = 20), 2)
  )
  specs <- track_layout(
    area = track("right", limits = c(0, 6)),
    ribbon = track("left", limits = c(0, 6)),
    box = track("right", limits = c(-1, 6)),
    violin = track("left", limits = c(-1, 6))
  )
  p <- ggideogram(track_kar, tracks = specs, show_names = FALSE) +
    geom_track_area(
      data = line, track = "area",
      ggplot2::aes(chr = Chr, position = Pos, value = Value),
      fill = "steelblue", alpha = 0.4
    ) +
    geom_track_ribbon(
      data = line, track = "ribbon",
      ggplot2::aes(chr = Chr, position = Pos, ymin = Low, ymax = High),
      fill = "orange", alpha = 0.4
    ) +
    geom_track_boxplot(
      data = distributions, track = "box",
      ggplot2::aes(
        chr = Chr, position = Pos, value = Value, fill = Class),
      width = 2
    ) +
    geom_track_violin(
      data = distributions, track = "violin",
      ggplot2::aes(
        chr = Chr, position = Pos, value = Value, fill = Class),
      width = 2
    )
  classes <- utils::tail(
    vapply(p$layers, function(layer) class(layer$geom)[1], character(1)), 4)

  expect_equal(
    unname(classes), c("GeomRibbon", "GeomRibbon", "GeomBoxplot", "GeomViolin"))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplotGrob(p))

  automatic <- ggideogram(
    track_kar,
    tracks = track_layout(area = track("right")),
    show_names = FALSE)
  expect_error(
    automatic + geom_track_area(
      data = line, track = "area",
      ggplot2::aes(chr = Chr, position = Pos, value = Value)),
    "explicit"
  )
})

test_that("track axes use resolved values and physical text and ticks", {
  specs <- track_layout(signal = track("right", limits = c(0, 5)))
  p <- ggideogram(track_kar, tracks = specs, show_names = FALSE) +
    geom_track_line(
      data = track_signal, track = "signal",
      ggplot2::aes(chr = Chr, position = Pos, value = Value)
    ) +
    geom_track_axis(
      "signal", chr = "A", breaks = c(0, 2.5, 5),
      tick_length = 2, size = 2.1)
  classes <- utils::tail(
    vapply(p$layers, function(layer) class(layer$geom)[1], character(1)), 3)

  expect_equal(unname(classes),
               c("GeomSegment", "GeomIdeogramTick",
                 "GeomIdeogramAxisText"))
  tick_layer <- p$layers[[length(p$layers) - 1L]]
  expect_equal(tick_layer$geom_params$tick_length, 2)
  expect_true(all(tick_layer$data$ny < 0))
  text_layer <- p$layers[[length(p$layers)]]
  expect_equal(text_layer$aes_params$size, 2.1)
  expect_equal(text_layer$geom_params$tick_length, 2)
  expect_equal(text_layer$geom_params$label_gap, 0.25)
  expect_equal(text_layer$data[c("nx", "ny")],
               tick_layer$data[c("nx", "ny")])
  expect_no_error(ggplot2::ggplotGrob(p))

  auto <- ggideogram(
    track_kar, tracks = track_layout(signal = track("right")),
    show_names = FALSE)
  expect_error(auto + geom_track_axis("signal", chr = "A"), "no resolved")
})

test_that("track scales and guides remain native ggplot2 components", {
  specs <- track_layout(signal = track("right", limits = c(0, 5)))
  p <- ggideogram(track_kar, tracks = specs, show_names = FALSE) +
    geom_track_point(
      data = track_signal, track = "signal",
      ggplot2::aes(
        chr = Chr, position = Pos, value = Value,
        colour = Sample, size = Value)
    ) +
    ggplot2::scale_size_continuous(range = c(0.8, 2.2)) +
    ggplot2::labs(colour = "Sample", size = "Signal") +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 9))
  built <- ggplot2::ggplot_build(p)
  point <- utils::tail(built$data, 1)[[1]]

  expect_equal(range(point$size), c(0.8, 2.2))
  expect_true(length(built$plot$guides$aesthetics) >= 1L)
  expect_no_error(ggplot2::ggplotGrob(p))
})
