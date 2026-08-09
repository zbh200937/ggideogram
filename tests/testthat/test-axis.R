axis_kar <- data.frame(
  Chr = c("A", "B"), Start = 0, End = c(1000, 800)
)

test_that("base-pair breaks and labels remain round and readable", {
  expect_equal(axis_breaks(1000), seq(0, 1000, 200))
  expect_equal(axis_breaks(2.48e8), seq(0, 2e8, 5e7))
  expect_equal(axis_breaks(0), 0)

  breaks <- axis_breaks(5e6)
  unit <- axis_unit(min(diff(breaks)), "auto")
  expect_equal(axis_labels(breaks, unit$div, unit$suffix)[2], "1 Mb")
  expect_equal(
    axis_labels(seq(0, 2e6, 5e5), 1e6, "Mb"),
    c("0.0 Mb", "0.5 Mb", "1.0 Mb", "1.5 Mb", "2.0 Mb")
  )
})

test_that("base-pair axes are additive components with clear outward ticks", {
  base <- ggideogram(axis_kar, show_names = FALSE)
  p <- base + geom_chr_axis(
    chr = "A", side = "left", gap = 0.4,
    tick_length = 2, size = 2.1
  )
  added <- utils::tail(p$layers, 3)

  expect_equal(
    unname(vapply(added, function(layer) class(layer$geom)[1], character(1))),
    c("GeomSegment", "GeomIdeogramTick", "GeomIdeogramAxisText")
  )
  tick <- added[[2]]
  body <- p$coordinates$layout$chrom[1, ]
  expect_equal(body$.x_min - unique(tick$data$x), 0.4)
  expect_true(all(tick$data$nx < 0))
  expect_equal(tick$geom_params$tick_length, 2)
  expect_equal(added[[3]]$aes_params$size, 2.1)
  expect_equal(added[[3]]$geom_params$tick_length, 2)
  expect_equal(added[[3]]$geom_params$label_gap, 0.25)
  expect_equal(added[[3]]$data[c("nx", "ny")], tick$data[c("nx", "ny")])
  expect_no_error(ggplot2::ggplotGrob(p))

  right <- base + geom_chr_axis(chr = "A", side = "right")
  right_tick <- right$layers[[length(right$layers) - 1L]]$data
  expect_true(all(right_tick$nx > 0))
  expect_error(base + geom_chr_axis(chr = "Z"), "Unknown axis chromosome")
})

test_that("horizontal chromosome sides rotate with their tracks and axes", {
  horizontal <- ggideogram(
    axis_kar[1, ], orientation = "horizontal", show_names = FALSE
  )
  right <- horizontal + geom_chr_axis(chr = "A", side = "right")
  left <- horizontal + geom_chr_axis(chr = "A", side = "left")
  right_tick <- right$layers[[length(right$layers) - 1L]]$data
  left_tick <- left$layers[[length(left$layers) - 1L]]$data

  expect_true(all(right_tick$ny < 0))
  expect_true(all(left_tick$ny > 0))
  expect_no_error(ggplot2::ggplotGrob(right))
  expect_no_error(ggplot2::ggplotGrob(left))
})

test_that("axis text clearance follows tick length, text size and direction", {
  base <- ggideogram(axis_kar[1, ], show_names = FALSE)
  short <- base + geom_chr_axis(
    chr = "A", side = "left", tick_length = 0.5,
    label_gap = 0.2, size = 2
  )
  long <- base + geom_chr_axis(
    chr = "A", side = "left", tick_length = 2,
    label_gap = 0.75, size = 3
  )
  short_text <- short$layers[[length(short$layers)]]
  long_text <- long$layers[[length(long$layers)]]

  expect_equal(short_text$geom_params$tick_length, 0.5)
  expect_equal(short_text$geom_params$label_gap, 0.2)
  expect_equal(long_text$geom_params$tick_length, 2)
  expect_equal(long_text$geom_params$label_gap, 0.75)
  expect_true(all(short_text$data$nx < 0))
  expect_true(all(short_text$data$ny == 0))
  expect_no_error(ggplot2::ggplotGrob(short))
  expect_no_error(ggplot2::ggplotGrob(long))
})
