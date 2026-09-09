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

test_that("bp axes are outside same-side tracks and include their spine in bounds", {
  for (orientation in c("vertical", "horizontal")) {
    for (side in c("left", "right")) {
      p <- ggideogram(axis_kar[1, ], orientation = orientation,
        tracks = track_layout(signal = track(side, width = 1.1, gap = 0.2)),
        axis = "A", axis_side = side, axis_gap = 0.3, show_names = FALSE)
      layout <- p$coordinates$layout
      tick <- p$layers[[length(p$layers) - 1L]]$data
      g <- layout$chrom[1, ]
      distance <- sqrt((tick$x[1] - g$.axis_start_x)^2 +
        (tick$y[1] - g$.axis_start_y)^2)
      expect_equal(distance, layout$chromosome_width / 2 + 1.1 + 0.2 + 0.3)
      expect_true(all(tick$x >= layout$bounds$x[1] & tick$x <= layout$bounds$x[2]))
      expect_true(all(tick$y >= layout$bounds$y[1] & tick$y <= layout$bounds$y[2]))
      expect_no_error(ggplot2::ggplotGrob(p))
    }
  }
})

test_that("direct marker lanes displace bp axes independent of addition order", {
  markers <- data.frame(Chr = "A", Pos = 300)
  for (orientation in c("vertical", "horizontal")) {
    for (side in c("left", "right")) {
      base <- ggideogram(axis_kar[1, ], orientation = orientation, show_names = FALSE)
      marker <- geom_chr_marker(data = markers,
        ggplot2::aes(chr = Chr, position = Pos), side = side, gap = 0.5, size = 1.8)
      axis <- geom_chr_axis(chr = "A", side = side, gap = 0.3)
      before <- base + axis
      first <- before + marker
      last <- base + marker + axis
      get_tick <- function(plot) Filter(function(layer)
        inherits(layer$geom, "GeomIdeogramTick"), plot$layers)[[1]]$data
      expect_equal(get_tick(first), get_tick(last))
      expect_equal(first$coordinates$layout$bounds, last$coordinates$layout$bounds)
      delta <- get_tick(first)[1, c("x", "y")] - get_tick(before)[1, c("x", "y")]
      expect_equal(sqrt(sum(delta^2)), 1)
      expect_no_error(ggplot2::ggplotGrob(first))
      expect_no_error(ggplot2::ggplotGrob(last))
      expect_null(before$coordinates$layout$marker_extent)
      expect_s3_class(utils::tail(first$layers, 1)[[1]]$geom, "GeomPoint")
    }
  }
})

test_that("long perimeter labels reserve physical space and remain theme-overridable", {
  p <- ggideogram(axis_kar[1, ], axis = TRUE, axis_side = "right",
    axis_breaks = c(0, 1000), axis_units = "bp", axis_size = 4,
    axis_tick_length = 3, axis_label_gap = 0.5)
  margin <- grid::convertUnit(p$theme$plot.margin, "mm", valueOnly = TRUE)
  text <- grid::textGrob("1000 bp", gp = grid::gpar(fontsize = 4 * ggplot2::.pt))
  expect_gte(margin[2], grid::convertWidth(grid::grobWidth(text), "mm", valueOnly = TRUE) + 5)
  overridden <- p + ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10))
  expect_equal(overridden$theme$plot.margin, ggplot2::margin(10, 10, 10, 10))
})
