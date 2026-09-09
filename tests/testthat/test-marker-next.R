marker_kar <- data.frame(
  Chr = c("A", "B"), Start = 0, End = c(1000, 700),
  CE_start = c(400, 250), CE_end = c(600, 400)
)

marker_data <- data.frame(
  Chr = c("A", "A", "A", "B"),
  Pos = c(300, 500, 700, 350),
  Type = c("circle", "box", "triangle", "circle"),
  Score = c(1, 2, 4, 3),
  stringsAsFactors = FALSE
)

test_that("chromosome markers are genuine GeomPoint layers", {
  layer <- geom_chr_marker(
    data = marker_data,
    ggplot2::aes(chr = Chr, position = Pos)
  )

  expect_s3_class(layer, "Layer")
  expect_s3_class(layer$geom, "GeomPoint")
  expect_false(inherits(layer$geom, "GeomIdeogramMarker"))
  expect_equal(layer$stat_params$gap, 0.25)

  p <- ggideogram(marker_kar, show_names = FALSE) + layer
  expect_no_error(ggplot2::ggplot_build(p))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("point aesthetics use native scales and guides", {
  p <- ggideogram(marker_kar, show_names = FALSE) +
    geom_chr_marker(
      data = marker_data,
      ggplot2::aes(
        chr = Chr, position = Pos, shape = Type,
        fill = Type, colour = Type, size = Score
      ),
      alpha = 0.8, stroke = 0.35
    ) +
    ggplot2::scale_shape_manual(
      values = c(circle = 21, box = 22, triangle = 24)
    ) +
    ggplot2::scale_size_continuous(range = c(0.8, 2.4)) +
    ggplot2::labs(shape = "Marker type", fill = "Marker type",
                  colour = "Marker type", size = "Score")

  built <- ggplot2::ggplot_build(p)
  point <- built$data[[length(built$data)]]
  expect_setequal(point$shape, c(21, 22, 24))
  expect_equal(range(point$size), c(0.8, 2.4))
  expect_true(all(point$alpha == 0.8))
  expect_true(all(point$stroke == 0.35))
  expect_true(length(built$plot$guides$aesthetics) >= 1L)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("constant marker size stays in ordinary ggplot2 units", {
  small <- ggideogram(marker_kar, show_names = FALSE) +
    geom_chr_marker(
      data = marker_data,
      ggplot2::aes(chr = Chr, position = Pos),
      size = 0.9
    )
  large <- ggideogram(marker_kar, show_names = FALSE) +
    geom_chr_marker(
      data = marker_data,
      ggplot2::aes(chr = Chr, position = Pos),
      size = 2.1
    )
  small_data <- utils::tail(ggplot2::ggplot_build(small)$data, 1)[[1]]
  large_data <- utils::tail(ggplot2::ggplot_build(large)$data, 1)[[1]]

  expect_true(all(small_data$size == 0.9))
  expect_true(all(large_data$size == 2.1))
  expect_false(any(c("size_px", "canvas") %in%
                     names(small$layers[[length(small$layers)]]$geom_params)))
})

test_that("marker side and gap are projected by the shared coordinate", {
  right <- ggideogram(marker_kar[1, ], show_names = FALSE) +
    geom_chr_marker(
      data = marker_data[1, ],
      ggplot2::aes(chr = Chr, position = Pos),
      side = "right", gap = 0.2
    )
  left <- ggideogram(marker_kar[1, ], show_names = FALSE) +
    geom_chr_marker(
      data = marker_data[1, ],
      ggplot2::aes(chr = Chr, position = Pos),
      side = "left", gap = 0.2
    )

  right_built <- ggplot2::ggplot_build(right)
  left_built <- ggplot2::ggplot_build(left)
  right_data <- utils::tail(right_built$data, 1)[[1]]
  left_data <- utils::tail(left_built$data, 1)[[1]]
  right_xy <- right$coordinates$transform(
    right_data, right_built$layout$panel_params[[1]])
  left_xy <- left$coordinates$transform(
    left_data, left_built$layout$panel_params[[1]])
  centre <- project_chr_point(
    right$coordinates$layout,
    data.frame(Chr = "A", Pos = 300), "Chr", "Pos")

  expect_gt(right_xy$x, 0.5)
  expect_lt(left_xy$x, 0.5)
  expect_equal(right_xy$y, left_xy$y)
  # Each side now reserves its own marker extent in the panel bounds.
  right_range <- right_built$layout$panel_params[[1]]$x.range
  left_range <- left_built$layout$panel_params[[1]]$x.range
  expect_equal(
    right_xy$x * diff(right_range) + right_range[1] -
      (left_xy$x * diff(left_range) + left_range[1]),
    2 * (right$coordinates$layout$chromosome_width / 2 + 0.2),
    tolerance = 1e-8
  )
  expect_equal(centre$.position, 300)
})

test_that("horizontal marker sides are a deterministic rotation", {
  right <- ggideogram(
    marker_kar[1, ], orientation = "horizontal", show_names = FALSE
  ) + geom_chr_marker(
    data = marker_data[1, ],
    ggplot2::aes(chr = Chr, position = Pos),
    side = "right", gap = 0.2
  )
  left <- ggideogram(
    marker_kar[1, ], orientation = "horizontal", show_names = FALSE
  ) + geom_chr_marker(
    data = marker_data[1, ],
    ggplot2::aes(chr = Chr, position = Pos),
    side = "left", gap = 0.2
  )
  right_built <- ggplot2::ggplot_build(right)
  left_built <- ggplot2::ggplot_build(left)
  right_xy <- right$coordinates$transform(
    utils::tail(right_built$data, 1)[[1]],
    right_built$layout$panel_params[[1]]
  )
  left_xy <- left$coordinates$transform(
    utils::tail(left_built$data, 1)[[1]],
    left_built$layout$panel_params[[1]]
  )

  expect_lt(right_xy$y, 0.5)
  expect_gt(left_xy$y, 0.5)
  expect_equal(right_xy$x, left_xy$x)
})

test_that("chromosome links are genuine GeomSegment layers", {
  layer <- geom_chr_link(
    data = marker_data,
    ggplot2::aes(chr = Chr, position = Pos),
    linewidth = 0.3
  )
  expect_s3_class(layer$geom, "GeomSegment")
  expect_false(layer$show.legend)

  p <- ggideogram(marker_kar, show_names = FALSE) + layer
  built <- ggplot2::ggplot_build(p)
  segment <- utils::tail(built$data, 1)[[1]]
  transformed <- p$coordinates$transform(
    segment, built$layout$panel_params[[1]])
  expect_true(all(transformed$xend > transformed$x))
  expect_equal(transformed$yend, transformed$y)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("position_chr_repel preserves anchors and bounds display positions", {
  dense <- data.frame(
    Chr = rep("A", 3), Pos = c(495, 500, 505), id = 1:3
  )
  repel <- position_chr_repel(100, units = "bp")
  p <- ggideogram(marker_kar[1, ], show_names = FALSE) +
    geom_chr_marker(
      data = dense,
      ggplot2::aes(chr = Chr, position = Pos, group = id),
      position = repel
    )
  built <- ggplot2::ggplot_build(p)
  point <- utils::tail(built$data, 1)[[1]]

  expect_equal(point$ideogram_anchor_position, dense$Pos)
  expect_equal(point$ideogram_position, dense$Pos)
  transformed_semantic <- apply_chr_repel_request(
    p$coordinates$layout, point)
  expect_equal(transformed_semantic$ideogram_anchor_position, dense$Pos)
  expect_gte(min(transformed_semantic$ideogram_position), 0)
  expect_lte(max(transformed_semantic$ideogram_position), 1000)
  expect_gte(min(diff(sort(transformed_semantic$ideogram_position))), 100)
})

test_that("repelled marker and link share display endpoints", {
  dense <- data.frame(Chr = rep("A", 3), Pos = c(495, 500, 505))
  repel <- position_chr_repel(0.8, units = "body_width")
  p <- ggideogram(marker_kar[1, ], show_names = FALSE) +
    geom_chr_link(
      data = dense, ggplot2::aes(chr = Chr, position = Pos),
      position = repel, gap = 0.3
    ) +
    geom_chr_marker(
      data = dense, ggplot2::aes(chr = Chr, position = Pos),
      position = repel, gap = 0.3, size = 1
    )
  built <- ggplot2::ggplot_build(p)
  link <- built$data[[length(built$data) - 1L]]
  point <- built$data[[length(built$data)]]
  link_xy <- p$coordinates$transform(
    link, built$layout$panel_params[[1]])
  point_xy <- p$coordinates$transform(
    point, built$layout$panel_params[[1]])

  expect_equal(link_xy$xend, point_xy$x)
  expect_equal(link_xy$yend, point_xy$y)
  expect_false(isTRUE(all.equal(link_xy$y, link_xy$yend)))
})

test_that("marker constructors reject ambiguous layout parameters", {
  expect_error(position_chr_repel(0), "positive")
  expect_error(position_chr_repel(1, units = "pixel"), "arg")
  expect_error(geom_chr_marker(gap = -1), "non-negative")
  expect_error(geom_chr_link(side = "top"), "arg")
})
