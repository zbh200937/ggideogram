interval_kar <- data.frame(
  Chr = c("A", "B"), Start = 0, End = c(1000, 800)
)
interval_data <- data.frame(
  Chr = c("A", "B"), Start = c(200, 100), End = c(450, 300),
  Class = c("qtl", "cnv")
)

test_that("chromosome intervals use standard GeomSegment semantics", {
  p <- ggideogram(interval_kar, show_names = FALSE) +
    geom_chr_interval(
      data = interval_data,
      ggplot2::aes(
        chr = Chr, start = Start, end = End, colour = Class
      ),
      linewidth = 1.1
    )
  layer <- p$layers[[length(p$layers)]]
  built <- ggplot2::ggplot_build(p)
  data <- utils::tail(built$data, 1)[[1]]
  coordinates <- p$coordinates$transform(
    data, built$layout$panel_params[[1]]
  )
  projected <- project_chr_interval(
    p$coordinates$layout, interval_data, "Chr", "Start", "End"
  )
  projected_start <- p$coordinates$transform(
    data.frame(x = projected$.x_start, y = projected$.y_start),
    built$layout$panel_params[[1]]
  )
  projected_end <- p$coordinates$transform(
    data.frame(x = projected$.x_end, y = projected$.y_end),
    built$layout$panel_params[[1]]
  )

  expect_s3_class(layer$geom, "GeomSegment")
  expect_equal(data$linewidth, c(1.1, 1.1))
  expect_equal(coordinates$x, projected_start$x)
  expect_equal(coordinates$y, projected_start$y)
  expect_equal(coordinates$xend, projected_end$x)
  expect_equal(coordinates$yend, projected_end$y)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("intervals can occupy declared beside tracks", {
  tracks <- track_layout(regions = track("right", width = 0.8, gap = 0.25))
  p <- ggideogram(
    interval_kar, tracks = tracks, show_names = FALSE
  ) + geom_chr_region(
    data = interval_data,
    ggplot2::aes(chr = Chr, start = Start, end = End),
    placement = "beside", track = "regions"
  )
  built <- ggplot2::ggplot_build(p)
  data <- utils::tail(built$data, 1)[[1]]
  coordinates <- p$coordinates$transform(
    data, built$layout$panel_params[[1]]
  )
  centre <- project_chr_point(
    p$coordinates$layout,
    data.frame(Chr = "A", Position = 200), "Chr", "Position"
  )
  centre <- p$coordinates$transform(
    data.frame(x = centre$.x, y = centre$.y),
    built$layout$panel_params[[1]]
  )

  expect_gt(coordinates$x[1], centre$x[1])
  expect_no_error(ggplot2::ggplotGrob(p))
  expect_error(
    geom_chr_interval(track = "regions"),
    "placement"
  )
  expect_error(
    ggplot2::ggplot_build(
      ggideogram(interval_kar) + geom_chr_interval(
        data = transform(interval_data[1, ], Start = 500, End = 100),
        ggplot2::aes(chr = Chr, start = Start, end = End)
      )
    ),
    "start <= end"
  )
})
