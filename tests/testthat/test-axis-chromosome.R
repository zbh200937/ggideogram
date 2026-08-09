axis_karyotype <- data.frame(
  Chr = c("A", "B", "C"),
  Start = c(0, 0, 0),
  End = c(100, 80, 50),
  CE_start = c(40, 28, 20),
  CE_end = c(60, 36, 30)
)

test_that("chromosome scales are native discrete position scales", {
  x <- scale_x_chromosome(axis_karyotype)
  y <- scale_y_chromosome(axis_karyotype)

  expect_s3_class(x, "ScaleDiscretePosition")
  expect_s3_class(y, "ScaleDiscretePosition")
  expect_s3_class(x$guide, "GuideChromosomeAxis")
  expect_s3_class(y$guide, "GuideChromosomeAxis")
  expect_identical(x$get_limits(), axis_karyotype$Chr)
  expect_identical(y$get_limits(), axis_karyotype$Chr)
})

test_that("native host geoms are not rewritten by a chromosome axis", {
  values <- data.frame(
    Chr = rep(axis_karyotype$Chr, each = 3),
    Value = seq_len(9)
  )
  columns <- ggplot2::ggplot(values, ggplot2::aes(.data$Chr, .data$Value)) +
    ggplot2::geom_col() +
    scale_x_chromosome(axis_karyotype)
  boxes <- ggplot2::ggplot(values, ggplot2::aes(.data$Value, .data$Chr)) +
    ggplot2::geom_boxplot() +
    scale_y_chromosome(axis_karyotype)
  line <- ggplot2::ggplot(
    values[seq(1, nrow(values), by = 3), ],
    ggplot2::aes(.data$Chr, .data$Value, group = 1)
  ) +
    ggplot2::geom_line() +
    scale_x_chromosome(axis_karyotype)

  expect_true(inherits(columns$layers[[1]]$geom, "GeomBar"))
  expect_true(inherits(boxes$layers[[1]]$geom, "GeomBoxplot"))
  expect_true(inherits(line$layers[[1]]$geom, "GeomLine"))
  expect_s3_class(ggplot2::ggplotGrob(columns), "gtable")
  expect_s3_class(ggplot2::ggplotGrob(boxes), "gtable")
  expect_s3_class(ggplot2::ggplotGrob(line), "gtable")
})

test_that("chromosome guides build at all four standard positions", {
  values <- data.frame(Chr = axis_karyotype$Chr, Value = c(2, 5, 3))

  for (position in c("bottom", "top")) {
    plot <- ggplot2::ggplot(
      values, ggplot2::aes(.data$Chr, .data$Value)) +
      ggplot2::geom_col() +
      scale_x_chromosome(axis_karyotype, position = position)
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  }
  for (position in c("left", "right")) {
    plot <- ggplot2::ggplot(
      values, ggplot2::aes(.data$Value, .data$Chr)) +
      ggplot2::geom_col() +
      scale_y_chromosome(axis_karyotype, position = position)
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  }
})

test_that("raw breaks own alignment while formatted labels remain free", {
  values <- data.frame(Chr = axis_karyotype$Chr, Value = c(2, 5, 3))
  plot <- ggplot2::ggplot(
    values, ggplot2::aes(.data$Chr, .data$Value)) +
    ggplot2::geom_col() +
    scale_x_chromosome(
      axis_karyotype,
      limits = c("C", "A", "B"),
      labels = function(chr) paste("Chromosome", chr)
    )
  built <- ggplot2::ggplot_build(plot)

  expect_identical(
    built$layout$panel_scales_x[[1]]$get_limits(),
    c("C", "A", "B")
  )
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
})

test_that("data-trained limits and unknown chromosome breaks are explicit", {
  values <- data.frame(Chr = c("A", "C"), Value = c(2, 3))
  plot <- ggplot2::ggplot(
    values, ggplot2::aes(.data$Chr, .data$Value)) +
    ggplot2::geom_col() +
    scale_x_chromosome(axis_karyotype, limits = NULL, drop = TRUE)
  built <- ggplot2::ggplot_build(plot)
  expect_identical(built$layout$panel_scales_x[[1]]$get_limits(), c("A", "C"))

  unknown <- rbind(values, data.frame(Chr = "Z", Value = 1))
  bad_plot <- ggplot2::ggplot(
    unknown, ggplot2::aes(.data$Chr, .data$Value)) +
    ggplot2::geom_col() +
    scale_x_chromosome(axis_karyotype, limits = NULL)
  expect_error(
    ggplot2::ggplotGrob(bad_plot),
    "Chromosome-axis break.*not found"
  )
})

test_that("axis chromosome length and centromere geometry are data-derived", {
  semantic <- as_ideogram_data(axis_karyotype)
  lengths <- semantic$karyotype$.end - semantic$karyotype$.start
  ratios <- lengths / max(lengths)
  expect_equal(ratios, c(1, 0.8, 0.5))

  profile <- chromosome_axis_profile(
    semantic$karyotype[1, , drop = FALSE],
    body_length_mm = 20,
    body_width_mm = 2,
    curve_points = 16
  )
  expect_equal(range(profile$long), c(0, 1))
  expect_true(all(profile$cross >= 0 & profile$cross <= 1))
  waist <- (axis_karyotype$CE_start[1] + axis_karyotype$CE_end[1]) /
    (2 * axis_karyotype$End[1])
  expect_true(any(abs(profile$long - waist) < 1e-12 &
                    abs(profile$cross - 0.5) < 1e-12))
})

test_that("axis physical dimensions reject relative and invalid units", {
  expect_error(
    guide_chromosome_axis(axis_karyotype, axis_length = grid::unit(1, "npc")),
    "absolute physical unit"
  )
  expect_error(
    guide_chromosome_axis(axis_karyotype, body_width = grid::unit(0, "mm")),
    "positive finite"
  )
  expect_error(
    guide_chromosome_axis(axis_karyotype, label_gap = grid::unit(-1, "mm")),
    "non-negative finite"
  )
})
