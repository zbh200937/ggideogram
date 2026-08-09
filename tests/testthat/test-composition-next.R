composition_kar <- data.frame(
  Chr = c("A", "B", "C"),
  Start = 0,
  End = c(100, 90, 70),
  CE_start = c(40, 30, 25),
  CE_end = c(55, 45, 35)
)

composition_signal <- data.frame(
  Chr = rep(composition_kar$Chr, each = 5),
  Position = rep(seq(10, 60, length.out = 5), 3),
  Value = rep(c(1, 3, 2, 5, 4), 3)
)

make_composition_ideogram <- function() {
  tracks <- track_layout(
    signal = track("right", width = 1.2, limits = c(0, 5))
  )
  markers <- data.frame(
    Chr = composition_kar$Chr,
    Position = c(30, 40, 50),
    Type = c("one", "two", "one")
  )
  ggideogram(
    composition_kar, ncol = 3, tracks = tracks,
    axis = "A", show_names = TRUE
  ) +
    geom_track_line(
      data = composition_signal,
      ggplot2::aes(
        chr = Chr, position = Position, value = Value, group = Chr
      ),
      track = "signal", linewidth = 0.35
    ) +
    geom_chr_marker(
      data = markers,
      ggplot2::aes(
        chr = Chr, position = Position, shape = Type, fill = Type
      ),
      size = 0.8
    ) +
    ggplot2::scale_shape_manual(values = c(one = 21, two = 24)) +
    ggplot2::labs(shape = "Marker", fill = "Marker")
}

test_that("ideogram is exactly a thin ggideogram preset", {
  direct <- ggideogram(composition_kar, ncol = 3)
  preset <- ideogram(composition_kar, ncol = 3)

  expect_identical(names(formals(ideogram)), "...")
  expect_equal(
    vapply(direct$layers, function(layer) class(layer$geom)[1], character(1)),
    vapply(preset$layers, function(layer) class(layer$geom)[1], character(1))
  )
  expect_equal(
    direct$coordinates$layout$chrom,
    preset$coordinates$layout$chrom
  )
  expect_false("ideogram_layout" %in% names(attributes(preset)))
})

test_that("standard plot composition keeps semantic and physical components", {
  skip_if_not_installed("patchwork")
  p <- make_composition_ideogram()
  ordinary <- ggplot2::ggplot(
    data.frame(group = c("a", "b"), value = c(2, 5)),
    ggplot2::aes(group, value, fill = group)
  ) + ggplot2::geom_col() + ggplot2::theme_minimal()

  across <- patchwork::wrap_plots(p, ordinary, widths = c(1.4, 1))
  stacked <- patchwork::wrap_plots(p, ordinary, ncol = 1)
  inset <- ordinary + patchwork::inset_element(
    p, left = 0.48, bottom = 0.42, right = 0.98, top = 0.98
  )

  expect_s3_class(patchwork::patchworkGrob(across), "gtable")
  expect_s3_class(patchwork::patchworkGrob(stacked), "gtable")
  expect_s3_class(patchwork::patchworkGrob(inset), "gtable")

  built <- ggplot2::ggplot_build(p)
  marker <- utils::tail(built$data, 1)[[1]]
  expect_true(all(marker$size == 0.8))
  expect_equal(
    p$layers[[which(vapply(
      p$layers, function(layer) inherits(layer$geom, "GeomText"),
      logical(1)
    ))[1]]]$aes_params$size,
    3.2
  )
})

test_that("cowplot, ggplotGrob and ggsave use the standard ggplot boundary", {
  skip_if_not_installed("cowplot")
  p <- make_composition_ideogram()
  ordinary <- ggplot2::ggplot(
    data.frame(x = 1:4, y = c(1, 4, 2, 3)), ggplot2::aes(x, y)
  ) + ggplot2::geom_line() + ggplot2::theme_minimal()
  combined <- cowplot::ggdraw(ordinary) +
    cowplot::draw_plot(p, x = 0.52, y = 0.48, width = 0.46, height = 0.48)

  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  expect_s3_class(ggplot2::ggplotGrob(combined), "gtable")

  wide <- tempfile(fileext = ".png")
  tall <- tempfile(fileext = ".png")
  on.exit(unlink(c(wide, tall)), add = TRUE)
  expect_no_error(ggplot2::ggsave(
    wide, p, width = 160, height = 70, units = "mm", dpi = 72
  ))
  expect_no_error(ggplot2::ggsave(
    tall, p, width = 80, height = 140, units = "mm", dpi = 72
  ))
  expect_gt(file.info(wide)$size, 0)
  expect_gt(file.info(tall)$size, 0)
})
