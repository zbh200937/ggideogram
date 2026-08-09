text_kar <- data.frame(Chr = "A", Start = 0, End = 1000)
text_data <- data.frame(
  Chr = "A", Position = c(200, 500, 800),
  Label = c("gene1", "gene2", "gene3"),
  Class = c("one", "two", "one")
)

test_that("chromosome text is a standard GeomText layer", {
  p <- ggideogram(text_kar, show_names = FALSE) +
    geom_chr_text(
      data = text_data,
      ggplot2::aes(
        chr = Chr, position = Position,
        label = Label, colour = Class
      ),
      size = 2.3,
      hjust = 0
    )
  layer <- p$layers[[length(p$layers)]]
  built <- ggplot2::ggplot_build(p)
  text <- utils::tail(built$data, 1)[[1]]

  expect_s3_class(layer$geom, "GeomText")
  expect_equal(text$label, text_data$Label)
  expect_true(all(text$size == 2.3))
  expect_true(length(built$plot$guides$aesthetics) >= 1L)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("text can align to a declared track and share marker links", {
  tracks <- track_layout(labels = track("right", width = 1.2, gap = 0.3))
  repel <- position_chr_repel(0.8, units = "body_width")
  p <- ggideogram(text_kar, tracks = tracks, show_names = FALSE) +
    geom_chr_link(
      data = text_data,
      ggplot2::aes(chr = Chr, position = Position),
      position = repel, track = "labels"
    ) +
    geom_chr_text(
      data = text_data,
      ggplot2::aes(chr = Chr, position = Position, label = Label),
      position = repel, track = "labels", hjust = 0, size = 2
    )
  built <- ggplot2::ggplot_build(p)
  link <- built$data[[length(built$data) - 1L]]
  text <- built$data[[length(built$data)]]
  link_xy <- p$coordinates$transform(
    link, built$layout$panel_params[[1]]
  )
  text_xy <- p$coordinates$transform(
    text, built$layout$panel_params[[1]]
  )

  expect_equal(link_xy$xend, text_xy$x)
  expect_equal(link_xy$yend, text_xy$y)
  expect_no_error(ggplot2::ggplotGrob(p))
})
