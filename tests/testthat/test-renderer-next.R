renderer_kar <- data.frame(
  Chr = c("A", "B"), Start = 0, End = c(1000, 700),
  CE_start = c(400, 250), CE_end = c(600, 400)
)

test_that("ggideogram returns a standard fixed-aspect ggplot", {
  p <- ggideogram(renderer_kar, ncol = 2, base_family = "sans")

  expect_s3_class(p, "ggplot")
  expect_s3_class(p$coordinates$layout, "ideogram_layout_v2")
  expect_equal(p$coordinates$ratio, 1)
  expect_false("ideogram_layout" %in% names(attributes(p)))
  expect_no_error(ggplot2::ggplot_build(p))
  grob <- ggplot2::ggplotGrob(p)
  expect_s3_class(grob, "gtable")
  expect_true(grob$respect)

  themed <- p + ggplot2::labs(title = "Chromosomes") +
    ggplot2::theme(plot.title = ggplot2::element_text(colour = "navy"))
  expect_identical(themed$labels$title, "Chromosomes")
  expect_no_error(ggplot2::ggplot_build(themed))
})

test_that("chromosome polygons follow caps, body width and centromere waist", {
  layout <- ideogram_layout(as_ideogram_data(renderer_kar),
                            chromosome_width = 2, max_chr_length = 20)
  polygon <- chromosome_polygon_data(layout, curve_points = 12)
  a <- polygon[polygon$.chr == "A", ]
  g <- layout$chrom[layout$chrom$.chr == "A", ]

  expect_equal(a$x[1], g$.axis_start_x, tolerance = 1e-12)
  expect_equal(a$y[1], g$.axis_start_y, tolerance = 1e-12)
  expect_equal(max(a$x) - min(a$x), 2, tolerance = 1e-8)
  waist_y <- mean(c(g$.centromere_start_y, g$.centromere_end_y))
  at_waist <- a[abs(a$y - waist_y) < 1e-10, ]
  expect_true(nrow(at_waist) >= 2)
  expect_equal(unique(at_waist$x), g$.axis_start_x, tolerance = 1e-10)
  expect_equal(unname(unlist(a[1, c("x", "y")])),
               unname(unlist(a[nrow(a), c("x", "y")])))
})

test_that("cytobands are clipped to the same chromosome silhouette", {
  bands <- data.frame(
    Chr = rep("A", 3), Start = c(0, 400, 600), End = c(400, 600, 1000),
    Stain = c("gneg", "acen", "gpos50")
  )
  semantic <- as_ideogram_data(renderer_kar[1, ], cytoband = bands)
  layout <- ideogram_layout(semantic, chromosome_width = 2,
                            max_chr_length = 20)
  body <- chromosome_polygon_data(layout, curve_points = 12)
  band <- cytoband_polygon_data(layout, curve_points = 12)

  expect_setequal(unique(band$.fill),
                  cytoband_colours(bands$Stain))
  expect_gte(min(band$x), min(body$x) - 1e-10)
  expect_lte(max(band$x), max(body$x) + 1e-10)
  expect_gte(min(band$y), min(body$y) - 1e-10)
  expect_lte(max(band$y), max(body$y) + 1e-10)

  p <- ggideogram(renderer_kar[1, ], cytoband = bands,
                  show_names = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(built$data), 3L)
  expect_setequal(unique(built$data[[2]]$fill),
                  cytoband_colours(bands$Stain))
})

test_that("chromosome names use standard text size and typographic gap", {
  p1 <- ggideogram(renderer_kar, name_size = 2.8, name_gap = 0.25)
  p2 <- ggideogram(renderer_kar, name_size = 2.8, name_gap = 1)
  name1 <- p1$layers[[length(p1$layers)]]
  name2 <- p2$layers[[length(p2$layers)]]

  expect_s3_class(name1$geom, "GeomText")
  expect_equal(name1$aes_params$size, 2.8)
  expect_equal(name1$data$y, name2$data$y)
  expect_gt(name2$aes_params$vjust, name1$aes_params$vjust)
})

test_that("local chromosome axes use physical ticks and standard text", {
  p <- ggideogram(
    renderer_kar, axis = "A", axis_side = "left",
    axis_tick_length = 2.25, axis_size = 2.1,
    axis_linewidth = 0.35, show_names = FALSE
  )
  classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  tick_index <- which(classes == "GeomIdeogramTick")
  text_index <- which(classes == "GeomIdeogramAxisText")

  expect_length(tick_index, 1L)
  expect_equal(formals(ggideogram)$axis_gap, 0.3)
  expect_equal(p$layers[[tick_index]]$geom_params$tick_length, 2.25)
  expect_equal(p$layers[[tick_index]]$aes_params$linewidth, 0.35)
  tick_data <- p$layers[[tick_index]]$data
  body <- p$coordinates$layout$chrom[
    p$coordinates$layout$chrom$.chr == "A", ]
  expect_equal(body$.x_min - unique(tick_data$x), 0.3)
  expect_true(all(tick_data$nx < 0))
  expect_length(text_index, 1L)
  expect_equal(p$layers[[text_index]]$aes_params$size, 2.1)
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("horizontal orientation renders the same semantic geometry", {
  vertical <- ggideogram(renderer_kar, orientation = "vertical",
                         show_names = FALSE)
  horizontal <- ggideogram(renderer_kar, orientation = "horizontal",
                           show_names = FALSE)
  bv <- ggplot2::ggplot_build(vertical)$data[[1]]
  bh <- ggplot2::ggplot_build(horizontal)$data[[1]]

  expect_equal(range(bh$x),
               diff(vertical$coordinates$layout$bounds$y) - rev(range(bv$y)))
  expect_no_error(ggplot2::ggplotGrob(horizontal))
})

test_that("new renderer rejects geometry that cannot form a chromosome body", {
  extreme <- data.frame(Chr = c("long", "tiny"), Start = 0,
                        End = c(1e6, 1))
  expect_error(ggideogram(extreme, max_chr_length = 20),
               "shorter than `chromosome_width`")
  expect_no_error(ggideogram(extreme, max_chr_length = 20,
                             scale_length = "per_chr"))
})

test_that("semantic input is not ambiguously remapped", {
  semantic <- as_ideogram_data(renderer_kar)
  expect_no_error(ggideogram(semantic))
  expect_error(
    ggideogram(semantic,
               mapping = ggplot2::aes(chr = .data$Chr,
                                      start = .data$Start, end = .data$End)),
    "cannot also supply"
  )
})
