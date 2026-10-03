view_kar <- data.frame(Chr = c("A", "B"), Start = 0, End = c(1000, 800),
                       CE_start = c(400, 300), CE_end = c(600, 450))

test_that("views retain source data and project original bp in both orientations", {
  source <- as_ideogram_data(view_kar)
  view <- chr_view(source, "A", 200, 500)
  expect_identical(view$karyotype, source$karyotype)
  expect_null(source$view)
  for (orientation in c("vertical", "horizontal")) {
    layout <- ideogram_layout(view, orientation = orientation, max_chr_length = 30)
    expect_equal(layout$chrom$.chr, "A")
    expect_equal(layout$chrom$.start, 200)
    expect_equal(layout$chrom$.end, 500)
    expect_equal(layout$chrom$.source_start, 0)
    expect_equal(layout$chrom$.source_end, 1000)
    points <- project_chr_point(layout, data.frame(Chr = "A", Pos = c(200, 350, 500)), "Chr", "Pos")
    long <- if (orientation == "vertical") points$.y else points$.x
    expect_equal(abs(diff(long)), c(15, 15))
    recovered <- unproject_chr_point(layout, points, "Chr", ".x", ".y")
    expect_equal(recovered$.position, c(200, 350, 500))
    p <- ggideogram(view, orientation = orientation, axis = TRUE,
      axis_breaks = c(200, 350, 500), axis_units = "bp")
    labels <- utils::tail(p$layers, 1)[[1]]$data$label
    expect_equal(labels, c("200 bp", "350 bp", "500 bp"))
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that("view selection retains anchors and clips crossing intervals explicitly", {
  view <- chr_view(view_kar, "A", 200, 500)
  points <- data.frame(Chr = c("A", "A", "A", "A", "B"), Pos = c(100, 200, 300, 501, 250))
  expect_equal(view_chr_data(points, view, position = "Pos")$Pos, c(200, 300))
  intervals <- data.frame(Chr = "A", Start = c(100, 450, 100, 501), End = c(250, 700, 199, 800))
  clipped <- view_chr_data(intervals, view)
  expect_equal(clipped$Start, c(201, 450))
  expect_equal(clipped$End, c(250, 500))
  expect_equal(clipped$.source_start, c(100, 450))
  expect_equal(clipped$.source_end, c(250, 700))
  expect_equal(intervals$Start, c(100, 450, 100, 501))
  p <- ggideogram(view) + geom_chr_interval(data = clipped,
    ggplot2::aes(chr = Chr, start = Start, end = End)) +
    geom_chr_marker(data = view_chr_data(points, view, position = "Pos"),
      ggplot2::aes(chr = Chr, position = Pos))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("cut ends are flat while real termini stay rounded", {
  view <- chr_view(view_kar, "A", 0, 300)
  layout <- ideogram_layout(view)
  g <- layout$chrom
  expect_false(g$.cut_start)
  expect_true(g$.cut_end)
  expect_equal(chromosome_halfwidth_next(g, c(0, g$.display_length), 1), c(0, 0.5))
  middle <- ideogram_layout(chr_view(view_kar, "A", 200, 300))$chrom
  expect_equal(chromosome_halfwidth_next(middle, c(0, middle$.display_length), 1), c(0.5, 0.5))
  # A partially visible centromere keeps its original waist position.
  waist <- ideogram_layout(chr_view(view_kar, "A", 450, 650), max_chr_length = 20)$chrom
  expect_equal(chromosome_halfwidth_next(waist, 5, 1), 0)
})

test_that("cytobands are clipped to the selected source interval", {
  bands <- data.frame(Chr = c("A", "A", "B"), Start = c(0, 300, 0),
    End = c(300, 700, 800), Stain = c("gneg", "gpos50", "gneg"))
  source <- as_ideogram_data(view_kar, cytoband = bands)
  view <- chr_view(source, "A", 200, 500)
  layout <- ideogram_layout(view)
  polygons <- cytoband_polygon_data(layout)
  expect_equal(unique(polygons$.chr), "A")
  expect_true(all(polygons$y >= 0 & polygons$y <= 40))
  expect_identical(view$cytoband, source$cytoband)
  expect_no_error(ggplot2::ggplotGrob(ggideogram(view)))
})

test_that("views reject impossible windows and invalid source observations", {
  expect_error(chr_view(view_kar, "Z", 0, 10), "one chromosome")
  expect_error(chr_view(view_kar, "A", 100, 100), "start < end")
  expect_error(chr_view(view_kar, "A", -1, 100), "outside")
  expect_error(chr_view(view_kar, "A", 0, 1001), "outside")
  view <- chr_view(view_kar, "A", 200, 500)
  expect_error(view_chr_data(data.frame(Chr = "A", Pos = 1001), view,
                             position = "Pos"), "outside")
  empty <- view_chr_data(data.frame(Chr = "B", Start = 1, End = 2), view)
  expect_equal(nrow(empty), 0)
})
