next_kar <- data.frame(
  Chr = c("A", "B", "C", "D"), Start = c(0, 10, 0, 100),
  End = c(1000, 510, 800, 500),
  CE_start = c(400, 200, 300, 250),
  CE_end = c(600, 300, 500, 350)
)

test_that("ideogram_data dispatches to a dimensionless layout", {
  layout <- ideogram_layout(as_ideogram_data(next_kar))

  expect_s3_class(layout, "ideogram_layout_v2")
  expect_identical(layout$units, "chromosome_width")
  expect_equal(layout$chromosome_width, 1)
  expect_equal(max(layout$chrom$.display_length), 40)
  expect_true(layout$length_comparable)
  expect_false(any(c("mpx", "px_per_bp", "flip_y", "chr_width") %in%
                     names(layout)))
  expect_false(any(grepl("canvas|dpi|a4|pixel|mm", names(layout),
                         ignore.case = TRUE)))
})

test_that("data frames and semantic inputs use the same single layout engine", {
  direct <- ideogram_layout(next_kar)
  semantic <- ideogram_layout(as_ideogram_data(next_kar))

  expect_s3_class(direct, "ideogram_layout_v2")
  expect_equal(direct$chrom, semantic$chrom)
  expect_equal(direct$bounds, semantic$bounds)
  expect_false(any(c("mpx", "bp2y", "chr_width") %in% names(direct)))
})

test_that("global length scale preserves chromosome comparisons", {
  layout <- ideogram_layout(as_ideogram_data(next_kar), max_chr_length = 25)
  expected <- (next_kar$End - next_kar$Start) / 1000 * 25

  expect_equal(layout$chrom$.display_length, expected)
  expect_equal(layout$chrom$.units_per_bp,
               rep(25 / 1000, nrow(next_kar)))
  expect_equal(layout$chrom$.axis_start_y - layout$chrom$.axis_end_y,
               expected)
})

test_that("multi-row layout is aligned and has only declared gaps", {
  layout <- ideogram_layout(
    as_ideogram_data(next_kar), ncol = 2,
    chromosome_width = 2, chromosome_gap = 0.5,
    max_chr_length = 20, row_gap = 3
  )

  expect_equal(layout$nrow, 2L)
  expect_equal(layout$ncol, 2L)
  expect_equal(layout$chrom$.row, c(0, 0, 1, 1))
  expect_equal(layout$chrom$.col, c(0, 1, 0, 1))
  expect_equal(diff(layout$chrom$.axis_start_x[1:2]), 2.5)
  expect_equal(layout$chrom$.axis_start_x[1], layout$chrom$.axis_start_x[3])

  first_bottom <- unique(layout$chrom$.axis_end_y[layout$chrom$.row == 0])
  second_top <- max(layout$chrom$.axis_start_y[layout$chrom$.row == 1])
  expect_equal(first_bottom - second_top, 3)
  expect_equal(layout$bounds$x, c(0, 4.5))
  expect_equal(diff(layout$bounds$y), sum(layout$row_height) + 3)
})

test_that("point projection is linear, bounded and reversible", {
  layout <- ideogram_layout(as_ideogram_data(next_kar), max_chr_length = 40)
  points <- data.frame(chromosome = c("A", "A", "B"),
                       bp = c(0, 500, 510), id = letters[1:3])
  projected <- project_chr_point(layout, points, "chromosome", "bp")

  a <- layout$chrom[layout$chrom$.chr == "A", ]
  b <- layout$chrom[layout$chrom$.chr == "B", ]
  expect_equal(projected$.x[1:2], rep(a$.axis_start_x, 2))
  expect_equal(projected$.y[1], a$.axis_start_y)
  expect_equal(projected$.y[2], mean(c(a$.axis_start_y, a$.axis_end_y)))
  expect_equal(projected$.y[3], b$.axis_end_y)
  expect_identical(projected$id, points$id)

  recovered <- unproject_chr_point(
    layout,
    transform(projected, plot_x = .x + c(0.4, -0.2, 0.7), plot_y = .y),
    "chromosome", "plot_x", "plot_y"
  )
  expect_equal(recovered$.position, points$bp)

  expect_error(project_chr_point(
    layout, data.frame(chr = "Z", pos = 1), "chr", "pos"),
    "Unknown chromosome")
  expect_error(project_chr_point(
    layout, data.frame(chr = "A", pos = 1001), "chr", "pos"),
    "outside its chromosome")
})

test_that("interval projection keeps biological endpoints", {
  layout <- ideogram_layout(as_ideogram_data(next_kar))
  intervals <- data.frame(chr = c("A", "B"), lo = c(100, 10),
                          hi = c(300, 510), label = c("x", "y"))
  projected <- project_chr_interval(layout, intervals, "chr", "lo", "hi")
  starts <- project_chr_point(
    layout, data.frame(chr = intervals$chr, pos = intervals$lo), "chr", "pos")
  ends <- project_chr_point(
    layout, data.frame(chr = intervals$chr, pos = intervals$hi), "chr", "pos")

  expect_equal(projected$.x_start, starts$.x)
  expect_equal(projected$.y_start, starts$.y)
  expect_equal(projected$.x_end, ends$.x)
  expect_equal(projected$.y_end, ends$.y)
  expect_identical(projected$label, intervals$label)
  expect_error(project_chr_interval(
    layout, data.frame(chr = "A", lo = 20, hi = 10), "chr", "lo", "hi"),
    "start <= end", fixed = TRUE)
})

test_that("horizontal orientation is a deterministic rotation", {
  data <- as_ideogram_data(next_kar)
  vertical <- ideogram_layout(data, ncol = 2, max_chr_length = 20)
  horizontal <- ideogram_layout(data, ncol = 2, max_chr_length = 20,
                                orientation = "horizontal")
  point <- data.frame(chr = c("A", "B", "C", "D"), pos = next_kar$Start)
  pv <- project_chr_point(vertical, point, "chr", "pos")
  ph <- project_chr_point(horizontal, point, "chr", "pos")

  expect_equal(ph$.x, diff(vertical$bounds$y) - pv$.y)
  expect_equal(ph$.y, diff(vertical$bounds$x) - pv$.x)
  expect_equal(horizontal$bounds$x, vertical$bounds$y)
  expect_equal(horizontal$bounds$y, vertical$bounds$x)
})

test_that("per-chromosome scaling is explicit and marked non-comparable", {
  layout <- ideogram_layout(as_ideogram_data(next_kar),
                            scale_length = "per_chr", max_chr_length = 17)
  expect_equal(layout$chrom$.display_length, rep(17, nrow(next_kar)))
  expect_false(layout$length_comparable)
  expect_identical(layout$scale_length, "per_chr")
})

test_that("dimensionless layout parameters are validated", {
  data <- as_ideogram_data(next_kar)
  expect_error(ideogram_layout(data, ncol = 1.5), "whole number")
  expect_error(ideogram_layout(data, chromosome_width = 0), "positive")
  expect_error(ideogram_layout(data, chromosome_gap = -1), "non-negative")
  expect_error(ideogram_layout(data, max_chr_length = Inf), "finite")
  expect_error(ideogram_layout(data, row_gap = NA_real_), "non-negative")
  expect_error(ideogram_layout(data, unknown = 1), "Unused argument")
})

test_that("dimensionless layout has an informative print method", {
  layout <- ideogram_layout(as_ideogram_data(next_kar), ncol = 2)
  out <- capture.output(print(layout))
  expect_match(out[1], "<ideogram_layout_v2>", fixed = TRUE)
  expect_true(any(grepl("rows/columns: 2/2", out, fixed = TRUE)))
  expect_true(any(grepl("global (comparable)", out, fixed = TRUE)))
})
