local_map <- function(chr = NULL, link_layout = "auto") {
  ggideogram(data.frame(Chr = c("A", "B", "C"), Start = 0, End = 100),
    orientation = "horizontal", ncol = 3, axis = TRUE, chr = chr,
    link_layout = link_layout)
}
local_edges <- function(b = c("B", "C")) {
  geom_chrlink(ggplot2::aes(chr1 = a, chr2 = b, position1 = p, position2 = q),
    data.frame(a = "A", b = b, p = 30, q = 70))
}

test_that("three-chromosome comparisons put the shared chromosome in the middle", {
  p <- local_map() + local_edges()
  expect_identical(p$coordinates$layout$chrom$.chr, c("B", "A", "C"))
  expect_identical((local_map(link_layout = "input") + local_edges())$coordinates$layout$chrom$.chr,
    c("A", "B", "C"))
  expect_identical((local_map(chr = c("C", "B", "A")) + local_edges())$coordinates$layout$chrom$.chr,
    c("C", "B", "A"))
  pairs <- Filter(function(x) inherits(x$ideogram_recipe, "ggideogram_pair_component"), p$layers)
  expect_s3_class(pairs[[1]]$geom, "GeomSegment")
  expect_equal(nrow(pairs[[1]]$data), 2)
  expect_equal(pairs[[1]]$data$.pair_start1, c(30, 30))
  expect_equal(microsynteny_side(p$coordinates$layout, c("B", "A", "C"), "axis"),
    c("left", "left", "right"))
})

test_that("annotations and later links share compact rebuilt source geometry", {
  loci <- data.frame(Chr = c("A", "B", "C"), Position = 30, Name = c("a", "b", "c"))
  labels <- geom_locus(ggplot2::aes(chr = Chr, position = Position, label = Name), loci,
    geom = "text", position = "spread")
  p <- local_map() + local_edges("B") + labels + local_edges("C")
  layout <- p$coordinates$layout
  expect_identical(layout$chrom$.chr, c("B", "A", "C"))
  expect_equal(annotation_extent(layout, c("B", "C"), c("right", "left")), c(0, 0))
  expect_equal(microsynteny_side(layout, c("B", "A", "C")), c("left", "right", "right"))
  fresh <- local_map() + local_edges("B") + local_edges("C") + labels
  expect_equal(layout$chrom, fresh$coordinates$layout$chrom)
  expect_equal(length(p$layers), length(fresh$layers))
  pairs <- Filter(function(x) inherits(x$ideogram_recipe, "ggideogram_pair_component"), p$layers)
  for (layer in pairs) {
    d <- layer$data
    source <- project_positions_checked(layout, d$.pair_chr1, d$.pair_start1)
    target <- project_positions_checked(layout, d$.pair_chr2, d$.pair_start2)
    expect_equal(d$.pair_x, source$x)
    expect_equal(d$.pair_xend, target$x)
  }
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("shared axes and middle icons leave both connection corridors clear", {
  loci <- data.frame(Chr = c("A", "B", "C"), Position = 30, Name = c("hub", "b", "c"))
  p <- local_map() + local_edges() +
    geom_locus(ggplot2::aes(chr = Chr, position = Position, label = Name), loci,
      geom = "text", position = "spread")
  axes <- Filter(function(x) !is.null(x$ideogram_axis_spec), p$layers)
  expect_equal(nrow(axes[[1]]$data), 1)
  expect_true("0 bp" %in% axes[[3]]$data$label)
  labels <- Filter(function(x) inherits(x$stat, "StatChrLabels"), p$layers)
  expect_false(any(vapply(labels, function(x) "A" %in% x$data$.label_chr, logical(1))))
  icon <- Filter(function(x) inherits(x$geom, "GeomPoint") && identical(x$aes_params$shape, 23), p$layers)
  expect_length(icon, 1)
  expect_equal(icon[[1]]$data$Chr, "A")
  pairs <- Filter(function(x) inherits(x$ideogram_recipe, "ggideogram_pair_component"), p$layers)[[1]]$data
  center <- p$coordinates$layout$chrom$.axis_start_y[2]
  expect_true(any(pairs$.pair_y > center))
  expect_true(any(pairs$.pair_y < center))
  expect_equal(abs(pairs$.pair_y - center), rep(p$coordinates$layout$chromosome_width / 2, 2))
  native <- ggideogram(data.frame(Chr = c("A", "B", "C"), Start = 0, End = c(100, 200, 400)),
    orientation = "horizontal", ncol = 3, axis = TRUE) + local_edges()
  layout <- native$coordinates$layout
  span <- layout$chrom$.axis_end_x - layout$chrom$.axis_start_x
  expect_equal(span / min(span), c(2, 1, 4))
  expect_equal(layout$chromosome_gap, 9)
  expect_equal(abs(diff(layout$chrom$.axis_start_y)) - layout$chromosome_width, c(9, 9))
  axis_text <- Filter(function(x) !is.null(x$ideogram_axis_spec) && "label" %in% names(x$data), native$layers)
  expect_true("400 bp" %in% axis_text[[1]]$data$label)
})


test_that("two chromosomes retain their own source-coordinate axes", {
  p <- ggideogram(data.frame(Chr = c("A", "B"), Start = c(1000, 2000), End = c(1400, 2200)),
    orientation = "horizontal", axis = TRUE, axis_units = "bp") +
    geom_chrlink(ggplot2::aes(chr1 = a, chr2 = b, position1 = p, position2 = q),
      data.frame(a = "A", b = "B", p = 1200, q = 2100))
  axes <- Filter(function(x) !is.null(x$ideogram_axis_spec), p$layers)
  expect_equal(nrow(axes[[1]]$data), 2)
  labels <- axes[[3]]$data
  expect_true(all(c("1000 bp", "1400 bp", "2000 bp", "2200 bp") %in% labels$label))
  expect_true(any(labels$ny > 0) && any(labels$ny < 0))
  expect_equal(p$coordinates$layout$chromosome_gap, 9)
})

test_that("automatic annotations filter complete source tables and retain mapped icon colours", {
  k <- data.frame(Chr = LETTERS[1:4], Start = 0, End = 100)
  d <- data.frame(Chr = LETTERS[1:4], Pos = 40, Label = letters[1:4], Group = letters[1:4])
  p <- ggideogram(k, chr = c("A", "B", "C"), orientation = "horizontal") + local_edges()
  point <- p + geom_locus(ggplot2::aes(chr = Chr, position = Pos), d)
  expect_no_error(ggplot2::ggplotGrob(point))
  q <- p + geom_locus(ggplot2::aes(chr = Chr, position = Pos, label = Label, colour = Group),
    d, geom = "text", position = "spread") +
    ggplot2::scale_colour_manual(values = c(a = "red", b = "blue", c = "green"))
  b <- ggplot2::ggplot_build(q)
  i <- which(vapply(q$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
  expect_length(i, 1)
  expect_equal(b$data[[i]]$chr, "B")
  expect_equal(b$data[[i]]$colour, "blue")
  expect_no_warning(ggplot2::ggplotGrob(q))
})
test_that("chromosome-scoped tracks share lanes and curved ribbons stay in the corridor", {
  k <- data.frame(Chr = c("A", "B", "C"), Start = 0, End = 100)
  d <- data.frame(Chr = rep(k$Chr, each = 2), Pos = c(10, 90), Value = c(0, 1))
  edges <- data.frame(a = c("A", "B"), b = c("B", "C"), start = 5, end = 15,
    start2 = 85, end2 = 95, strand = c("+", "-"))
  signal <- function(keys, side, width) geom_track(chr = keys, side = side, width = width,
    data = d, ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_line())
  p <- ggideogram(k, orientation = "horizontal", max_chr_length = 100, chromosome_gap = 2,
    tracks = list(outer_a = signal("A", "left", 5), inner_a = signal("A", "right", 1),
      inner_b = signal("B", "left", 2), inner_c = signal("C", "left", 3),
      shared = signal(c("B", "C"), "left", 1))) +
    geom_chrlink(ggplot2::aes(chr1 = a, chr2 = b, start1 = start, end1 = end,
      start2 = start2, end2 = end2, orientation = strand), edges, curvature = 1)
  for (plot in list(p, p + geom_track(track = "inner_b", width = 5))) {
    layout <- ideogram_plot_layout(plot)
    spec <- function(id) track_table_row(layout, id)
    expect_equal(spec("outer_a")$low_offset, spec("inner_b")$low_offset)
    expect_equal(spec("inner_b")$low_offset, spec("inner_c")$low_offset)
    expect_gt(abs(spec("shared")$low_offset), max(abs(spec("inner_b")$high_offset),
      abs(spec("inner_c")$high_offset)))
    layer <- Filter(function(layer) inherits(layer$geom, "GeomPolygon"), plot$layers)
    ribbon <- utils::tail(layer, 1)[[1]]$data
    for (i in seq_len(nrow(edges))) {
      rows <- ribbon[ribbon$.pair_id == i, ]
      y <- layout$chrom$.axis_start_y[match(c(edges$a[i], edges$b[i]), layout$chrom$.chr)]
      near <- y[1] - layout$chromosome_width / 2 - annotation_extent(layout, edges$a[i], "right")
      far <- y[2] + layout$chromosome_width / 2 + annotation_extent(layout, edges$b[i], "left")
      expect_equal(range(rows$.pair_y), c(far, near))
      expect_equal(unique(rows$.pair_start1), edges$start[i])
      expect_equal(unique(rows$.pair_orientation), edges$strand[i])
    }
    expect_no_warning(ggplot2::ggplotGrob(plot))
  }
})
