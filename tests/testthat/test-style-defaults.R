test_that("bp labels preserve a local window's fractional unit offset", {
  expect_equal(axis_labels(c(2500, 4500, 14500), 1e3, "kb"),
    c("2.5 kb", "4.5 kb", "14.5 kb"))
  expect_equal(axis_labels(c(1000001, 2000001), 1e6, "Mb"),
    c("1.000001 Mb", "2.000001 Mb"))
  expect_equal(axis_labels(c(0, 50000000), 1e6, "Mb"), c("0 Mb", "50 Mb"))
  k <- data.frame(Chr = "A", Start = 0, End = 20000)
  p <- ggideogram(chr_view(k, "A", 2500, 15000), orientation = "horizontal", axis = TRUE)
  text <- Filter(function(x) !is.null(x$ideogram_axis_spec) &&
    inherits(x$geom, "GeomText"), p$layers)[[1]]$data
  expect_equal(as.numeric(sub(" kb", "", text$label)) * 1000,
    seq(2500, 14500, by = 2000))
})

test_that("vertical bp axes get default space and retain explicit spacing", {
  k <- data.frame(Chr = c("A", "B", "C"), Start = 0, End = c(100, 90, 80))
  base <- ggideogram(k)
  p <- ggideogram(k, axis = TRUE)
  expect_equal(ideogram_plot_layout(base)$chromosome_gap, 3)
  expect_equal(ideogram_plot_layout(p)$chromosome_gap, 8)
  explicit <- ggideogram(k, axis = TRUE, chromosome_gap = 4)
  expect_equal(ideogram_plot_layout(explicit)$chromosome_gap, 4)
  updated <- p + geom_track(track = "signal")
  expect_equal(ideogram_plot_layout(updated)$chromosome_gap, 8)
})

test_that("gene lanes reserve text margins and chromosome-name clearance", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Start = c(10, 40, 70), End = c(25, 55, 90),
    Gene = c("GENE_LONG_001", "GENE_LONG_002", "GENE_LONG_003"), Type = "exon", Strand = "+")
  track <- geom_track(track = "genes", data = d,
    ggplot2::aes(chr = Chr, start = Start, end = End, gene = Gene, type = Type, strand = Strand),
    layers = geom_genemodel())
  expect_equal(track$width, 6)
  expect_equal(geom_track(width = 4, layers = geom_genemodel())$width, 4)
  for (orientation in c("horizontal", "vertical")) {
    p <- ggideogram(k, orientation = orientation, base_family = "serif") + track
    label <- Filter(function(x) isTRUE(x$ideogram_gene_label), p$layers)[[1]]
    name <- Filter(function(x) !is.null(x$ideogram_name_gap), p$layers)[[1]]
    expect_s3_class(label$geom, "GeomText")
    expect_identical(label$aes_params$family, "serif")
    glyph <- grid::textGrob(d$Gene[1], gp = grid::gpar(fontfamily = "serif",
      fontsize = label$aes_params$size * ggplot2::.pt))
    extent <- if (orientation == "horizontal") {
      grid::grobWidth(glyph) * label$aes_params$hjust
    } else grid::grobHeight(glyph) * (1 - label$aes_params$vjust)
    required <- grid::convertUnit(extent, "mm", valueOnly = TRUE)
    at <- if (orientation == "horizontal") 4 else 1
    expect_gte(grid::convertUnit(p$theme$plot.margin[at], "mm", valueOnly = TRUE), required)
    expect_gte(grid::convertUnit(name$geom_params$name_axis_clearance[1], "mm", valueOnly = TRUE), required)
  }
})

test_that("track text inherits font family and rebuilt axes keep their margins", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 80), Value = c(0, 1))
  for (orientation in c("horizontal", "vertical", "circular")) {
    p <- ggideogram(k, orientation = orientation, base_family = "serif") +
      geom_track(track = "signal", data = d,
        ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_line(),
        label = "Signal", axis = TRUE)
    for (plot in list(p, p + geom_track(track = "signal", width = 4))) {
      text <- Filter(function(x) inherits(x$geom, "GeomText") &&
        !is.null(x$ideogram_scope), plot$layers)
      expect_true(all(vapply(text, function(x) identical(x$aes_params$family, "serif"), logical(1))))
      axes <- Filter(function(x) grepl("^axis:", x$ideogram_scope_slot %||% ""), text)
      expect_gt(length(axes), 0)
      expect_true(all(vapply(axes, function(x) identical(x$aes_params$size, 2.2), logical(1))))
      expect_gte(grid::convertUnit(plot$theme$legend.box.spacing, "mm", valueOnly = TRUE), 6)
      if (orientation == "horizontal") {
        axis <- Filter(function(x) inherits(x$stat, "StatTrackAxis"), text)[[1]]
        glyph <- grid::textGrob(axis$data$label[1], gp = grid::gpar(fontfamily = "serif",
          fontsize = axis$aes_params$size * ggplot2::.pt))
        required <- grid::grobWidth(glyph) + grid::unit(axis$geom_params$tick_length +
          axis$geom_params$label_gap * axis$aes_params$size, "mm")
        expect_gte(grid::convertUnit(plot$theme$plot.margin[2], "mm", valueOnly = TRUE),
          grid::convertUnit(required, "mm", valueOnly = TRUE))
      }
    }
    override <- p + geom_track(track = "signal", label = list(family = "mono"),
      axis = list(family = "mono", size = 3.1)) +
      ggplot2::theme(legend.box.spacing = grid::unit(9, "mm"))
    text <- Filter(function(x) inherits(x$geom, "GeomText") &&
      !is.null(x$ideogram_scope), override$layers)
    expect_true(all(vapply(text, function(x) identical(x$aes_params$family, "mono"), logical(1))))
    axes <- Filter(function(x) grepl("^axis:", x$ideogram_scope_slot %||% ""), text)
    expect_true(all(vapply(axes, function(x) identical(x$aes_params$size, 3.1), logical(1))))
    expect_equal(grid::convertUnit(override$theme$legend.box.spacing, "mm", valueOnly = TRUE), 9)
  }
})

test_that("gene label style is independent of blocks and survives relayout", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Start = 10, End = 90, Gene = "NAC001",
    Type = "exon", Strand = "+")
  for (orientation in c("horizontal", "vertical", "circular")) {
    p <- ggideogram(k, orientation = orientation) +
      geom_track(track = "genes", data = d,
        ggplot2::aes(chr = Chr, start = Start, end = End, gene = Gene,
          type = Type, strand = Strand),
        layers = geom_genemodel(fill = "steelblue", label_fontface = "italic",
          label_family = "serif", label_colour = "navy"))
    q <- p + geom_track(track = "genes", width = 8)
    for (plot in list(p, q)) {
      label <- Filter(function(x) isTRUE(x$ideogram_gene_label), plot$layers)[[1]]
      expect_identical(label$aes_params$fontface, "italic")
      expect_identical(label$aes_params$family, "serif")
      expect_identical(label$aes_params$colour, "navy")
      expect_no_warning(ggplot2::ggplotGrob(plot))
    }
  }
})

test_that("circular text reservation keeps explicit legend gaps during component updates", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 80), Value = c(0, 1), Group = c("a", "b"))
  automatic <- ggideogram(k, orientation = "circular", axis = TRUE,
    axis_size = 5, name_size = 4)
  expect_gt(grid::convertUnit(automatic$theme$legend.box.spacing, "mm", valueOnly = TRUE), 6)
  themes <- list(ggplot2::theme(legend.box.spacing = grid::unit(1, "mm")),
    theme_ideogram(legend.box.spacing = grid::unit(1, "mm")))
  for (settings in themes) {
    p <- automatic + settings
    q <- p + geom_track(track = "signal", data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Group),
      geom = ggplot2::geom_point(), axis = TRUE, label = "Signal")
    r <- q + geom_chr_axis(chr = "A", size = 4)
    expect_equal(grid::convertUnit(q$theme$legend.box.spacing, "mm", valueOnly = TRUE), 1)
    expect_equal(grid::convertUnit(r$theme$legend.box.spacing, "mm", valueOnly = TRUE), 1)
    expect_no_warning(ggplot2::ggplotGrob(r))
  }
})
