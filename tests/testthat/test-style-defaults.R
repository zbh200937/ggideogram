test_that("explicit bp labels retain precision and automatic ticks use round coordinates", {
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
    c(5000, 10000, 15000))
})

test_that("complete themes do not expose internal layout axes", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  for (orientation in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(k, orientation = orientation, axis = TRUE) +
      ggplot2::theme_classic() + geom_track(track = "signal", width = 2)
    gt <- ggplot2::ggplotGrob(p)
    slots <- grepl("^(axis|xlab|ylab)-", gt$layout$name)
    expect_true(all(vapply(gt$grobs[slots], inherits, logical(1), "zeroGrob")))
    expect_true(any(vapply(p$layers, function(layer)
      !is.null(layer$ideogram_axis_spec), logical(1))))
  }
})

test_that("local automatic axes stay round and explicit breaks stay exact", {
  k <- data.frame(Chr = "A", Start = 0, End = 1e8)
  view <- chr_view(k, "A", 79007836, 81007836)
  axis_text <- function(p) Filter(function(layer)
    inherits(layer$geom, "GeomIdeogramAxisText"), p$layers)[[1]]
  for (orientation in c("vertical", "horizontal")) {
    p <- ggideogram(view, orientation = orientation, axis = TRUE, axis_units = "Mb")
    labels <- axis_text(p)
    expect_equal(labels$data$label, c("79.5 Mb", "80.0 Mb", "80.5 Mb", "81.0 Mb"))
    expect_true(labels$geom_params$check_overlap)
    explicit <- c(79007836, 80007836)
    q <- ggideogram(view, orientation = orientation, axis = TRUE,
      axis_units = "Mb", axis_breaks = explicit, reverse_chr = "A")
    expect_equal(axis_text(q)$data$label, c("79.007836 Mb", "80.007836 Mb"))
    expect_false(axis_text(q)$geom_params$check_overlap)
    expect_equal(p$coordinates$layout$data$karyotype, as_ideogram_data(k)$karyotype)
  }
})

test_that("endpoint names align with bodies by default and tracks by request", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  name_data <- function(p) {
    b <- ggplot2::ggplot_build(p)
    i <- which(vapply(b$plot$layers, function(l) inherits(l$stat, "StatChromosomeName"), logical(1)))
    b$data[[i]]
  }
  for (orientation in c("vertical", "horizontal")) {
    p <- ggideogram(k, orientation = orientation,
      tracks = list(labels = geom_track(side = "right", width = 16, gap = 0)))
    g <- p$coordinates$layout$chrom
    axis <- if (orientation == "vertical") "x" else "y"
    expect_equal(name_data(p)[[axis]], g[[paste0(".axis_start_", axis)]])
    q <- ggideogram(k, orientation = orientation, name_align = "tracks",
      tracks = list(labels = geom_track(side = "right", width = 16, gap = 0)))
    expect_equal(abs(name_data(q)[[axis]] - name_data(p)[[axis]]), rep(8.1, 2))
  }
})

test_that("vertical bp axes reserve room in narrow chromosome gaps", {
  k <- data.frame(Chr = c("A", "B", "C"), Start = 0, End = c(100, 90, 80))
  base <- ggideogram(k)
  p <- ggideogram(k, axis = TRUE)
  expect_equal(ideogram_plot_layout(base)$chromosome_gap, 3)
  expect_equal(ideogram_plot_layout(p)$chromosome_gap, 8)
  explicit <- ggideogram(k, axis = TRUE, chromosome_gap = 4)
  expect_equal(ideogram_plot_layout(explicit)$chromosome_gap, 4.3)
  tight <- ggideogram(k, axis = TRUE, chromosome_gap = 4, min_gap = 0)
  expect_equal(ideogram_plot_layout(tight)$chromosome_gap, 4.3)
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
      expect_true(all(vapply(axes, function(x)
        identical(x$aes_params$size, ideogram_text_size("bp") / 3), logical(1))))
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

test_that("default bodies are wider and automatic layout uses annotation space", {
  k <- data.frame(Chr = paste0("Chr", 1:15), Start = 0, End = 100)
  tracks <- list(labels = geom_track(width = 16))
  p <- ggideogram(k, tracks = tracks)
  l <- p$coordinates$layout
  expect_equal(l$chromosome_width, 1.2)
  expect_gt(l$nrow, 1)
  expect_lt(abs(log(diff(l$bounds$x) / diff(l$bounds$y))),
    abs(log((15 * 17.2 + 14 * 3) / 40)))
  explicit <- ggideogram(k, tracks = tracks, ncol = 15, chromosome_width = 2)
  expect_equal(explicit$coordinates$layout$ncol, 15)
  expect_equal(explicit$coordinates$layout$chromosome_width, 2)
  horizontal <- ggideogram(k[1:2, ], orientation = "horizontal", tracks = tracks,
    layout_aspect = 3)
  h <- horizontal$coordinates$layout
  expect_equal(horizontal$coordinates$aspect(list()), 1 / 3)
  expect_equal(h$max_chr_length, 40)
  fixed <- ggideogram(k[1:2, ], orientation = "horizontal", tracks = tracks,
    layout_aspect = 3, max_chr_length = 40)
  expect_equal(fixed$coordinates$layout$max_chr_length, 40)
})
