test_that("closed locus intervals retain one base and agree with continuous projection", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- as_chr_features(data.frame(Chr = "A", Start = c(1, 11, 100),
    End = c(1, 20, 100)))
  original <- d
  for (orientation in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(k, orientation = orientation, reverse_chr = "A",
      show_names = FALSE) + geom_locus(geom = "interval", data = d,
        ggplot2::aes(chr = Chr, start = Start, end = End), linewidth = .3)
    b <- ggplot2::ggplot_build(p)
    rows <- utils::tail(b$data, 1)[[1]]
    actual <- p$coordinates$transform(rows, b$layout$panel_params[[1]])
    expected <- project_chr_interval(p, transform(d, Start = Start - 1),
      "Chr", "Start", "End")
    lo <- p$coordinates$transform(data.frame(x = expected$.x_start, y = expected$.y_start),
      b$layout$panel_params[[1]])
    hi <- p$coordinates$transform(data.frame(x = expected$.x_end, y = expected$.y_end),
      b$layout$panel_params[[1]])
    expect_equal(rows$start, d$Start)
    expect_equal(rows$end, d$End)
    expect_equal(actual$x, lo$x)
    expect_equal(actual$y, lo$y)
    expect_equal(actual$xend, hi$x)
    expect_equal(actual$yend, hi$y)
    expect_true(all((actual$x - actual$xend)^2 + (actual$y - actual$yend)^2 > 0))
    expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, "GeomSegment")
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(d, original)
})

test_that("closed interval views retain whole source bases and clip precise geometry", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  d <- data.frame(Chr = c(rep("A", 5), "B"), Start = c(1, 1, 11, 20, 21, 1),
    End = c(10, 15, 11, 20, 30, 100), ID = letters[1:6])
  original <- d
  for (window in list(c(10, 20), c(10.5, 19.5))) {
    view <- chr_view(k, "A", window[1], window[2])
    clipped <- view_chr_data(d, view)
    expect_identical(clipped$ID, c("b", "c", "d"))
    expect_equal(clipped$Start, c(11, 11, 20))
    expect_equal(clipped$End, c(15, 11, 20))
    expect_equal(clipped$.source_start, c(1, 11, 20))
    expect_equal(clipped$.source_end, c(15, 11, 20))
    for (orientation in c("horizontal", "circular")) {
      p <- ggideogram(view, orientation = orientation, show_names = FALSE,
        tracks = list(intervals = geom_track(side = "right", chr = "A", data = d,
          ggplot2::aes(chr = Chr, start = Start, end = End),
          layers = list(geom_locus(geom = "interval"))))) +
        geom_locus(geom = "interval", data = clipped,
          ggplot2::aes(chr = Chr, start = Start, end = End))
      layers <- Filter(function(l) inherits(l$stat, "StatChrInterval"), p$layers)
      expect_identical(layers[[1]]$data$ID, clipped$ID)
      expect_equal(layers[[1]]$data$.source_start, clipped$.source_start)
      b <- ggplot2::ggplot_build(p)
      rows <- Filter(function(x) "ideogram_interval" %in% names(x), b$data)
      expect_length(rows, 2)
      for (x in rows) {
        expect_equal(x$ideogram_position, pmax(c(10, 10, 19), window[1]))
        expect_equal(x$ideogram_end_position, pmin(c(15, 11, 20), window[2]))
      }
      expect_no_warning(ggplot2::ggplotGrob(p))
    }
    touching <- ggideogram(view, show_names = FALSE) +
      geom_locus(geom = "interval", data = d[c(1, 5), ],
        ggplot2::aes(chr = Chr, start = Start, end = End))
    expect_no_warning(b <- ggplot2::ggplot_build(touching))
    expect_equal(nrow(utils::tail(b$data, 1)[[1]]), 0)
    expect_no_warning(ggplot2::ggplotGrob(touching))
  }
  expect_identical(d, original)
  expect_equal(nrow(view_chr_data(d[3, ], chr_view(k, "A", 10.5, 10.6))), 1)
  invalid <- transform(d[1, ], Start = 0)
  expect_error(view_chr_data(invalid, chr_view(k, "A", 0, 20)), "closed integer")
  expect_warning(ggplot2::ggplot_build(ggideogram(k) +
    geom_locus(geom = "interval", data = invalid,
      ggplot2::aes(chr = Chr, start = Start, end = End))), "positive integer")
})

test_that("gene and transcript identity does not depend on character separators", {
  k <- data.frame(Chr = c("A.B", "A"), Start = 0, End = 100)
  d <- data.frame(Chr = rep(k$Chr, each = 2), Start = c(11, 16, 61, 66),
    End = c(30, 20, 80, 70), Model = rep(c("C", "B.C"), each = 2),
    Type = rep(c("exon", "CDS"), 2), Strand = rep(c("+", "-"), each = 2))
  original <- d
  for (mode in c("gene", "transcript")) {
    mapping <- ggplot2::aes(chr = Chr, start = Start, end = End,
      type = Type, strand = Strand)
    mapping[[mode]] <- ggplot2::aes(model = Model)$model
    p <- ggideogram(k, orientation = "horizontal", show_names = FALSE,
      tracks = list(models = geom_track(width = 3))) +
      geom_genemodel(mode = mode, data = d, mapping = mapping,
        track = "models", labels = FALSE)
    backbone <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)[[1]]$data
    expect_equal(nrow(backbone), 2)
    at <- match(k$Chr, backbone$.model_chr)
    expect_equal(backbone$.model_id[at], c("C", "B.C"))
    expect_equal(backbone$.model_start[at], c(10, 60))
    expect_equal(backbone$.model_end[at], c(30, 80))
    expect_equal(backbone$.model_strand[at], c("+", "-"))
    blocks <- Filter(function(l) inherits(l$geom, "GeomRect"), p$layers)
    expect_equal(sum(vapply(blocks, function(l) nrow(l$data), integer(1))), 4)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(d, original)
})
