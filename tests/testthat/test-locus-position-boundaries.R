test_that("native locus nudges survive projection and retain source anchors", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = 50, Start = 40, End = 60, Label = "gene")
  for (orientation in c("vertical", "horizontal", "circular")) {
    for (reverse in c(FALSE, TRUE)) {
      for (geom in c("point", "text", "interval")) {
        mapping <- switch(geom,
          point = ggplot2::aes(chr = Chr, position = Pos),
          text = ggplot2::aes(chr = Chr, position = Pos, label = Label),
          interval = ggplot2::aes(chr = Chr, start = Start, end = End))
        base <- ggideogram(k, show_names = FALSE, orientation = orientation,
          reverse_chr = if (reverse) "A" else character())
        identity <- base + geom_locus(mapping, d, geom = geom)
        nudged <- base + geom_locus(mapping, d, geom = geom,
          position = ggplot2::position_nudge(x = 0.4, y = 0.2))
        source <- ggplot2::ggplot_build(identity)$data[[length(identity$layers)]]
        actual_build <- ggplot2::ggplot_build(nudged)
        actual <- actual_build$data[[length(nudged$layers)]]
        layout <- actual_build$plot$coordinates$layout
        if (!all(source$ideogram_projected %in% TRUE) || is.null(source$ideogram_projected)) {
          source <- transform_ideogram_semantics(layout, source)
        }
        expect_equal(actual$x - source$x, 0.4)
        expect_equal(actual$y - source$y, 0.2)
        if (geom == "interval") {
          expect_equal(actual$xend - source$xend, 0.4)
          expect_equal(actual$yend - source$yend, 0.2)
          expect_equal(actual$start, d$Start)
          expect_equal(actual$end, d$End)
        } else {
          expect_equal(actual$position, d$Pos)
          expect_equal(actual$ideogram_anchor_position, d$Pos)
        }
        expect_s3_class(nudged$layers[[length(nudged$layers)]]$position, "PositionNudge")
        expect_no_error(ggplot2::ggplotGrob(nudged))
      }
    }
  }
})

test_that("a native Position evaluates displayed coordinates after a late track change", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(30, 70), Label = c("first", "second"))
  for (orientation in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(k, orientation = orientation, show_names = FALSE,
      tracks = list(annotations = geom_track(width = 2))) +
      geom_locus(ggplot2::aes(chr = Chr, position = Pos), d,
        track = "annotations", geom = ggplot2::geom_point(
          position = ggplot2::position_jitter(width = .2, height = .1, seed = 9)))
    q <- p + geom_track(track = "annotations", width = 4)
    actual <- function(plot) {
      b <- ggplot2::ggplot_build(plot)
      i <- which(vapply(plot$layers, function(x) inherits(x$stat, "StatChrMarker"), logical(1)))
      b$data[[i]]
    }
    first <- actual(p); second <- actual(q)
    expected <- offset_chr_points_signed(q$coordinates$layout, d$Chr,
      project_positions_checked(q$coordinates$layout, d$Chr, d$Pos), rep(2.7, nrow(d)))
    identity <- transform_ideogram_semantics(p$coordinates$layout, first)
    expect_gt(max(abs(first$x - identity$x) + abs(first$y - identity$y)), 0)
    expect_equal(second$x - expected$x, first$x - identity$x)
    expect_equal(second$y - expected$y, first$y - identity$y)
    expect_equal(second$ideogram_anchor_position, d$Pos)
    expect_no_error(ggplot2::ggplotGrob(q))
  }
})

test_that("native rectangle widths are cropped at displayed endpoints", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 80), Value = c(1, 1))
  for (orientation in c("vertical", "horizontal", "circular")) {
    for (reverse in c(FALSE, TRUE)) {
      p <- ggideogram(chr_view(k, "A", 20, 80), orientation = orientation,
        reverse_chr = if (reverse) "A" else character(), show_names = FALSE,
        tracks = list(counts = geom_track(data = d,
          ggplot2::aes(chr = Chr, x = Pos, y = Value), limits = c(0, 1),
          geom = ggplot2::geom_col(width = 30))))
      b <- ggplot2::ggplot_build(p)
      i <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomCol"), logical(1)))
      raw <- b$data[[i]]
      expect_equal(raw$xmin - raw$ideogram_position_shift, c(5, 65))
      expect_equal(raw$xmax - raw$ideogram_position_shift, c(35, 95))
      drawn <- transform_ideogram_track(b$plot$coordinates$layout, raw)
      g <- b$plot$coordinates$layout$chrom
      long_bounds <- if (orientation == "vertical") range(g$.axis_start_y, g$.axis_end_y) else
        range(g$.axis_start_x, g$.axis_end_x)
      extent <- if (orientation == "vertical") range(drawn$ymin, drawn$ymax) else
        range(drawn$xmin, drawn$xmax)
      expect_equal(extent, long_bounds)
      expect_equal(raw$position, d$Pos)
      expect_equal(raw$value, d$Value)
      expect_s3_class(p$layers[[i]]$geom, "GeomCol")
      expect_no_error(ggplot2::ggplotGrob(p))
    }
  }
})

test_that("overlay tile clip off keeps the complete native rectangle", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = 20, Value = .5)
  for (orientation in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(chr_view(k, "A", 20, 80), orientation = orientation,
      show_names = FALSE, tracks = list(tiles = geom_track(side = "overlay",
        data = d, ggplot2::aes(chr = Chr, x = Pos, y = Value),
        limits = c(0, 1), clip = "off", geom = ggplot2::geom_tile(width = 30, height = 1))))
    b <- ggplot2::ggplot_build(p)
    i <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomTile"), logical(1)))
    raw <- b$data[[i]]
    drawn <- transform_ideogram_track(b$plot$coordinates$layout, raw)
    at <- project_positions_raw(p$coordinates$layout$chrom, "A", c(5, 35))
    expect_equal(if (orientation == "vertical") range(drawn$ymin, drawn$ymax) else
      range(drawn$xmin, drawn$xmax), if (orientation == "vertical") range(at$y) else range(at$x))
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that("overlay columns crop sector boundaries unless clipping is disabled", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = 20, Value = 1)
  for (clip in c("auto", "off")) {
    p <- ggideogram(chr_view(k, "A", 20, 80), orientation = "circular",
      show_names = FALSE, tracks = list(counts = geom_track(side = "overlay",
        data = d, ggplot2::aes(chr = Chr, x = Pos, y = Value),
        limits = c(0, 1), clip = clip, geom = ggplot2::geom_col(width = 30))))
    b <- ggplot2::ggplot_build(p)
    i <- which(vapply(p$layers, function(x) inherits(x$geom, "GeomCol"), logical(1)))
    drawn <- transform_ideogram_track(b$plot$coordinates$layout, b$data[[i]])
    expect_equal(drawn$xmin, if (clip == "auto") 0 else
      project_positions_raw(p$coordinates$layout$chrom, "A", 5)$x)
    expect_equal(b$data[[i]]$position, d$Pos)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
  p <- ggideogram(chr_view(k, "A", 20, 80), orientation = "circular", show_names = FALSE,
    tracks = list(tiles = geom_track(side = "overlay", limits = c(0, 1)))) +
    geom_track_tile(data = d, ggplot2::aes(chr = Chr, position = Pos, value = Value),
      track = "tiles", clip = "off", width = 30, height = 0)
  b <- ggplot2::ggplot_build(p)
  raw <- b$data[[length(p$layers)]]
  drawn <- transform_ideogram_track(b$plot$coordinates$layout, raw)
  expect_equal(drawn$xmin, project_positions_raw(p$coordinates$layout$chrom, "A", 5)$x)
})

test_that("nudged locus leaders keep their scientific start anchor", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = 50)
  for (orientation in c("vertical", "horizontal", "circular")) {
    base <- ggideogram(k, orientation = orientation, show_names = FALSE)
    identity <- base + geom_locus(ggplot2::aes(chr = Chr, position = Pos), d, geom = "link")
    nudged <- base + geom_locus(ggplot2::aes(chr = Chr, position = Pos), d, geom = "link",
      position = ggplot2::position_nudge(x = .4, y = .2))
    original <- ggplot2::ggplot_build(identity)$data[[length(identity$layers)]]
    if (is.null(original$ideogram_projected)) {
      original <- transform_ideogram_semantics(identity$coordinates$layout, original)
    }
    actual <- ggplot2::ggplot_build(nudged)$data[[length(nudged$layers)]]
    expect_equal(actual$x, original$x)
    expect_equal(actual$y, original$y)
    expect_equal(actual$xend - original$xend, .4)
    expect_equal(actual$yend - original$yend, .2)
    expect_equal(actual$ideogram_anchor_position, d$Pos)
    expect_no_error(ggplot2::ggplotGrob(nudged))
  }
})
