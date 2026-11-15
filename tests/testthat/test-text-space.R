text_space_fixture <- function(orientation, label = NULL, axis = FALSE, reverse = FALSE) {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  d <- data.frame(Chr = rep(k$Chr, each = 3), Pos = rep(c(20, 40, 60), 2),
    Value = rep(0:2, 2))
  ggideogram(k, orientation = orientation, reverse_chr = if (reverse) k$Chr else character(),
    tracks = list(signal = geom_track(
    data = d, ggplot2::aes(chr = Chr, x = Pos, y = Value),
    geom = ggplot2::geom_line(), label = label, axis = axis)))
}

test_that("removing track text restores the same space as a fresh construction", {
  for (orientation in c("horizontal", "vertical", "circular")) {
    p <- text_space_fixture(orientation, "A much longer track title",
      list(position = "start", breaks = 0:2))
    margin <- p$theme$plot.margin
    q <- p + geom_track(track = "signal", label = NULL, axis = FALSE)
    fresh <- text_space_fixture(orientation)
    for (field in c("plot.margin", "legend.box.spacing")) {
      expect_equal(grid::convertUnit(q$theme[[field]], "mm", valueOnly = TRUE),
        grid::convertUnit(fresh$theme[[field]], "mm", valueOnly = TRUE))
    }
    for (field in c("plot.title", "plot.subtitle")) {
      expect_equal(grid::convertUnit(q$theme[[field]]$margin, "mm", valueOnly = TRUE),
        grid::convertUnit(fresh$theme[[field]]$margin, "mm", valueOnly = TRUE))
    }
    expect_identical(p$theme$plot.margin, margin)
    expect_identical(q$coordinates$layout$data, p$coordinates$layout$data)
  }
})

test_that("explicit theme spacing survives removal of automatic text space", {
  p <- text_space_fixture("horizontal", "A much longer track title") +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 13, 14, 15, unit = "mm"),
      plot.title = ggplot2::element_text(colour = "red",
        margin = ggplot2::margin(1, 2, 18, 4, unit = "mm")))
  q <- p + geom_track(track = "signal", label = NULL)
  expect_equal(grid::convertUnit(q$theme$plot.margin, "mm", valueOnly = TRUE), c(12, 13, 14, 15))
  expect_equal(grid::convertUnit(q$theme$plot.title$margin, "mm", valueOnly = TRUE), c(1, 2, 18, 4))
  expect_identical(q$theme$plot.title$colour, "red")
})

test_that("track titles clear the rendered tick labels at either source endpoint", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (orientation in c("horizontal", "vertical")) for (reverse in c(FALSE, TRUE)) {
    p <- text_space_fixture(orientation, "Signal",
      list(position = if (reverse) "end" else "start", breaks = 0:2), reverse = reverse)
    for (mm in list(c(120, 70), c(180, 120))) {
      grid::pushViewport(grid::viewport(width = grid::unit(mm[1], "mm"),
        height = grid::unit(mm[2], "mm")))
      b <- ggplot2::ggplot_build(p)
      draw <- function(index) {
        layer <- b$plot$layers[[index]]
        grob <- do.call(layer$geom$draw_panel, c(list(data = b$data[[index]],
          panel_params = b$layout$panel_params[[1]], coord = b$layout$coord),
          layer$computed_geom_params[intersect(names(layer$computed_geom_params), layer$geom$parameters())]))
        if (inherits(grob, "ideogram_track_axis_text")) grid::makeContent(grob)$children[[1]] else grob
      }
      title <- draw(which(vapply(p$layers, function(l)
        !is.null(l$ideogram_track_title), logical(1))))
      ticks <- draw(which(vapply(p$layers, function(l)
        inherits(l$stat, "StatTrackAxis") && l$stat_params$axis_part == 3L, logical(1))))
      dimension <- if (orientation == "horizontal") "width" else "height"
      extent <- vapply(seq_along(ticks$label), function(i) {
        g <- grid::textGrob(ticks$label[i], gp = ticks$gp)
        grid::convertUnit(if (dimension == "width") grid::grobWidth(g) else
          grid::grobHeight(g), "mm", valueOnly = TRUE)
      }, numeric(1))
      if (orientation == "horizontal") {
        title_near <- grid::convertX(title$x, "mm", valueOnly = TRUE)
        tick_far <- grid::convertX(ticks$x, "mm", valueOnly = TRUE) - extent
        expect_lt(title_near, min(tick_far))
      } else {
        title_near <- grid::convertY(title$y, "mm", valueOnly = TRUE)
        tick_far <- grid::convertY(ticks$y, "mm", valueOnly = TRUE) + extent
        expect_gt(title_near, max(tick_far))
      }
      grid::popViewport()
    }
  }
})

test_that("track titles share their start edge even with different axis and label settings", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (orientation in c("horizontal", "vertical")) {
    k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
    p <- ggideogram(k, orientation = orientation, ncol = 2, reverse_chr = "B",
      tracks = list(
        density = geom_track(chr = "A", limits = c(0, 1000), label = "Density",
          axis = list(position = "start", gap = 0.5, breaks = c(0, 1000), size = 2)),
        coverage = geom_track(chr = "B", limits = c(0, 1), axis = TRUE,
          label = list(text = "Coverage", size = 3.1, gap = 0.75))))
    q <- p + geom_track(track = "density", axis = list(gap = 1.2))
    for (plot in list(p, q)) {
      built <- ggplot2::ggplot_build(plot)
      titles <- which(vapply(plot$layers, function(layer)
        !is.null(layer$ideogram_track_title), logical(1)))
      positions <- vapply(titles, function(i) {
        layer <- built$plot$layers[[i]]
        params <- layer$computed_geom_params
        grob <- do.call(layer$geom$draw_panel, c(list(data = built$data[[i]],
          panel_params = built$layout$panel_params[[1]], coord = built$layout$coord),
          params[intersect(names(params), layer$geom$parameters())]))
        if (orientation == "horizontal") grid::convertX(grob$x, "mm", TRUE)
        else grid::convertY(grob$y, "mm", TRUE)
      }, numeric(1))
      expect_length(positions, 2)
      expect_equal(unname(positions[1]), unname(positions[2]))
    }
  }
})

test_that("adding an axis does not mutate the source plot's name spacing", {
  p <- text_space_fixture("horizontal", "Signal")
  name <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomChromosomeName"), logical(1)))
  before <- p$layers[[name]]$geom_params
  q <- p + geom_track_axis("signal", position = "start", breaks = 0:2)
  expect_identical(p$layers[[name]]$geom_params, before)
  expect_no_warning(ggplot2::ggplotGrob(q))
})

test_that("adding a selected bp axis preserves an unrelated track title's space", {
  p <- text_space_fixture("horizontal", "Signal", list(position = "start", breaks = 0:2))
  title_space <- function(plot) {
    title <- Filter(function(layer) !is.null(layer$ideogram_track_title), plot$layers)[[1]]
    title$geom_params$text_clearance
  }
  q <- p + geom_chr(chr = "B", component = "axis", side = "right", breaks = c(0, 80))
  expect_equal(title_space(q), title_space(p))
  expect_no_warning(ggplot2::ggplotGrob(q))
})

test_that("a centered title near the panel edge remains inside the exported image", {
  skip_if_not_installed("ragg")
  p <- ggideogram(data.frame(Chr = c("A", "B"), Start = 0, End = 100),
    show_names = FALSE, tracks = list(
      signal = geom_track(side = "left", width = 1,
        label = list(text = "Transcriptional abundance", colour = "#CC00CC")),
      annotations = geom_track(width = 25)))
  for (width in c(100, 160)) {
    capture <- ragg::agg_capture(width = width, height = 100,
      units = "mm", res = 100, background = "white")
    print(p)
    pixels <- capture()
    grDevices::dev.off()
    rgb <- grDevices::col2rgb(pixels)
    title <- matrix(rgb[1, ] > 150 & rgb[3, ] > 150 & rgb[2, ] < 100,
      nrow = nrow(pixels))
    columns <- which(title, arr.ind = TRUE)[, "col"]
    expect_gt(length(columns), 0)
    expect_gt(min(columns), 1)
    expect_lt(max(columns), ncol(pixels))
  }
})

test_that("constructing ideograms, tracks and spread/repel labels never opens a graphics device", {
  while (grDevices::dev.cur() > 1L) grDevices::dev.off()
  old_opt <- options(device = function(...) stop("Unexpected graphics device opened"))
  old_wd <- getwd()
  ro_dir <- tempfile("ideogram_readonly_")
  dir.create(ro_dir)
  Sys.chmod(ro_dir, "0555")
  on.exit({
    setwd(old_wd)
    Sys.chmod(ro_dir, "0755")
    unlink(ro_dir, recursive = TRUE)
    options(old_opt)
  }, add = TRUE)
  setwd(ro_dir)

  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100e6, 80e6))
  d <- data.frame(Chr = rep(k$Chr, each = 3), Pos = rep(c(20e6, 40e6, 60e6), 2),
    Value = rep(c(1, 5, 9), 2))
  loci <- data.frame(Chr = c("A", "A", "B"), Pos = c(15e6, 18e6, 30e6),
    Gene = c("LOC_A001", "LOC_A002", "LOC_B001"))

  expect_no_warning({
    p <- ggideogram(k, orientation = "horizontal", axis = TRUE, base_family = "Arial") +
      geom_track(track = "sig", data = d,
        ggplot2::aes(chr = Chr, x = Pos, y = Value),
        geom = ggplot2::geom_line(), label = "Signal", axis = TRUE) +
      geom_locus(data = loci,
        ggplot2::aes(chr = Chr, position = Pos, label = Gene),
        geom = "text", position = "spread") +
      ggplot2::labs(title = "Title", subtitle = "Subtitle")
    built <- ggplot2::ggplot_build(p)
  })
  expect_identical(unname(grDevices::dev.cur()), 1L)
  expect_length(list.files(ro_dir, all.files = TRUE, no.. = TRUE), 0L)
  expect_s3_class(built, "ggplot_built")

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    expect_no_warning({
      p_repel <- ggideogram(k, base_family = "Arial") +
        geom_locus(data = loci,
          ggplot2::aes(chr = Chr, position = Pos, label = Gene),
          geom = "text", position = "repel")
    })
    expect_identical(unname(grDevices::dev.cur()), 1L)
    expect_s3_class(p_repel, "ggplot")
  }

  # Preserve element_blank() on plot.title across track rebuilds
  p_blank <- ggideogram(k) +
    ggplot2::theme(plot.title = ggplot2::element_blank()) +
    geom_track(track = "sig", data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value),
      geom = ggplot2::geom_line(), label = "Signal")
  expect_true(inherits(p_blank$theme$plot.title, "element_blank"))
  expect_identical(unname(grDevices::dev.cur()), 1L)
})

test_that("axis titles and captions clear outward-extending bp ticks and tick labels", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  k <- data.frame(Chr = c("Chr1", "Chr2"), Start = 0, End = c(25e6, 20e6))
  p <- ggideogram(k, orientation = "horizontal", axis = TRUE, axis_side = "right",
    axis_tick_length = 1.5, axis_label_gap = 0.25, axis_size = 2.8) +
    ggplot2::labs(x = "Genomic position (Mb)", y = "Chromosomes", caption = "Source") +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10))

  tick_glyph <- grid::textGrob("20 Mb", gp = grid::gpar(fontsize = 2.8 * ggplot2::.pt))
  required_bottom <- 1.5 + 0.25 * 2.8 +
    grid::convertHeight(grid::grobHeight(tick_glyph), "mm", valueOnly = TRUE)

  xlab_margin <- grid::convertUnit(
    ggplot2::calc_element("axis.title.x.bottom", p$theme)$margin,
    "mm", valueOnly = TRUE
  )
  caption_margin <- grid::convertUnit(
    ggplot2::calc_element("plot.caption", p$theme)$margin,
    "mm", valueOnly = TRUE
  )
  ylab_margin <- grid::convertUnit(
    ggplot2::calc_element("axis.title.y.left", p$theme)$margin,
    "mm", valueOnly = TRUE
  )
  expect_gte(xlab_margin[1], required_bottom - 1e-6)
  expect_gte(caption_margin[1], required_bottom - 1e-6)
  expect_gt(ylab_margin[2], 3)

  gt <- ggplot2::ggplotGrob(p)
  xlab_row <- gt$layout$t[which(gt$layout$name == "xlab-b")]
  expect_gt(grid::convertHeight(gt$heights[xlab_row], "mm", valueOnly = TRUE), required_bottom + 2)
})

test_that("near-boundary linear axes, nameless circular axes and tick-only value axes keep clearance", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  # Row 2 has a slightly wider left marker than Row 1's left bp axis
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100e6, 80e6))
  m <- data.frame(Chr = "B", Pos = 40e6)
  p_near <- ggideogram(k, ncol = 1, axis = "A", axis_side = "left", axis_gap = 0.3,
    axis_size = 3.5) +
    geom_locus(data = m, ggplot2::aes(chr = Chr, position = Pos), side = "left", gap = 0.5)
  left_margin <- grid::convertUnit(p_near$theme$plot.margin, "mm", valueOnly = TRUE)[4]
  expect_gt(left_margin, 6)

  # Circular outer bp axis with show_names = FALSE and split geom_chr(component = 'axis')
  p_circ <- ggideogram(k, orientation = "circular", show_names = FALSE) +
    geom_chr(chr = "A", component = "axis", side = "outer", size = 3.5, gap = 0.3) +
    geom_chr(chr = "B", component = "axis", side = "outer", size = 3.5, gap = 0.6)
  circ_margins <- grid::convertUnit(p_circ$theme$plot.margin, "mm", valueOnly = TRUE)
  expect_equal(p_circ$coordinates$layout$circular_name_axis$gap, 0.6)
  expect_true(all(circ_margins >= 3.5))

  # Tick-only value axis at start still clears chromosome name
  d <- data.frame(Chr = "A", Pos = c(20e6, 80e6), Value = c(0, 10))
  p_ticks <- ggideogram(k[1, ], orientation = "vertical",
    tracks = list(sig = geom_track(data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value),
      geom = ggplot2::geom_line(),
      axis = list(position = "start", labels = NULL, tick_length = 2.5))))
  name_layer <- Filter(function(l) !is.null(l$ideogram_name_gap), p_ticks$layers)[[1]]
  expect_gt(name_layer$geom_params$name_end_gap[1], 0)
  expect_gte(
    grid::convertUnit(name_layer$geom_params$name_axis_clearance[1], "mm", valueOnly = TRUE),
    2.5 - 1e-6
  )
})

test_that("explicit small plot.margin survives track rebuilds and complete themes suppress layout grids", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  d <- data.frame(Chr = c("A", "B"), Pos = c(20, 40), Value = c(1, 2))
  p_early <- ggideogram(k, axis = TRUE) +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 2, 3, 4, unit = "mm")) +
    geom_track(track = "sig", data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_point())
  expect_equal(
    grid::convertUnit(p_early$theme$plot.margin, "mm", valueOnly = TRUE),
    c(1, 2, 3, 4)
  )

  p_min <- ggideogram(k) + ggplot2::theme_minimal()
  gt <- ggplot2::ggplotGrob(p_min)
  panel_grob <- gt$grobs[[which(gt$layout$name == "panel")]]
  bg_children <- vapply(
    panel_grob$children[[1]]$children,
    function(g) class(g)[1],
    character(1)
  )
  expect_false(any(grepl("polyline", bg_children)))
})
