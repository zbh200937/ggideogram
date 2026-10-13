responsive_fixture <- function() {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 70))
  ggideogram(k, orientation = "horizontal", show_names = FALSE,
    fill = "#4477AA", colour = "#4477AA", linewidth = 0) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
}

responsive_capture <- function(grob, width, height, slot = c(1, 1)) {
  capture <- ragg::agg_capture(width = width, height = height,
    units = "mm", res = 100, background = "white")
  on.exit(grDevices::dev.off())
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(width = slot[1], height = slot[2]))
  grid::grid.draw(grob)
  grid::popViewport()
  pixels <- toupper(capture()) == "#4477AA"
  index <- which(pixels, arr.ind = TRUE)
  c(width = diff(range(index[, 2])) / ncol(pixels),
    height = diff(range(index[, 1])) / nrow(pixels)) / slot
}

test_that("the same grob fills different devices and allocated child viewports", {
  skip_if_not_installed("ragg")
  p <- responsive_fixture()
  original <- p$coordinates$layout
  g <- ggplot2::ggplotGrob(p)
  for (size in list(c(150, 60), c(80, 130))) {
    occupied <- responsive_capture(g, size[1], size[2])
    expect_gt(occupied[["width"]], 0.94)
    expect_gt(occupied[["height"]], 0.80)
  }
  for (slot in list(c(0.9, 0.35), c(0.35, 0.9))) {
    occupied <- responsive_capture(g, 180, 140, slot)
    expect_gt(occupied[["width"]], 0.94)
    expect_gt(occupied[["height"]], 0.80)
  }
  expect_identical(p$coordinates$layout, original)
})

test_that("native polygons regenerate circular caps in the actual panel", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file, width = 8, height = 8)
  on.exit({ grDevices::dev.off(); unlink(file) })
  for (orientation in c("horizontal", "vertical")) {
    p <- ggideogram(data.frame(Chr = "A", Start = 0, End = 100),
      orientation = orientation, show_names = FALSE, curve_points = 24)
    built <- ggplot2::ggplot_build(p)
    panel <- built$layout$panel_params[[1]]
    coord <- built$layout$coord
    body <- p$layers[[1]]$geom$draw_panel(built$data[[1]], panel, coord, curve_points = 24)
    for (mm in list(c(160, 50), c(50, 160))) {
      grid::pushViewport(grid::viewport(width = grid::unit(mm[1], "mm"),
        height = grid::unit(mm[2], "mm")))
      drawn <- grid::makeContent(body)$children[[1]]
      xy <- cbind(as.numeric(drawn$x) * mm[1], as.numeric(drawn$y) * mm[2])
      tip <- xy[1, ]
      edge <- xy[24, ]
      centre <- if (orientation == "horizontal") c(edge[1], tip[2]) else c(tip[1], edge[2])
      distances <- sqrt(rowSums(sweep(xy[1:24, , drop = FALSE], 2, centre)^2))
      expect_equal(distances, rep(distances[1], 24), tolerance = 1e-8)
      grid::popViewport()
    }
  }
})

test_that("internal fill masks share the responsive chromosome boundary", {
  skip_if_not_installed("ragg")
  p <- responsive_fixture()
  d <- data.frame(Chr = c("A", "B"), Start = 1, End = c(100, 70))
  filled <- p + geom_track(track = "body", side = "overlay", width = 1.2) +
    geom_chr(component = "fill", track = "body", data = d,
      ggplot2::aes(chr = Chr, start = Start, end = End), fill = "#4477AA")
  for (size in list(c(150, 60), c(80, 130))) {
    plain <- responsive_capture(ggplot2::ggplotGrob(p), size[1], size[2])
    masked <- responsive_capture(ggplot2::ggplotGrob(filled), size[1], size[2])
    expect_equal(masked, plain, tolerance = 0.01)
  }
})
