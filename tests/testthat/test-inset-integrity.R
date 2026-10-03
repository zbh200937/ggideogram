inset_integrity_fixture <- function(component = NULL) {
  karyotype <- data.frame(Chr = "A", Start = 0, End = 100)
  if (is.null(component)) {
    component <- ggplot2::ggplot(
      data.frame(x = 1:3, y = c(2, 4, 1),
        Class = c("child-A", "child-B", "child-A")),
      ggplot2::aes(x, y, colour = Class)
    ) + ggplot2::geom_point(size = 1.2) +
      ggplot2::labs(colour = "Child group") +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        legend.background = ggplot2::element_rect(fill = "white", colour = NA))
  }
  data <- data.frame(Chr = "A", Pos = 50)
  data$Plot <- list(component)
  base <- ggideogram(karyotype, show_names = FALSE, max_chr_length = 40,
    fill = "#4477AA", colour = "#4477AA",
    tracks = track_layout(inset = track("right", width = 8)))
  inset <- function(clip = "on", hjust = NULL) geom_chr_inset(
    data = data, ggplot2::aes(chr = Chr, position = Pos, plot = Plot),
    track = "inset", width = grid::unit(0.55, "npc"),
    height = grid::unit(0.32, "npc"), clip = clip, hjust = hjust)
  list(base = base, child = component, inset = inset)
}

inset_integrity_body_pixels <- function(plot) {
  capture <- ragg::agg_capture(
    width = 100, height = 120, units = "mm", res = 150, background = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  sum(toupper(capture()) == "#4477AA")
}

test_that("complete child guides cannot overflow a clipped inset into the host", {
  skip_if_not_installed("ragg")
  fixture <- inset_integrity_fixture()
  original <- fixture$child
  baseline <- inset_integrity_body_pixels(fixture$base)
  clipped <- fixture$base + fixture$inset("on")
  overflow <- fixture$base + fixture$inset("off")
  expect_equal(inset_integrity_body_pixels(clipped), baseline)
  expect_lt(inset_integrity_body_pixels(overflow), baseline)
  expect_identical(fixture$child, original)
  expect_equal(length(ggplot2::ggplot_build(original)$plot$guides$aesthetics), 1L)
  expect_length(ggplot2::ggplot_build(clipped)$plot$guides$aesthetics, 0L)
})

test_that("clipped grobs retain the inset boundary when nested clipping is off", {
  skip_if_not_installed("ragg")
  child <- grid::grobTree(
    grid::rectGrob(x = -0.2, width = 1.5,
      gp = grid::gpar(fill = "white", col = NA)),
    vp = grid::viewport(clip = "off"))
  fixture <- inset_integrity_fixture(child)
  baseline <- inset_integrity_body_pixels(fixture$base)
  expect_equal(
    inset_integrity_body_pixels(fixture$base + fixture$inset("on")), baseline)
  expect_lt(
    inset_integrity_body_pixels(fixture$base + fixture$inset("off")), baseline)
})

test_that("inset clipping preserves normalized viewport collision checks", {
  fixture <- inset_integrity_fixture(grid::rectGrob())
  expect_error(ggplot2::ggplotGrob(
    fixture$base + fixture$inset() + fixture$inset()), "Inset viewport collision")
  expect_error(ggplot2::ggplotGrob(
    fixture$base + fixture$inset(hjust = 1)), "Inset viewport collision")
})
