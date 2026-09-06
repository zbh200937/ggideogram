test_that("track axes follow ranges registered after axis addition", {
  kar <- data.frame(Chr = "A", Start = 0, End = 100)
  values <- data.frame(Chr = "A", Pos = c(20, 80), Value = c(0, 10))
  for (orientation in c("vertical", "horizontal")) {
    p <- ggideogram(kar, orientation = orientation, show_names = FALSE,
      tracks = track_layout(signal = track())) +
      geom_track_point(ggplot2::aes(chr = Chr, position = Pos, value = Value),
        values, track = "signal") +
      geom_track_axis("signal", breaks = c(0, 10))
    updated <- p + geom_track_point(
      ggplot2::aes(chr = Chr, position = Pos, value = Value),
      transform(values, Value = c(0, 100)), track = "signal")
    built <- ggplot2::ggplot_build(updated)
    text <- built$data[[length(built$data) - 1L]]
    expected <- project_chr_track(updated$coordinates$layout,
      transform(values, Pos = 100), "Chr", "Pos", "Value", "signal")
    expect_equal(text$label, c("0", "10"))
    expect_equal(text$x, expected$.x)
    expect_equal(text$y, expected$.y)
    expect_no_error(ggplot2::ggplotGrob(updated))
    # Building a derived plot must not change the original plot's range.
    original <- ggplot2::ggplot_build(p)
    old_text <- original$data[[length(original$data)]]
    expect_false(isTRUE(all.equal(old_text[c("x", "y")], text[c("x", "y")])))
  }
})

test_that("distribution tracks preserve discrete aesthetic groups", {
  kar <- data.frame(Chr = "A", Start = 0, End = 100)
  values <- data.frame(Chr = "A", Pos = 50,
    Value = c(1, 2, 3, 7, 8, 9), Class = rep(c("a", "b"), each = 3))
  base <- ggideogram(kar, tracks = track_layout(signal = track(limits = c(0, 10))))
  p <- base + geom_track_boxplot(
    ggplot2::aes(chr = Chr, position = Pos, value = Value, fill = Class),
    values, track = "signal")
  expect_warning(built <- ggplot2::ggplot_build(p), NA)
  boxes <- built$data[[length(built$data)]]
  expect_equal(nrow(boxes), 2L)
  expect_equal(length(unique(boxes$fill)), 2L)
  expect_equal(sort(boxes$xmiddle), c(1.4, 2.0))
  violin <- base + geom_track_violin(
    ggplot2::aes(chr = Chr, position = Pos, value = Value, fill = Class),
    values, track = "signal")
  expect_warning(built <- ggplot2::ggplot_build(violin), NA)
  expect_equal(length(unique(built$data[[length(built$data)]]$group)), 2L)
})

test_that("inset collision checks span layers and reset on each build", {
  kar <- data.frame(Chr = "A", Start = 0, End = 100)
  data <- data.frame(Chr = "A", Pos = 50)
  data$Plot <- list(grid::rectGrob())
  inset <- function(overlap = "error", pos = 50) geom_chr_inset(
    ggplot2::aes(chr = Chr, position = Pos, plot = Plot),
    transform(data, Pos = pos), width = grid::unit(0.1, "npc"),
    height = grid::unit(0.1, "npc"), overlap = overlap)
  base <- ggideogram(kar)
  expect_error(ggplot2::ggplotGrob(base + inset() + inset()), "collision")
  expect_no_error(ggplot2::ggplotGrob(base + inset("allow") + inset("allow")))
  separate <- base + inset(pos = 20) + inset(pos = 80)
  expect_no_error(ggplot2::ggplotGrob(separate))
  built <- ggplot2::ggplot_build(separate)
  expect_no_error(ggplot2::ggplot_gtable(built))
  expect_no_error(ggplot2::ggplot_gtable(built))
})

test_that("inset intervals validate both genomic endpoints", {
  kar <- data.frame(Chr = "A", Start = 0, End = 100)
  for (ends in list(c(-50, 150), c(-10, 50), c(50, 110))) {
    data <- data.frame(Chr = "A", Start = ends[1], End = ends[2])
    data$Plot <- list(grid::rectGrob())
    p <- ggideogram(kar) + geom_chr_inset(
      ggplot2::aes(chr = Chr, start = Start, end = End, plot = Plot), data,
      width = grid::unit(0.1, "npc"), height = grid::unit(0.1, "npc"))
    expect_error(ggplot2::ggplotGrob(p), "outside|range")
  }
})
