track_label_fixture <- function(orientation, limits = c(0, 10)) {
  k <- data.frame(Genome = c("wild", "rice"), Assembly = c("v2", "v1"),
    Chr = "1", Start = c(10, 0), End = c(90, 100))
  model <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start,
    end = End, genome = Genome, assembly = Assembly))
  keys <- model$karyotype$.chr
  data <- data.frame(Chr = rep(keys, each = 2), Pos = rep(c(30, 60), 2),
    Value = c(1, 2, 4, 5), Label = paste0("gene", 1:4), ID = 1:4)
  args <- list(data = model, orientation = orientation, reverse_chr = keys[2],
    max_chr_length = 35, show_names = FALSE,
    tracks = track_layout(labels = track("right", width = 2, limits = limits)))
  if (orientation != "circular") args$ncol <- 1
  list(base = do.call(ggideogram, args), data = data)
}

track_label_layer <- function(data, geom, position = "identity") {
  geom_chr_track(data = data, geom = geom, track = "labels", position = position,
    ggplot2::aes(chr = Chr, position = Pos, value = Value,
      label = Label, group = ID), size = 2.5)
}

track_label_nodes <- function(grob) {
  result <- list(grob)
  children <- c(grob$grobs, if (!is.null(grob$children)) as.list(grob$children) else list())
  for (child in children) result <- c(result, track_label_nodes(child))
  result
}

track_label_tree <- function(plot) {
  trees <- Filter(function(g) inherits(g, "textrepeltree") || inherits(g, "labelrepeltree"),
    track_label_nodes(ggplot2::ggplotGrob(plot)))
  stopifnot(length(trees) == 1L)
  trees[[1]]
}

track_label_expected <- function(plot, data) {
  point <- project_chr_track(plot$coordinates$layout, data,
    "Chr", "Pos", "Value", "labels")
  b <- ggplot2::ggplot_build(plot)
  plot$coordinates$transform(data.frame(x = point$.x, y = point$.y),
    b$layout$panel_params[[1]])
}

test_that("native text Position uses the final track layout and retains source fields", {
  for (orientation in c("vertical", "horizontal", "circular")) {
    fixture <- track_label_fixture(orientation)
    p <- fixture$base + track_label_layer(fixture$data, ggplot2::geom_text,
      ggplot2::position_nudge(x = 5, y = 0.5))
    layer <- utils::tail(p$layers, 1)[[1]]
    b <- ggplot2::ggplot_build(p)
    data <- utils::tail(b$data, 1)[[1]]
    point <- p$coordinates$transform(data, b$layout$panel_params[[1]])
    expected <- track_label_expected(p,
      transform(fixture$data, Pos = Pos + 5, Value = Value + 0.5))
    expect_s3_class(layer$geom, "GeomText")
    expect_equal(point$x, expected$x)
    expect_equal(point$y, expected$y)
    expect_equal(data$position, fixture$data$Pos)
    expect_equal(data$value, fixture$data$Value)
    expect_identical(layer$data$ID, fixture$data$ID)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that("third-party text receives standard x/y without unrequested drawing nudges", {
  skip_if_not_installed("ggrepel")
  for (orientation in c("vertical", "horizontal", "circular")) {
    fixture <- track_label_fixture(orientation)
    p <- fixture$base + track_label_layer(fixture$data, ggrepel::geom_text_repel)
    layer <- utils::tail(p$layers, 1)[[1]]
    tree <- track_label_tree(p)
    expected <- track_label_expected(p, fixture$data)
    expect_s3_class(layer$geom, "GeomTextRepel")
    expect_equal(tree$data$x, expected$x)
    expect_equal(tree$data$y, expected$y)
    expect_equal(tree$data$nudge_x, rep(0, nrow(fixture$data)))
    expect_equal(tree$data$nudge_y, rep(0, nrow(fixture$data)))
    expect_identical(layer$data$ID, fixture$data$ID)
  }
  fixture <- track_label_fixture("horizontal")
  p <- fixture$base + track_label_layer(fixture$data, ggrepel::geom_label_repel)
  tree <- track_label_tree(p)
  expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, "GeomLabelRepel")
  expect_equal(tree$data$nudge_x, rep(0, nrow(fixture$data)))
  expect_equal(tree$data$nudge_y, rep(0, nrow(fixture$data)))
})

test_that("later shared-range training changes label projection through the final layout", {
  skip_if_not_installed("ggrepel")
  fixture <- track_label_fixture("vertical", limits = NULL)
  text <- track_label_layer(fixture$data, ggrepel::geom_text_repel)
  wide <- data.frame(Chr = fixture$data$Chr[c(1, 3)], Pos = 50, Value = c(0, 10))
  point <- geom_track_point(data = wide, track = "labels",
    ggplot2::aes(chr = Chr, position = Pos, value = Value))
  early <- fixture$base + text
  after <- early + point
  reverse <- fixture$base + point + text
  first <- track_label_tree(early)$data
  final <- track_label_tree(after)$data
  reordered <- track_label_tree(reverse)$data
  expected <- track_label_expected(after, fixture$data)
  expect_false(isTRUE(all.equal(first$x, final$x)))
  expect_equal(final$x, expected$x)
  expect_equal(final$y, expected$y)
  expect_equal(final$x, reordered$x)
  expect_equal(final$y, reordered$y)
  expect_equal(final$nudge_x, rep(0, nrow(fixture$data)))
  expect_equal(final$nudge_y, rep(0, nrow(fixture$data)))
  expect_equal(after$coordinates$layout$track_ranges[["track:labels"]]$limits, c(0, 10))
  expect_equal(utils::tail(ggplot2::ggplot_build(reverse)$data, 1)[[1]]$position, fixture$data$Pos)
})

test_that("explicit repel Position preserves original anchors and only the requested nudge", {
  skip_if_not_installed("ggrepel")
  for (orientation in c("vertical", "horizontal", "circular")) {
    fixture <- track_label_fixture(orientation)
    p <- fixture$base + track_label_layer(fixture$data, ggrepel::geom_text_repel,
      ggrepel::position_nudge_repel(x = 5, y = 0.5))
    tree <- track_label_tree(p)
    expected <- track_label_expected(p, fixture$data)
    shifted <- track_label_expected(p,
      transform(fixture$data, Pos = Pos + 5, Value = Value + 0.5))
    expect_equal(tree$data$x, expected$x)
    expect_equal(tree$data$y, expected$y)
    expect_equal(tree$data$nudge_x, shifted$x - expected$x)
    expect_equal(tree$data$nudge_y, shifted$y - expected$y)
    b <- utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]
    expect_equal(b$position, fixture$data$Pos)
    expect_equal(b$value, fixture$data$Value)
  }
})
