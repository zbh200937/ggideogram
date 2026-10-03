staged_track_fixture <- function() {
  values <- c(1, 1.2, 1.5, 2, 2.7, 3, 3.2, 4, 5, 6, 7, 8)
  data <- data.frame(Chr = rep(c("A", "B"), each = 24), Pos = 50,
    Value = rep(values, 4), Class = rep(rep(c("a", "b"), each = 12), 2))
  data$Value[data$Class == "b"] <- data$Value[data$Class == "b"] + 2
  list(karyotype = data.frame(Chr = c("A", "B"), Start = 0, End = 100),
    data = data)
}

staged_track_layer_index <- function(plot, geom) {
  which(vapply(plot$layers, function(layer) inherits(layer$geom, geom), logical(1)))
}

test_that("implicit violin groups defer statistics and scaled aesthetics", {
  fixture <- staged_track_fixture()
  source <- fixture$data
  stage <- ggplot2::stage
  mappings <- list(
    ggplot2::aes(fill = sqrt(ggplot2::after_stat(scaled)),
      colour = ggplot2::after_scale(fill)),
    ggplot2::aes(fill = stage(after_stat = sqrt(scaled)),
      colour = ggplot2::after_scale(fill)))
  for (mapping in mappings) {
    reference_mapping <- combine_track_mapping(mapping,
      ggplot2::aes(x = Pos, y = Value))
    reference <- ggplot2::ggplot_build(
      ggplot2::ggplot(source[source$Chr == "A", ], reference_mapping) +
        ggplot2::geom_violin(width = 5))$data[[1]]
    for (orientation in c("horizontal", "vertical", "circular")) {
      track_mapping <- combine_track_mapping(mapping,
        ggplot2::aes(chr = Chr, x = Pos, y = Value))
      p <- ggideogram(fixture$karyotype, orientation = orientation,
        show_names = FALSE, tracks = list(distribution = geom_track(
          data = source, mapping = track_mapping,
          geom = ggplot2::geom_violin(width = 5))))
      expect_no_warning(build <- ggplot2::ggplot_build(p))
      actual <- build$data[[staged_track_layer_index(p, "GeomViolin")]]
      expect_equal(length(unique(actual$group)), 2L)
      for (group in unique(actual$group)) {
        rows <- actual[actual$group == group, ]
        # Projection may reverse the density grid; its native statistical
        # values and colour scale must still agree with the ordinary plot.
        expect_equal(sort(rows$scaled), sort(reference$scaled), tolerance = 1e-9)
        expect_equal(sort(rows$fill), sort(reference$fill))
        expect_equal(rows$colour, rows$fill)
      }
      expect_no_warning(ggplot2::ggplotGrob(p))
    }
  }
  expect_identical(fixture$data, source)
})

test_that("stage start preserves implicit discrete groups before native styling", {
  fixture <- staged_track_fixture()
  opacity <- .6
  stage <- ggplot2::stage
  mapping <- ggplot2::aes(
    fill = stage(start = Class,
      after_scale = scales::alpha(fill, opacity)),
    colour = ggplot2::after_scale(fill))
  for (geom in c("violin", "boxplot")) {
    native <- switch(geom, violin = ggplot2::geom_violin(width = 5),
      boxplot = ggplot2::geom_boxplot(width = 5))
    reference <- ggplot2::ggplot_build(
      ggplot2::ggplot(fixture$data[fixture$data$Chr == "A", ],
        combine_track_mapping(mapping, ggplot2::aes(x = Pos, y = Value))) +
        native)$data[[1]]
    for (orientation in c("horizontal", "vertical", "circular")) {
      p <- ggideogram(fixture$karyotype, orientation = orientation,
        show_names = FALSE, tracks = list(distribution = geom_track(
          data = fixture$data,
          mapping = combine_track_mapping(mapping,
            ggplot2::aes(chr = Chr, x = Pos, y = Value)), geom = native)))
      expect_no_warning(build <- ggplot2::ggplot_build(p))
      index <- staged_track_layer_index(p,
        if (geom == "violin") "GeomViolin" else "GeomBoxplot")
      actual <- build$data[[index]]
      expect_equal(length(unique(actual$group)), 4L)
      expect_equal(nrow(actual), 2L * nrow(reference))
      expect_setequal(actual$fill, reference$fill)
      expect_setequal(actual$colour, reference$colour)
      expect_no_warning(ggplot2::ggplotGrob(p))
    }
  }
})

test_that("explicit native groups still accept combined delayed mappings", {
  fixture <- staged_track_fixture()
  p <- ggideogram(fixture$karyotype, show_names = FALSE,
    tracks = list(distribution = geom_track(data = fixture$data,
      ggplot2::aes(chr = Chr, x = Pos, y = Value, group = Class,
        fill = ggplot2::after_stat(scaled), colour = ggplot2::after_scale(fill)),
      geom = ggplot2::geom_violin(width = 5))))
  expect_no_warning(build <- ggplot2::ggplot_build(p))
  actual <- build$data[[staged_track_layer_index(p, "GeomViolin")]]
  expect_equal(length(unique(actual$group)), 4L)
  expect_equal(actual$colour, actual$fill)
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("distribution widths and bandwidths retain source units", {
  k <- data.frame(Chr = "A", Start = 0, End = 1e6)
  d <- data.frame(Chr = "A", Pos = 5e5,
    Value = c(seq(1, 3, length.out = 15), seq(10, 15, length.out = 15)))
  native <- ggplot2::geom_violin(width = 1e5, bw = 1, trim = FALSE)
  reference <- ggplot2::ggplot_build(ggplot2::ggplot(d,
    ggplot2::aes(x = Pos, y = Value)) + native)$data[[1]]
  for (orientation in c("horizontal", "vertical", "circular")) {
    for (geom in list(native, ggplot2::geom_boxplot(width = 1e5))) {
      p <- ggideogram(k, orientation = orientation, reverse_chr = "A", show_names = FALSE,
        tracks = list(dist = geom_track(data = d,
          ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = geom)))
      built <- ggplot2::ggplot_build(p)$data[[staged_track_layer_index(p, class(geom$geom)[1])]]
      expected <- project_chr_point(p, data.frame(Chr = "A", Pos = c(45e4, 55e4)), "Chr", "Pos")
      axis <- if (orientation == "vertical") "y" else "x"
      expect_equal(range(built[[paste0(axis, "min")]], built[[paste0(axis, "max")]]),
        range(expected[[paste0(".", axis)]]))
      if (inherits(geom$geom, "GeomViolin")) {
        expect_equal(built$scaled, reference$scaled, tolerance = 1e-10)
        expect_equal(p$coordinates$layout$track_ranges[["track:dist"]]$limits, range(reference$y))
      }
      expect_no_warning(ggplot2::ggplotGrob(p))
    }
  }
})

test_that("box summaries precede nonlinear track transforms", {
  k <- data.frame(Chr = "A", Start = 0, End = 1e6)
  d <- data.frame(Chr = "A", Pos = 5e5, Value = c(1, 2, 4, 8, 16, 32))
  p <- ggideogram(k, orientation = "horizontal", show_names = FALSE,
    tracks = list(dist = geom_track(data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value), transform = "sqrt",
      geom = ggplot2::geom_boxplot(width = 1e5))))
  built <- ggplot2::ggplot_build(p)$data[[staged_track_layer_index(p, "GeomBoxplot")]]
  q <- stats::quantile(d$Value, c(.25, .5, .75))
  expected <- project_chr_track(p, data.frame(Chr = "A", Pos = 5e5, Value = q),
    "Chr", "Pos", "Value", "dist")
  expect_equal(unname(unlist(built[c("lower", "middle", "upper")])), expected$.y)
})
