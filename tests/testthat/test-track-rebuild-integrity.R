test_that("track reconstruction keeps interleaved native and domain layer order", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = 50, Value = 2, ID = "gene")
  p <- ggideogram(k, show_names = FALSE, orientation = "horizontal",
    tracks = list(s = geom_track(data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value),
      geom = ggplot2::geom_point(colour = "red")))) +
    ggplot2::geom_point(data = data.frame(x = 5, y = 1),
      ggplot2::aes(x, y), colour = "green") +
    geom_locus(data = d, ggplot2::aes(chr = Chr, position = Pos, label = ID),
      geom = "text", colour = "black")
  colours <- function(plot) unname(vapply(Filter(function(layer) {
    inherits(layer$geom, "GeomPoint") || inherits(layer$geom, "GeomText")
  }, plot$layers), function(layer) layer$aes_params$colour, character(1)))
  q <- p + geom_track(track = "s", geom = ggplot2::geom_point(colour = "blue"))
  expect_identical(colours(q), c("red", "green", "black", "blue"))
  r <- q + geom_track(track = "s", width = 4, label = "Signal", axis = TRUE)
  expect_identical(utils::head(colours(r), 4), colours(q))
  expect_no_warning(ggplot2::ggplotGrob(r))
  expect_identical(colours(p), c("red", "green", "black"))

  regions <- data.frame(Chr = "A", Start = 21, End = 80)
  p <- ggideogram(k, show_names = FALSE,
    tracks = list(body = geom_track(side = "overlay"))) +
    geom_chr(component = "fill", track = "body", data = regions,
      ggplot2::aes(chr = Chr, start = Start, end = End), fill = "pink") +
    geom_track(track = "body", data = d,
      layers = list(geom_locus(ggplot2::aes(chr = Chr, position = Pos)),
        geom_locus(ggplot2::aes(chr = Chr, position = Pos, label = ID), geom = "text")))
  fill <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomIdeogramRect"), logical(1)))
  marks <- which(vapply(p$layers, function(layer) inherits(layer$stat, "StatChrMarker"), logical(1)))
  expect_true(all(fill < marks))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("track reconstruction retains chromosome metadata and attachment registries", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  meta <- data.frame(Chr = c("A", "B"), Group = c("control", "treated"))
  loci <- data.frame(Chr = "A", Pos = 20, ID = "gene_a")
  p <- ggideogram(k) + chr_data(meta, name = "meta") +
    locus_data(loci, name = "loci") +
    chr_data(data.frame(Chr = c("A", "A"), Batch = c(1, 2)),
      name = "batches", relationship = "many-to-many")
  registry <- lapply(c("meta", "batches"), function(name) get_ideogram_attachment(p, name))
  locus <- get_ideogram_attachment(p, "loci", type = "locus")
  q <- p + geom_track(track = "signal", data = data.frame(Chr = "A", Pos = 30, Value = 2),
    ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_point())
  q <- q + geom_track(track = "signal", width = 4)
  for (i in seq_along(registry)) {
    expect_identical(get_ideogram_attachment(q, c("meta", "batches")[i]), registry[[i]])
  }
  expect_identical(get_ideogram_attachment(q, "loci", type = "locus"), locus)
  expect_identical(q$coordinates$layout$data$karyotype$Group, meta$Group)
  expect_identical(q$coordinates$layout$chrom$Group, meta$Group)
  expect_no_warning(ggplot2::ggplotGrob(q))
})

test_that("track updates retain effective native multi-scale prototypes", {
  skip_if_not_installed("ggnewscale")
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 50, 80), Value = c(1, 2, 3),
    Group = c("g1", "g2", "g1"))
  p <- ggideogram(k, show_names = FALSE, tracks = list(one = geom_track(data = d,
    ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Group),
    geom = ggplot2::geom_point(size = 3)))) +
    ggplot2::scale_colour_manual(values = c(g1 = "red", g2 = "blue")) +
    ggnewscale::new_scale_colour()
  p <- p + geom_track(track = "two", data = d,
    ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Value),
    geom = ggplot2::geom_point(size = 2)) +
    ggplot2::scale_colour_gradient(low = "black", high = "yellow")
  snapshot <- ggplot2::ggplot_build(p)
  q <- p + geom_track(track = "one", width = 4)
  result <- ggplot2::ggplot_build(q)
  for (i in seq_along(p$layers)) {
    fields <- intersect(c("colour", "colour_ggnewscale_1", "position", "value"),
      names(snapshot$data[[i]]))
    expect_equal(result$data[[i]][fields], snapshot$data[[i]][fields])
    expect_identical(names(q$layers[[i]]$mapping), names(p$layers[[i]]$mapping))
  }
  expect_no_warning(ggplot2::ggplotGrob(q))
  q <- q + ggnewscale::new_scale_colour() +
    geom_track(track = "three", data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Group),
      geom = ggplot2::geom_point()) +
    ggplot2::scale_colour_manual(values = c(g1 = "purple", g2 = "orange"))
  expect_no_warning(ggplot2::ggplotGrob(q + geom_track(track = "two", width = 3)))

  p <- ggideogram(k, show_names = FALSE, tracks = list(one = geom_track(data = d,
    ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Group),
    geom = ggplot2::geom_point()))) +
    ggplot2::scale_colour_manual(values = c(g1 = "red", g2 = "blue")) +
    ggnewscale::new_scale_colour() +
    ggplot2::geom_point(data = data.frame(x = 4, y = 2, z = 1),
      ggplot2::aes(x, y, colour = z)) +
    ggplot2::scale_colour_gradient(low = "black", high = "yellow")
  expect_no_warning(ggplot2::ggplotGrob(p + geom_track(track = "one", width = 4)))
})
