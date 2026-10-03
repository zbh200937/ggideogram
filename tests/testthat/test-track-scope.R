track_scope_fixture <- function() {
  k <- data.frame(Genome = c("rice", "wild"), Assembly = c("IRGSP", "W1943"),
    Chr = "1", Start = 0, End = c(100, 90))
  model <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End,
    genome = Genome, assembly = Assembly))
  d <- data.frame(Chr = rep(model$karyotype$.chr, each = 3),
    Pos = rep(c(20, 50, 80), 2), Value = c(1, 3, 2, 4, 6, 5),
    Series = rep(c("rice", "wild"), each = 3), ID = paste0("locus", 1:6))
  list(model = model, data = d, keys = model$karyotype$.chr)
}

track_scope_indices <- function(plot, geom) {
  which(vapply(plot$layers, function(layer) inherits(layer$geom, geom), logical(1)))
}

track_scope_xy <- function(build, index) {
  build$plot$coordinates$transform(build$data[[index]], build$layout$panel_params[[1]])
}

track_scope_expected <- function(plot, data, track) {
  point <- project_chr_track(plot, data, "Chr", "Pos", "Value", track)
  build <- ggplot2::ggplot_build(plot)
  plot$coordinates$transform(data.frame(x = point$.x, y = point$.y),
    build$layout$panel_params[[1]])
}

test_that("named track scopes share data and preserve local native layer settings", {
  f <- track_scope_fixture()
  source <- f$data
  local <- transform(f$data[c(1, 4), ], Alt = c(8, 10), Category = c("low", "high"))
  glyph <- function(data, params, size) grid::circleGrob(r = .3,
    gp = grid::gpar(col = "purple", fill = NA))
  native <- ggplot2::geom_point(ggplot2::aes(y = Alt, colour = Category), data = local,
    position = ggplot2::position_nudge(x = 5, y = -.25), size = 2.4, stroke = .31,
    show.legend = FALSE, key_glyph = glyph)
  scopes <- list(
    signal = geom_track(side = "left", width = 1.3, data = f$data,
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Series),
      layers = list(ggplot2::geom_line(linewidth = .37), native)),
    counts = geom_track(side = "right", width = 1.1, data = f$data,
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value, fill = Series),
      layers = list(ggplot2::geom_col(width = 4,
        position = ggplot2::position_dodge2(width = 3, preserve = "single"),
        show.legend = c(fill = TRUE))))
  )
  p <- ggideogram(f$model, tracks = scopes, ncol = 1, show_names = FALSE,
    reverse_chr = f$keys[2]) +
    ggplot2::scale_colour_manual(values = c(rice = "red", wild = "blue",
      low = "orange", high = "purple")) +
    ggplot2::scale_fill_manual(values = c(rice = "red", wild = "blue"))
  b <- ggplot2::ggplot_build(p)
  line <- track_scope_indices(p, "GeomLine")[1]
  point <- track_scope_indices(p, "GeomPoint")[1]
  col <- track_scope_indices(p, "GeomCol")[1]
  expected <- track_scope_expected(p, f$data, "signal")
  actual <- track_scope_xy(b, line)
  expect_equal(actual$x, expected$x)
  expect_equal(actual$y, expected$y)
  expected <- track_scope_expected(p,
    transform(local, Pos = Pos + 5, Value = Alt - .25), "signal")
  actual <- track_scope_xy(b, point)
  expect_equal(actual$x, expected$x)
  expect_equal(actual$y, expected$y)
  expect_equal(b$data[[point]]$position, local$Pos)
  expect_equal(b$data[[point]]$value, local$Alt)
  expect_identical(p$layers[[point]]$data$ID, local$ID)
  expect_equal(p$layers[[point]]$aes_params$size, 2.4)
  expect_equal(p$layers[[point]]$aes_params$stroke, .31)
  expect_false(p$layers[[point]]$show.legend)
  key <- p$layers[[point]]$geom$draw_key(list(), list(), 5)
  expect_s3_class(key, "circle")
  expect_equal(key$gp$col, "purple")
  expect_equal(p$layers[[col]]$position$width, 3)
  expect_equal(p$layers[[col]]$position$preserve, "single")
  expect_equal(p$layers[[col]]$show.legend, c(fill = TRUE))
  expect_identical(f$data, source)
  expect_equal(native$data, local)
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("nine supported native geoms share chromosome scopes in all orientations", {
  f <- track_scope_fixture()
  d <- f$data
  d$Low <- d$Value - .3
  d$High <- d$Value + .4
  distribution <- do.call(rbind, lapply(c(-.2, 0, .2, .4),
    function(delta) transform(d, Value = Value + delta)))
  distribution$Bin <- interaction(distribution$Chr, distribution$Pos)
  native <- list(
    point = ggplot2::geom_point(size = 1.5),
    line = ggplot2::geom_line(linewidth = .3),
    col = ggplot2::geom_col(width = 4),
    area = ggplot2::geom_area(stat = "identity", position = "identity",
      outline.type = "both"),
    ribbon = ggplot2::geom_ribbon(ggplot2::aes(ymin = Low, ymax = High),
      outline.type = "both"),
    boxplot = ggplot2::geom_boxplot(ggplot2::aes(group = Bin),
      data = distribution, width = 4, coef = 2, outlier.shape = 23),
    violin = ggplot2::geom_violin(ggplot2::aes(group = Bin),
      data = distribution, width = 4, adjust = .7),
    tile = ggplot2::geom_tile(width = 4, height = .1),
    text = ggplot2::geom_text(ggplot2::aes(label = ID), size = 2)
  )
  classes <- c(point = "GeomPoint", line = "GeomLine", col = "GeomCol",
    area = "GeomArea", ribbon = "GeomRibbon", boxplot = "GeomBoxplot",
    violin = "GeomViolin", tile = "GeomTile", text = "GeomText")
  for (orientation in c("vertical", "horizontal", "circular")) {
    scopes <- lapply(names(native), function(name) geom_track(side = "right",
      width = .9, data = d,
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value, group = Series),
      layers = list(native[[name]])))
    names(scopes) <- names(native)
    p <- ggideogram(f$model, orientation = orientation, tracks = scopes,
      reverse_chr = f$keys[2], show_names = FALSE, radius = 30)
    for (name in names(classes)) {
      expect_true(any(vapply(p$layers,
        function(layer) inherits(layer$geom, classes[[name]]), logical(1))),
        info = paste(orientation, name))
    }
    box <- p$layers[[track_scope_indices(p, "GeomBoxplot")[1]]]
    violin <- p$layers[[track_scope_indices(p, "GeomViolin")[1]]]
    expect_s3_class(box$stat, "StatBoxplot")
    expect_equal(box$stat_params$coef, 2)
    expect_s3_class(violin$stat, "StatYdensity")
    expect_equal(violin$stat_params$adjust, .7)
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
})

test_that("native area preflight includes values interpolated at staggered samples", {
  f <- track_scope_fixture()
  d <- data.frame(Chr = rep(f$keys, each = 6),
    Pos = rep(c(20, 40, 60, 30, 50, 70), 2),
    Value = rep(c(1, 2, 3, 2, 3, 4), 2),
    Series = rep(rep(c("one", "two"), each = 3), 2))
  source <- d
  native <- ggplot2::geom_area()
  mapping <- ggplot2::aes(chr = Chr, x = Pos, y = Value, fill = Series, group = Series)
  trained <- native_range_track_values(native, d, mapping)
  expect_equal(range(trained$value), c(0, 6.5))
  expect_setequal(trained$chr, f$keys)
  for (orientation in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(f$model, orientation = orientation, reverse_chr = f$keys[2],
      show_names = FALSE, tracks = list(area = geom_track(data = d,
        mapping = mapping, layers = list(native))))
    expect_equal(p$coordinates$layout$track_ranges[["track:area"]]$limits, c(0, 6.5))
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(d, source)
})

test_that("native area alignment and stacking precede the final shared projection", {
  f <- track_scope_fixture()
  first <- transform(f$data, Series = "one")
  second <- transform(f$data, Series = "two", Pos = Pos + 5, Value = Value + .7)
  d <- rbind(first, second)
  source <- d
  native <- ggplot2::geom_area(outline.type = "upper", colour = "black", linewidth = .29)
  params <- native$geom_params
  palette <- ggplot2::scale_fill_manual(values = c(one = "red", two = "blue"))
  references <- lapply(f$keys, function(chr) ggplot2::ggplot_build(
    ggplot2::ggplot(d[d$Chr == chr, ], ggplot2::aes(Pos, Value, fill = Series, group = Series)) +
      native + palette)$data[[1]])
  wide <- data.frame(Chr = f$keys, Pos = 50, Value = 20)
  for (orientation in c("vertical", "horizontal", "circular")) {
    before <- ggideogram(f$model, orientation = orientation, reverse_chr = f$keys[2],
      show_names = FALSE, tracks = list(area = geom_track(side = "right", width = 1.8,
        data = d, mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value,
          fill = Series, group = Series), layers = list(native)))) + palette
    expected_limits <- range(unlist(lapply(references,
      function(data) c(data$ymin, data$ymax))))
    expect_equal(before$coordinates$layout$track_ranges[["track:area"]]$limits,
      expected_limits)
    expect_no_warning(ggplot2::ggplotGrob(before))
    after <- before + geom_track(track = "area", data = wide,
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value),
      layers = list(ggplot2::geom_point(size = 1.3, show.legend = FALSE)))
    index <- track_scope_indices(after, "GeomArea")[1]
    layer <- after$layers[[index]]
    expect_s3_class(layer$geom, "GeomArea")
    expect_s3_class(layer$stat, "StatAlign")
    expect_s3_class(layer$position, "PositionStack")
    expect_equal(layer$geom_params$outline.type, "upper")
    expect_equal(layer$aes_params$linewidth, .29)
    built <- ggplot2::ggplot_build(after)$data[[index]]
    expect_true(all(built$ideogram_projected))
    for (j in seq_along(f$keys)) {
      actual <- built[built$ideogram_chr == f$keys[j], ]
      raw <- references[[j]]
      expect_equal(nrow(actual), nrow(raw))
      expected <- function(value) project_chr_track(after,
        data.frame(Chr = f$keys[j], Pos = raw$x, Value = value),
        "Chr", "Pos", "Value", "area")
      lower <- expected(raw$ymin)
      upper <- expected(raw$ymax)
      center <- expected(raw$y)
      expect_equal(actual$position, raw$x)
      expect_equal(actual$x, center$.x)
      expect_equal(actual$y, center$.y)
      if (orientation == "vertical") {
        expect_equal(actual$xmin, lower$.x)
        expect_equal(actual$xmax, upper$.x)
      } else {
        expect_equal(actual$ymin, lower$.y)
        expect_equal(actual$ymax, upper$.y)
      }
    }
    expect_no_warning(ggplot2::ggplotGrob(after))
  }
  expect_identical(d, source)
  expect_equal(native$geom_params, params)
  expect_s3_class(native$geom, "GeomArea")
  expect_s3_class(native$stat, "StatAlign")
})

test_that("native ribbon bounds retain their original upper and lower semantics", {
  f <- track_scope_fixture()
  d <- transform(f$data, Low = Value - .3, High = Value + .4)
  native <- ggplot2::geom_ribbon(outline.type = "upper", colour = "black", linewidth = .27)
  for (orientation in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(f$model, orientation = orientation, reverse_chr = f$keys[2],
      show_names = FALSE, tracks = list(ribbon = geom_track(side = "left", width = 1.7,
        reverse = TRUE, data = d,
        mapping = ggplot2::aes(chr = Chr, x = Pos, ymin = Low, ymax = High, group = Series),
        layers = list(native))))
    index <- track_scope_indices(p, "GeomRibbon")[1]
    layer <- p$layers[[index]]
    expect_s3_class(layer$geom, "GeomRibbon")
    expect_s3_class(layer$stat, "StatIdentity")
    expect_s3_class(layer$position, "PositionIdentity")
    expect_equal(layer$geom_params$outline.type, "upper")
    built <- ggplot2::ggplot_build(p)$data[[index]]
    expect_true(all(built$ideogram_projected))
    lower <- project_chr_track(p, d, "Chr", "Pos", "Low", "ribbon")
    upper <- project_chr_track(p, d, "Chr", "Pos", "High", "ribbon")
    expect_equal(built$x, lower$.x)
    expect_equal(built$y, lower$.y)
    if (orientation == "vertical") {
      expect_equal(built$xmin, lower$.x)
      expect_equal(built$xmax, upper$.x)
    } else {
      expect_equal(built$ymin, lower$.y)
      expect_equal(built$ymax, upper$.y)
    }
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
})

test_that("statistical scopes derive fixed limits from all layers before projection", {
  f <- track_scope_fixture()
  distribution <- data.frame(Chr = rep(f$keys, each = 4), Pos = 50,
    Value = rep(c(1, 2, 4, 8), 2), ID = rep(1:4, 2))
  wide <- data.frame(Chr = rep(f$keys, each = 2), Pos = rep(c(30, 70), 2),
    Value = rep(c(0, 20), 2))
  box <- ggplot2::geom_boxplot(width = 4, coef = 2)
  point <- ggplot2::geom_point(data = wide, size = 1.6)
  make_plot <- function(layers) ggideogram(f$model, orientation = "horizontal",
    show_names = FALSE, tracks = list(dist = geom_track(side = "right", width = 2,
      data = distribution, mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value),
      layers = layers)))
  p <- make_plot(list(box, point))
  q <- make_plot(list(point, box))
  build <- ggplot2::ggplot_build(p)
  reverse <- ggplot2::ggplot_build(q)
  a <- build$data[[track_scope_indices(p, "GeomBoxplot")[1]]]
  b <- reverse$data[[track_scope_indices(q, "GeomBoxplot")[1]]]
  expected <- project_chr_track(p, data.frame(Chr = f$keys, Pos = 50, Value = 3),
    "Chr", "Pos", "Value", "dist")
  expect_equal(a$middle, expected$.y)
  expect_equal(a$middle, b$middle)
  expect_equal(a$x, b$x)
  expect_no_warning(ggplot2::ggplotGrob(p))
  before <- a
  outside <- data.frame(Chr = f$keys[1], Pos = 50, Value = 100)
  expect_error(p + geom_track(track = "dist", data = outside,
    mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value),
    layers = list(ggplot2::geom_point())))
  expect_equal(ggplot2::ggplot_build(p)$data[[track_scope_indices(p, "GeomBoxplot")[1]]],
    before)
})

track_scope_arabidopsis <- function() {
  path <- system.file("extdata", "arabidopsis-first-genes.gff3", package = "ggideogram")
  features <- read_chr_features(path)
  transcripts <- features[features$Type == "mRNA", ]
  features$Transcript <- ifelse(features$Type == "mRNA", features$ID, features$Parent)
  models <- features[features$Transcript %in% transcripts$ID, ]
  models$Gene <- transcripts$Parent[match(models$Transcript, transcripts$ID)]
  list(features = features, models = models,
    karyotype = data.frame(Chr = "1", Start = 0, End = 30427671))
}

test_that("complete source GFF data feeds windowed native and domain scopes together", {
  f <- track_scope_arabidopsis()
  source <- f$models
  exons <- f$features[f$features$Type == "exon", ]
  counts <- bin_genome(f$features[f$features$Type == "gene", ],
    f$karyotype, window = 1000, method = "count")
  counts$Mid <- (counts$Start + counts$End) / 2
  view <- chr_view(f$karyotype, "1", 7000, 9000)
  scopes <- list(
    genes = geom_track(side = "right", width = 4, data = f$models,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End,
        gene = Gene, type = Type, strand = Strand, fill = Strand),
      layers = list(geom_genemodel(mode = "gene", labels = FALSE))),
    exons = geom_track(side = "overlay", width = .6, data = exons,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, fill = Strand),
      layers = list(geom_chr(component = "fill"))),
    count = geom_track(side = "left", width = 1.5, data = counts,
      mapping = ggplot2::aes(chr = Chr, x = Mid, y = Value),
      layers = list(ggplot2::geom_col(width = 800), ggplot2::geom_point(size = 1.4)))
  )
  for (orientation in c("horizontal", "circular")) {
    p <- ggideogram(view, orientation = orientation, tracks = scopes,
      reverse_chr = "1", show_names = FALSE) +
      ggplot2::scale_fill_manual(values = c("+" = "blue", "-" = "orange"))
    point <- p$layers[[track_scope_indices(p, "GeomPoint")[1]]]
    selected <- counts[counts$Mid >= 7000 & counts$Mid <= 9000, ]
    expect_equal(point$data$Mid, selected$Mid)
    expect_equal(point$data$Value, selected$Value)
    fills <- Filter(function(layer) inherits(layer$geom, "GeomIdeogramRect"), p$layers)
    expect_length(fills, 1)
    keep <- exons$Start - 1 < 9000 & exons$End > 7000
    expect_equal(fills[[1]]$data$.source_start, exons$Start[keep])
    expect_equal(fills[[1]]$data$.source_end, exons$End[keep])
    expect_true(any(vapply(p$layers,
      function(layer) inherits(layer$geom, "GeomSegment"), logical(1))))
    expect_no_warning(ggplot2::ggplotGrob(p))
  }
  expect_identical(f$models, source)
  invalid <- data.frame(Chr = "1", Pos = 30427672, Value = 1)
  expect_error(ggideogram(view, tracks = list(bad = geom_track(data = invalid,
    mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value),
    layers = list(ggplot2::geom_point())))))
})

track_scope_grobs <- function(grob) {
  out <- list(grob)
  children <- c(grob$grobs,
    if (!is.null(grob$children)) as.list(grob$children) else list())
  for (child in children) out <- c(out, track_scope_grobs(child))
  out
}

test_that("third-party scope labels coexist with complete insets and external guides", {
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("patchwork")
  f <- track_scope_fixture()
  d <- f$data[c(1, 4), ]
  child <- ggplot2::ggplot(data.frame(Group = c("low", "high"), N = c(2, 5)),
    ggplot2::aes(Group, N, fill = Group)) + ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(low = "grey70", high = "grey30"),
      name = "Child groups")
  child_before <- ggplot2::ggplot_build(child)$data
  p <- ggideogram(f$model, ncol = 1, reverse_chr = f$keys[2], show_names = FALSE,
    tracks = list(
      labels = geom_track(side = "right", width = 2, data = d,
        mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value,
          label = ID, colour = Series),
        layers = list(ggrepel::geom_text_repel(seed = 11, max.overlaps = Inf,
          size = 2.7), ggplot2::geom_point(size = 1.6))),
      inset = geom_track(side = "right", width = 12))) +
    ggplot2::scale_colour_manual(values = c(rice = "red", wild = "blue"),
      name = "Genome")
  inset <- data.frame(Chr = f$keys[1], Pos = 50)
  inset$Plot <- list(child)
  p <- p + geom_locus_inset(child, data = inset,
    mapping = ggplot2::aes(chr = Chr, position = Pos), track = "inset",
    width = grid::unit(.18, "npc"), height = grid::unit(.15, "npc"), clip = "on")
  trees <- Filter(function(grob) inherits(grob, "textrepeltree"),
    track_scope_grobs(ggplot2::ggplotGrob(p)))
  expect_length(trees, 1)
  expected <- track_scope_expected(p, d, "labels")
  expect_equal(trees[[1]]$data$x, expected$x)
  expect_equal(trees[[1]]$data$y, expected$y)
  expect_equal(trees[[1]]$data$nudge_x, c(0, 0))
  expect_equal(trees[[1]]$data$nudge_y, c(0, 0))
  expect_s3_class(p$layers[[track_scope_indices(p, "GeomTextRepel")[1]]]$geom,
    "GeomTextRepel")
  summary <- ggplot2::ggplot(d, ggplot2::aes(Chr, Value, colour = Series)) +
    ggplot2::geom_point(size = 1.6) + scale_x_chromosome(f$model) +
    ggplot2::scale_colour_manual(values = c(rice = "red", wild = "blue"),
      name = "Genome")
  combined <- patchwork::wrap_plots(p, summary) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  expect_no_warning(patchwork::patchworkGrob(combined))
  expect_equal(ggplot2::ggplot_build(child)$data, child_before)
})

test_that("adding a track reprojects the complete source recipe together", {
  f <- track_scope_fixture()
  models <- data.frame(Chr = rep(f$keys, each = 2),
    Start = rep(c(21, 51), 2), End = rep(c(35, 65), 2),
    Gene = rep(c("g1", "g2"), each = 2), Type = "exon", Strand = c("+", "+", "-", "-"))
  marks <- f$data[c(1, 4), ]
  nodes <- chr_nodes(marks, ggplot2::aes(id = ID, chr = Chr, position = Pos))
  edges <- data.frame(From = marks$ID[1], To = marks$ID[2], Pair = "ortholog")
  child <- ggplot2::ggplot(data.frame(x = 1:3, y = c(1, 4, 2)),
    ggplot2::aes(x, y)) + ggplot2::geom_line(linewidth = .31) +
    ggplot2::labs(x = "Sample", y = "Expression", title = "Inset")
  insets <- data.frame(Chr = f$keys[2], Pos = 50)
  insets$Plot <- list(child)
  scopes <- list(
    signal = geom_track(side = "left", width = 1.4, data = f$data,
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value, colour = Series),
      layers = list(ggplot2::geom_line(linewidth = .37), ggplot2::geom_point(size = 1.6))),
    genes = geom_track(side = "right", width = 2),
    inset = geom_track(side = "right", width = 12))
  extra <- geom_track(side = "left", width = 5, data = f$data,
    mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value, fill = Series),
    layers = list(ggplot2::geom_col(width = 4)), label = "Counts")
  make_plot <- function(tracks) {
    ggideogram(f$model, tracks = tracks, ncol = 2, reverse_chr = f$keys[2],
      show_names = TRUE, axis = TRUE, axis_breaks = c(0, 50, 90)) +
      geom_genemodel(data = models, mode = "gene", track = "genes",
        ggplot2::aes(chr = Chr, start = Start, end = End,
          gene = Gene, type = Type, strand = Strand), labels = FALSE) +
      geom_chrlink(data = edges, nodes = nodes,
        ggplot2::aes(from = From, to = To), linewidth = .23, curvature = 0) +
      geom_locus(data = marks, ggplot2::aes(chr = Chr, position = Pos, colour = Series),
        side = "left", gap = .5, size = 2.1, stroke = .28) +
      geom_locus_inset(child, data = insets,
        mapping = ggplot2::aes(chr = Chr, position = Pos), track = "inset",
        width = grid::unit(.1, "npc"), height = grid::unit(.18, "npc"), clip = "on") +
      ggplot2::scale_colour_manual(values = c(rice = "red", wild = "blue"), name = "Genome") +
      ggplot2::labs(title = "Whole recipe", caption = "Original source coordinates") +
      ggplot2::theme(legend.position = "bottom", text = ggplot2::element_text(size = 9))
  }
  before <- make_plot(scopes)
  snapshot <- ggplot2::ggplot_build(before)
  after <- before + geom_track(track = "counts", side = "left", width = 5,
    data = f$data, mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value, fill = Series),
    layers = list(ggplot2::geom_col(width = 4)), label = "Counts")
  reference <- make_plot(c(scopes, list(counts = extra)))
  built <- ggplot2::ggplot_build(after)
  expected <- ggplot2::ggplot_build(reference)
  selectors <- list(
    native = function(layer) inherits(layer$geom, "GeomLine"),
    gene = function(layer) is.data.frame(layer$data) && ".model_id" %in% names(layer$data),
    link = function(layer) is.data.frame(layer$data) && ".pair_id" %in% names(layer$data),
    marker = function(layer) inherits(layer$stat, "StatChrMarker"),
    inset = function(layer) inherits(layer$geom, "GeomChrInset"))
  fields <- c("x", "y", "xend", "yend", "xmin", "xmax", "ymin", "ymax",
    "chr", "position", "value", "ideogram_chr", "ideogram_position", "label",
    "colour", "fill", "alpha", "size", "linewidth", "shape", "stroke")
  for (name in names(selectors)) {
    select <- selectors[[name]]
    a <- which(vapply(after$layers, select, logical(1)))
    b <- which(vapply(reference$layers, select, logical(1)))
    expect_true(length(a) > 0, info = name)
    expect_length(a, length(b))
    for (j in seq_along(a)) {
      actual <- track_scope_xy(built, a[j])
      target <- track_scope_xy(expected, b[j])
      keep <- intersect(fields, intersect(names(actual), names(target)))
      expect_equal(actual[keep], target[keep], info = paste("Reprojection", name, j))
    }
  }
  a <- track_scope_indices(before, "GeomLine")[1]
  b <- track_scope_indices(after, "GeomLine")[1]
  expect_false(isTRUE(all.equal(track_scope_xy(snapshot, a)$x,
    track_scope_xy(built, b)$x)))
  expect_equal(built$data[[b]]$position, snapshot$data[[a]]$position)
  expect_equal(built$data[[b]]$value, snapshot$data[[a]]$value)
  expect_equal(after$theme$legend.position, "bottom")
  expect_equal(after$theme$text$size, 9)
  expect_equal(after$labels$title, "Whole recipe")
  expect_equal(after$scales$get_scales("colour")$name, "Genome")
  expect_equal(ggplot2::ggplot_build(before)$data[[a]], snapshot$data[[a]])
  expect_no_warning(ggplot2::ggplotGrob(after))
  inset_grobs <- Filter(function(grob) inherits(grob, "gTree") &&
    !is.null(grob$vp) && !is.null(grob$vp$mask),
    track_scope_grobs(ggplot2::ggplotGrob(after)))
  expect_true(any(vapply(inset_grobs,
    function(grob) inherits(grob$vp$mask, "GridMask"), logical(1))))
})

test_that("scope layers reject unsupported Stats and complete plots explicitly", {
  f <- track_scope_fixture()
  mapping <- ggplot2::aes(chr = Chr, x = Pos, y = Value)
  expect_error(ggideogram(f$model, tracks = list(signal = geom_track(
    data = f$data, mapping = mapping,
    layers = list(ggplot2::stat_summary(geom = "line", fun = mean))))))
  expect_error(ggideogram(f$model, tracks = list(signal = geom_track(
    data = f$data, mapping = mapping,
    layers = list(ggplot2::ggplot(f$data, ggplot2::aes(Pos, Value)) +
      ggplot2::geom_point())))))
  expect_error(ggideogram(f$model, tracks = list(signal = geom_track(
    data = f$data, mapping = mapping,
    geom = ggplot2::geom_point(inherit.aes = FALSE)))))
})

test_that("one track can be resized restyled and switched off without duplicating contents", {
  f <- track_scope_fixture()
  for (orientation in c("horizontal", "vertical", "circular")) {
    p <- ggideogram(f$model, orientation = orientation, show_names = FALSE,
      reverse_chr = f$keys[2], tracks = list(signal = geom_track(side = "left",
        data = f$data, ggplot2::aes(chr = Chr, x = Pos, y = Value),
        layers = list(ggplot2::geom_line(linewidth = .37), ggplot2::geom_point(size = 1.6)),
        label = "Signal", axis = list(breaks = c(1, 6)))))
    q <- p + geom_track(track = "signal", width = 3, gap = .7,
      reverse = TRUE, limits = c(0, 10), label = list(size = 3.1, fontface = "bold", gap = .5),
      axis = list(size = 2.4))
    expect_length(track_scope_indices(q, "GeomLine"), 1)
    expect_length(track_scope_indices(q, "GeomPoint"), 1)
    point <- track_scope_indices(q, "GeomPoint")
    build <- ggplot2::ggplot_build(q)
    actual <- track_scope_xy(build, point)
    projected <- track_scope_expected(q, f$data, "signal")
    built <- build$data[[point]]
    expect_equal(actual$x, projected$x)
    expect_equal(actual$y, projected$y)
    expect_equal(built$position, f$data$Pos)
    expect_equal(q$layers[[point]]$aes_params$size, 1.6)
    title <- Filter(function(layer) !is.null(layer$ideogram_track_title), q$layers)[[1]]
    expect_identical(title$data$label, "Signal")
    expect_equal(title$aes_params$size, 3.1)
    expect_equal(title$aes_params$fontface, "bold")
    axes <- Filter(function(layer) inherits(layer$stat, "StatTrackAxis"), q$layers)
    expect_equal(axes[[1]]$stat_params$axis_spec$breaks, c(1, 6))
    expect_equal(axes[[1]]$stat_params$axis_spec$size, 2.4)
    expect_no_warning(ggplot2::ggplotGrob(q))
    hidden <- q + geom_track(track = "signal", axis = FALSE, label = NULL)
    expect_false(any(vapply(hidden$layers, function(layer)
      inherits(layer$stat, "StatTrackAxis") || !is.null(layer$ideogram_track_title), logical(1))))
    expect_length(track_scope_indices(hidden, "GeomPoint"), 1)
    expect_equal(p$coordinates$layout$tracks$width, 2)
    expect_no_warning(ggplot2::ggplotGrob(hidden))
    edited <- q + geom_track(track = "signal", replace = TRUE,
      layers = list(ggplot2::geom_line(colour = "#4477AA", linewidth = .6),
        ggplot2::geom_point(size = 2.2)))
    expect_length(track_scope_indices(edited, "GeomLine"), 1)
    expect_length(track_scope_indices(edited, "GeomPoint"), 1)
    expect_equal(edited$layers[[track_scope_indices(edited, "GeomPoint")]]$aes_params$size, 2.2)
    expect_equal(edited$coordinates$layout$tracks$width, 3)
    expect_no_warning(ggplot2::ggplotGrob(edited))
  }
})

test_that("track ownership follows chromosome selection and overlay edits", {
  f <- track_scope_fixture()
  base <- ggideogram(f$model, orientation = "horizontal", show_names = FALSE,
    tracks = list(signal = geom_track(data = f$data,
      ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_point(),
      axis = TRUE, label = "Signal")))
  selected <- base + geom_track(track = "signal", chr = f$keys[2])
  point <- selected$layers[[track_scope_indices(selected, "GeomPoint")]]
  expect_identical(unique(point$data$Chr), f$keys[2])
  axes <- Filter(function(layer) inherits(layer$stat, "StatTrackAxis"), selected$layers)
  expect_identical(axes[[1]]$stat_params$chromosomes, f$keys[2])
  all <- selected + geom_track(track = "signal", chr = NULL)
  expect_equal(nrow(all$layers[[track_scope_indices(all, "GeomPoint")]]$data), nrow(f$data))
  body <- base + geom_track(track = "body", side = "overlay", width = .4,
    data = f$data, ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_point())
  shifted <- body + geom_track(track = "body", offset = .2)
  expect_equal(shifted$coordinates$layout$tracks$offset[2], .2)
  expect_no_warning(ggplot2::ggplotGrob(shifted))
})

test_that("rebinding track data respects layer overrides and retrains fixed ranges", {
  f <- track_scope_fixture()
  own <- transform(f$data[1, ], Value = 8, ID = "own")
  p <- ggideogram(f$model, orientation = "horizontal", show_names = FALSE,
    tracks = list(signal = geom_track(data = f$data,
      ggplot2::aes(chr = Chr, x = Pos, y = Value),
      layers = list(ggplot2::geom_point(), ggplot2::geom_point(data = own)))))
  replacement <- transform(f$data, Value = Value * 2, Pos = Pos + 2)
  q <- p + geom_track(track = "signal", data = replacement)
  points <- track_scope_indices(q, "GeomPoint")
  expect_equal(q$layers[[points[1]]]$data, replacement)
  expect_equal(q$layers[[points[2]]]$data, own)
  expect_no_warning(ggplot2::ggplotGrob(q))
  box <- ggideogram(f$model, orientation = "horizontal", show_names = FALSE,
    tracks = list(dist = geom_track(data = f$data,
      ggplot2::aes(chr = Chr, x = Pos, y = Value), geom = ggplot2::geom_boxplot(width = 4))))
  changed <- box + geom_track(track = "dist", limits = c(0, 20), transform = "sqrt")
  expect_equal(changed$coordinates$layout$track_ranges[["track:dist"]]$limits, c(0, 20))
  expect_equal(changed$coordinates$layout$tracks$transform_name, "sqrt")
  expect_no_warning(ggplot2::ggplotGrob(changed))
  automatic <- changed + geom_track(track = "dist", limits = NULL, transform = "identity")
  expect_equal(automatic$coordinates$layout$track_ranges[["track:dist"]]$limits, range(f$data$Value))
  expect_no_warning(ggplot2::ggplotGrob(automatic))
})

test_that("track chromosome selection follows validation of the complete source intervals", {
  f <- track_scope_fixture()
  d <- data.frame(Chr = f$keys, Start = 11, End = c(30, 95),
    Gene = c("a", "b"), Type = "exon", Strand = "+")
  for (component in list(geom_genemodel(labels = FALSE), geom_chr(component = "fill"))) {
    expect_error(ggideogram(f$model, show_names = FALSE, tracks = list(body = geom_track(
      side = "overlay", chr = f$keys[1], data = d,
      ggplot2::aes(chr = Chr, start = Start, end = End, gene = Gene, type = Type, strand = Strand),
      layers = list(component)))), "source chromosome bounds")
  }
})

test_that("native interval contents respect the track's selection and visible window", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  d <- data.frame(Chr = c("A", "A", "B"), Start = c(10, 40, 20), End = c(90, 50, 70),
    ID = c("crossing", "inside", "other"))
  p <- ggideogram(chr_view(k, "A", 20, 80), orientation = "horizontal", show_names = FALSE,
    tracks = list(intervals = geom_track(chr = "A", data = d,
      ggplot2::aes(chr = Chr, start = Start, end = End),
      layers = list(geom_locus(geom = "interval")))))
  interval <- Filter(function(layer) inherits(layer$stat, "StatChrInterval"), p$layers)[[1]]
  expect_identical(interval$data$ID, c("crossing", "inside"))
  expect_equal(interval$data$.source_start, c(10, 40))
  expect_equal(interval$data$.source_end, c(90, 50))
  built <- Filter(function(data) "ideogram_interval" %in% names(data), ggplot2::ggplot_build(p)$data)[[1]]
  expect_equal(built$start, c(21, 40))
  expect_equal(built$end, c(80, 50))
  expect_no_warning(ggplot2::ggplotGrob(p))
  invalid <- d; invalid$End[3] <- 101
  expect_error(ggideogram(k, show_names = FALSE, tracks = list(intervals = geom_track(
    chr = "A", data = invalid, ggplot2::aes(chr = Chr, start = Start, end = End),
    layers = list(geom_locus(geom = "interval"))))), "source chromosome bounds")
})

test_that("overlay tiles and locus layers can be copied and adjusted together", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 50, 80), Value = c(2, 5, 8))
  for (orientation in c("horizontal", "vertical", "circular")) {
    p <- ggideogram(k, orientation = orientation, show_names = FALSE,
      tracks = list(body = geom_track(side = "overlay", width = .6, data = d,
        ggplot2::aes(chr = Chr, x = Pos, y = Value), limits = c(0, 10),
        geom = ggplot2::geom_tile(width = 8, height = 1)),
        markers = geom_track(data = d,
          ggplot2::aes(chr = Chr, position = Pos), layers = list(geom_locus(size = 1.7)))))
    tile <- track_scope_indices(p, "GeomTile")
    point <- track_scope_indices(p, "GeomPoint")
    expect_length(tile, 1)
    expect_length(point, 1)
    expect_s3_class(p$layers[[tile]]$geom, "GeomChrTrackTile")
    expect_equal(p$layers[[point]]$aes_params$size, 1.7)
    expect_no_warning(ggplot2::ggplotGrob(p))
    unclipped <- p + geom_track(track = "body", clip = "off")
    expect_false(inherits(unclipped$layers[[tile]]$geom, "GeomChrTrackTile"))
    expect_equal(unclipped$layers[[point]]$data, d)
    expect_no_warning(ggplot2::ggplotGrob(unclipped))
  }
})
