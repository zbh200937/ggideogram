test_that("cytoband import, centromeres and polygons share exact bp boundaries", {
  ucsc <- data.frame(chr = "A", start = c(0, 1, 40, 60),
    end = c(1, 40, 60, 100), name = letters[1:4],
    stain = c("gneg", "gpos50", "acen", "gneg"))
  bands <- read_cytoband(ucsc)
  k <- cytoband_karyotype(bands)
  expect_equal(k$CE_start, 40)
  expect_equal(k$CE_end, 60)
  p <- ggideogram(k, cytoband = bands, orientation = "horizontal")
  polygons <- cytoband_polygon_data(p$coordinates$layout)
  bp <- unproject_chr_point(p, polygons, ".chr", "x", "y")$.position
  bounds <- t(vapply(split(bp, polygons$.band), range, numeric(2)))
  expect_equal(unname(bounds), cbind(ucsc$start, ucsc$end))
  expect_identical(p$coordinates$layout$data$cytoband$Start, bands$Start)
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("locus prototypes retain or explicitly override native settings", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 20, 80))
  mapping <- ggplot2::aes(chr = Chr, position = Pos)
  native <- ggplot2::geom_point(size = 2, alpha = .8,
    position = ggplot2::position_nudge(x = 1))
  for (geom in list(native, ggplot2::geom_point, "point")) {
    p <- ggideogram(k, show_names = FALSE) + geom_locus(mapping, d,
      geom = geom, size = 8, alpha = .4, position = ggplot2::position_nudge(x = 3))
    layer <- utils::tail(p$layers, 1)[[1]]
    expect_equal(layer$aes_params$size, 8)
    expect_equal(layer$aes_params$alpha, .4)
    expect_equal(layer$position$x, 3)
    expect_equal(utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]$size, rep(8, 3))
  }
  retained <- geom_locus(mapping, d, geom = native)$component
  expect_equal(retained$position$x, 1)
  expect_equal(retained$aes_params$size, 2)
  expect_equal(native$aes_params$size, 2)
  expect_equal(native$position$x, 1)
  expect_error(geom_locus(mapping, d,
    geom = ggplot2::geom_point(stat = "summary", fun = mean)), "require stat = 'identity'")
  expect_error(geom_locus(mapping, d, geom = "point", stat = "summary"), "require stat = 'identity'")
})

test_that("chromosome component styles are applied or diagnosed", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  base <- ggideogram(k, show_names = FALSE)
  p <- base + geom_chr(alpha = .4, linetype = "dashed")
  built <- ggplot2::ggplot_build(p)$data
  expect_equal(unique(built[[3]]$alpha), .4)
  expect_equal(unique(built[[4]]$linetype), "dashed")
  constructor <- ggplot2::ggplot(k) + geom_chr(alpha = .4, linetype = "dashed", names = FALSE)
  expect_equal(ggplot2::ggplot_build(constructor)$data, built[3:4])
  expect_error(base + geom_chr(linewidht = 2), "Unsupported.*linewidht")
  expect_error(base + geom_chr(ggplot2::aes(fill = Chr)), "mapped aesthetics: fill")
  expect_no_warning(ggplot2::ggplotGrob(base + geom_chr(
    component = c("body", "name"), alpha = .4, size = 3)))
})

test_that("sorted interval windows match direct overlapping observations", {
  set.seed(41)
  k <- data.frame(Chr = "A", Start = 0, End = 173)
  d <- data.frame(Chr = "A", Start = sample(1:173, 80, replace = TRUE),
    Score = rep(c(NA, -3, 7, 9), 20), Group = rep(c("a", "b"), 40))
  d$End <- pmin(173, d$Start + sample(0:80, 80, replace = TRUE))
  for (method in c("coverage", "weighted_mean")) for (na.rm in c(FALSE, TRUE)) {
    out <- bin_genome(d, k, window = 17, method = method, group = "Group",
      value = if (method == "weighted_mean") "Score" else NULL, na.rm = na.rm)
    for (i in seq_len(nrow(out))) {
      hit <- d$Group == out$Group[i] & d$Start <= out$End[i] & d$End >= out$Start[i]
      rows <- d[hit, ]
      covered <- lapply(seq_len(nrow(rows)), function(j)
        seq.int(max(out$Start[i], rows$Start[j]), min(out$End[i], rows$End[j])))
      expected <- if (method == "coverage") length(unique(unlist(covered))) / out$Width[i] else
        stats::weighted.mean(rows$Score, lengths(covered), na.rm = na.rm)
      if (is.nan(expected)) expected <- NA_real_
      expect_equal(out$Value[i], expected)
      expect_equal(out$N[i], nrow(rows))
      expect_equal(out$N_valid[i], if (method == "coverage") nrow(rows) else sum(!is.na(rows$Score)))
    }
  }
})

test_that("GFF3 example models retain shared exon relationships", {
  e <- new.env()
  sys.source(system.file("examples", "gene-model-data.R", package = "ggideogram"), e)
  d <- data.frame(Chr = "A", Type = c("mRNA", "mRNA", "exon"),
    ID = c("t1", "t2", "e1"), Parent = c("g1", "g1", "t1,t2"),
    Start = c(1, 1, 20), End = c(100, 100, 40))
  models <- e$gene_model_data(d)
  shared <- models[models$Type == "exon", ]
  expect_equal(shared$Transcript, c("t1", "t2"))
  expect_equal(shared$Parent, rep("t1,t2", 2))
  expect_equal(shared$Gene, rep("g1", 2))
  expect_equal(shared$Start, rep(20, 2))
})

test_that("text layout rejects overridden statistics consistently", {
  for (geom in list("text", ggplot2::geom_text, ggplot2::geom_text())) {
    expect_error(geom_locus(geom = geom, position = "spread", stat = "summary"),
      "requires stat = 'identity'")
    text <- geom_locus(geom = geom, position = "spread", size = 5)$component
    expect_equal(text$params$size, 5)
  }
})
