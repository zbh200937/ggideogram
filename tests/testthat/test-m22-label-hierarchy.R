test_that("M22 decouples per-track title clearance, auto-rotates discrete track titles, and keeps Chr names snug", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  kar <- data.frame(Chr = "Chr2", Start = 0, End = 20e6)
  sig <- data.frame(Chr = "Chr2", Pos = seq(1e6, 19e6, by = 3e6), Val = c(2, 8, 35, 12, 5, 9, 3))
  lines <- paste0("MAGIC.", c(294, 658, 482, 310, 51, 290), " (55 cm)")
  bins <- data.frame(
    Chr = "Chr2",
    Start = rep(c(1, 10000001), each = length(lines)),
    End = rep(c(10000000, 20000000), each = length(lines)),
    Line = factor(rep(lines, 2), levels = lines),
    State = rep(c("Col", "Ler"), length(lines))
  )

  make_binmap_plot <- function(ril_label) {
    ggideogram(kar, orientation = "horizontal", tracks = list(
      gwas = geom_track(
        side = "left", width = 2.6, limits = c(0, 40),
        data = sig, mapping = ggplot2::aes(chr = Chr, x = Pos, y = Val),
        geom = ggplot2::geom_line(),
        label = "Height GWAS",
        axis = list(breaks = c(0, 20, 40))
      ),
      ril = geom_track(
        side = "right", width = 6.8,
        data = bins,
        mapping = ggplot2::aes(chr = Chr, xmin = Start - 1, xmax = End, y = Line, fill = State),
        geom = ggplot2::geom_rect(),
        label = ril_label,
        axis = list(position = "start", size = 2.2)
      )
    ))
  }

  p_auto <- make_binmap_plot("Founders & RILs")
  p_horiz <- make_binmap_plot(list(text = "Founders & RILs", angle = 0))

  titles_auto <- Filter(function(l) !is.null(l$ideogram_track_title), p_auto$layers)
  expect_length(titles_auto, 2L)
  expect_equal(titles_auto[[1]]$aes_params$angle, 0)
  expect_equal(titles_auto[[2]]$aes_params$angle, 90)
  expect_equal(
    grid::convertUnit(titles_auto[[1]]$geom_params$text_clearance, "mm", valueOnly = TRUE),
    0
  )
  expect_gt(
    grid::convertUnit(titles_auto[[2]]$geom_params$text_clearance, "mm", valueOnly = TRUE),
    15
  )

  # Auto 90-deg title on wide discrete track saves substantial left margin vs angle = 0
  # without falsely inflating top/bottom plot margins.
  margins_auto <- grid::convertUnit(p_auto$theme$plot.margin, "mm", valueOnly = TRUE)
  left_auto <- margins_auto[4]
  left_horiz <- grid::convertUnit(p_horiz$theme$plot.margin[4], "mm", valueOnly = TRUE)
  expect_lt(left_auto + 12, left_horiz)
  expect_lte(margins_auto[1], 5)
  expect_lte(margins_auto[3], 5)

  # Rendered Chr2 name stays snug beside chromosome capsule rather than jumping 25+ mm left
  grid::pushViewport(grid::viewport(width = grid::unit(180, "mm"), height = grid::unit(120, "mm")))
  b <- ggplot2::ggplot_build(p_auto)
  name_i <- which(vapply(b$plot$layers, function(l) inherits(l$geom, "GeomChromosomeName"), logical(1)))
  ril_title_i <- which(vapply(b$plot$layers, function(l)
    identical(l$ideogram_track_title$track, "ril"), logical(1)))

  draw_layer <- function(idx) {
    layer <- b$plot$layers[[idx]]
    do.call(layer$geom$draw_panel, c(
      list(data = b$data[[idx]], panel_params = b$layout$panel_params[[1]], coord = b$plot$coordinates),
      layer$computed_geom_params[intersect(names(layer$computed_geom_params), layer$geom$parameters())]
    ))
  }
  name_grob <- grid::makeContent(draw_layer(name_i))$children[[1]]
  ril_grob <- draw_layer(ril_title_i)
  chr_x_mm <- grid::convertX(name_grob$x[1], "mm", valueOnly = TRUE)
  ril_x_mm <- grid::convertX(ril_grob$x[1], "mm", valueOnly = TRUE)
  expect_gt(chr_x_mm, ril_x_mm + 15)
  grid::popViewport()

  # Rebuild updates angle and position cleanly
  p_top <- p_auto + geom_track(track = "gwas", label = list(position = "top"))
  gwas_top <- Filter(function(l) identical(l$ideogram_track_title$track, "gwas"), p_top$layers)[[1]]
  expect_identical(gwas_top$ideogram_track_title$position, "top")
  expect_equal(gwas_top$geom_params$title_hjust, 0)
  expect_equal(gwas_top$data$nx, 0)
  expect_equal(gwas_top$data$ny, 1)
  expect_no_warning(ggplot2::ggplotGrob(p_top))

  expect_error(geom_track(label = list(text = "Bad", position = "bottom")), "position")
  expect_error(geom_track(label = list(text = "Bad", angle = NA_real_)), "angle")
})

test_that("M22 geom_genemodel packs non-overlapping genes into compact lanes with local labels", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  kar <- chr_view(data.frame(Chr = "Chr2", Start = 0, End = 20e6), "Chr2", 11196000, 11231000)
  starts <- seq(11196500, by = 3200, length.out = 10)
  genes <- data.frame(
    Chr = "Chr2",
    Start = starts,
    End = starts + 1800,
    Type = "CDS",
    Strand = rep(c("+", "-"), 5),
    Gene = sprintf("AT2G263%d0", 0:9)
  )

  p_packed <- ggideogram(kar, orientation = "horizontal", tracks = list(
    models = geom_track(
      side = "right", width = 5.4,
      data = genes,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, type = Type, strand = Strand, gene = Gene),
      label = "TAIR10 genes",
      layers = list(geom_genemodel(mode = "gene", label_fontface = "italic"))
    )
  ))
  p_unpacked <- ggideogram(kar, orientation = "horizontal", tracks = list(
    models = geom_track(
      side = "right", width = 5.4,
      data = genes,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, type = Type, strand = Strand, gene = Gene),
      label = "TAIR10 genes",
      layers = list(geom_genemodel(mode = "gene", pack = FALSE, label_position = "start"))
    )
  ))

  seg_packed <- Filter(function(l) inherits(l$geom, "GeomSegment"), p_packed$layers)[[1]]$data
  seg_unpacked <- Filter(function(l) inherits(l$geom, "GeomSegment"), p_unpacked$layers)[[1]]$data
  lanes_packed <- length(unique(round(seg_packed$.gy, 6)))
  lanes_unpacked <- length(unique(round(seg_unpacked$.gy, 6)))
  expect_equal(lanes_unpacked, 10L)
  expect_lte(lanes_packed, 3L)

  # Local "above" labels stay inside [x_min, x_max] and require zero left gutter clearance
  lbl_layer <- Filter(function(l) isTRUE(l$ideogram_gene_label), p_packed$layers)[[1]]
  expect_false(isTRUE(lbl_layer$ideogram_gene_label_start))
  title_packed <- Filter(function(l) !is.null(l$ideogram_track_title), p_packed$layers)[[1]]
  expect_equal(
    grid::convertUnit(title_packed$geom_params$text_clearance, "mm", valueOnly = TRUE),
    0
  )
  left_packed <- grid::convertUnit(p_packed$theme$plot.margin[4], "mm", valueOnly = TRUE)
  left_unpacked <- grid::convertUnit(p_unpacked$theme$plot.margin[4], "mm", valueOnly = TRUE)
  expect_lt(left_packed + 12, left_unpacked)

  # "beside" and vertical orientations build cleanly
  p_beside <- ggideogram(kar, orientation = "horizontal", tracks = list(
    models = geom_track(
      side = "right", width = 5.0,
      data = genes,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, type = Type, strand = Strand, gene = Gene),
      layers = list(geom_genemodel(mode = "gene", label_position = "beside"))
    )
  ))
  p_vert <- ggideogram(kar, orientation = "vertical", tracks = list(
    models = geom_track(
      side = "right", width = 5.0,
      data = genes,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, type = Type, strand = Strand, gene = Gene),
      layers = list(geom_genemodel(mode = "gene", label_position = "above"))
    )
  ))
  expect_no_warning(ggplot2::ggplotGrob(p_packed))
  expect_no_warning(ggplot2::ggplotGrob(p_beside))
  expect_no_warning(ggplot2::ggplotGrob(p_vert))

  expect_error(
    geom_genemodel(mode = "gene", pack = TRUE, lane_order = genes$Gene),
    "lane_order"
  )
})

test_that("M22 packed gene models and rotated/top track titles open zero graphics devices", {
  while (grDevices::dev.cur() > 1L) grDevices::dev.off()
  old_opt <- options(device = function(...) stop("Unexpected graphics device opened"))
  on.exit(options(old_opt), add = TRUE)

  kar <- data.frame(Chr = "Chr2", Start = 0, End = 10000)
  genes <- data.frame(
    Chr = "Chr2", Start = c(500, 4200, 7600), End = c(1800, 5600, 9100),
    Type = "exon", Strand = c("+", "-", "+"), Gene = c("G1", "G2", "G3")
  )
  p <- ggideogram(kar, orientation = "horizontal", tracks = list(
    top_hdr = geom_track(
      side = "left", width = 2,
      data = data.frame(Chr = "Chr2", Pos = c(1000, 8000), Val = c(1, 5)),
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Val),
      geom = ggplot2::geom_line(),
      label = list(text = "Header", position = "top")
    ),
    genes = geom_track(
      side = "right", width = 5,
      data = genes,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, type = Type, strand = Strand, gene = Gene),
      label = list(text = "Genes", angle = 90),
      layers = list(geom_genemodel(mode = "gene"))
    )
  ))
  built <- ggplot2::ggplot_build(p)
  expect_identical(unname(grDevices::dev.cur()), 1L)
  expect_s3_class(built, "ggplot_built")
})
