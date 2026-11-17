test_that("log10 and non-zero limits work for geom_col and geom_area across value scales", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 120))
  d <- data.frame(
    Chr = c("A", "A", "B", "B"),
    Pos = c(20, 70, 30, 90),
    Val = c(10, 100, 50, 1000)
  )
  for (vscale in c("per_track", "per_chr", "global")) {
    p_log <- ggideogram(k, orientation = "horizontal", tracks = list(
      col_log = geom_track(
        data = d, ggplot2::aes(chr = Chr, x = Pos, y = Val),
        geom = ggplot2::geom_col(width = 10),
        transform = "log10", value_scale = vscale, axis = TRUE
      ),
      area_log = geom_track(
        data = d, ggplot2::aes(chr = Chr, x = Pos, y = Val),
        geom = ggplot2::geom_area(),
        transform = "log10", value_scale = vscale, axis = TRUE
      )
    ))
    p_lim <- ggideogram(k, orientation = "horizontal", tracks = list(
      col_lim = geom_track(
        data = d, ggplot2::aes(chr = Chr, x = Pos, y = Val),
        geom = ggplot2::geom_col(width = 10),
        limits = c(5, 1500), value_scale = vscale, axis = TRUE
      )
    ))
    expect_no_warning(ggplot2::ggplotGrob(p_log))
    expect_no_warning(ggplot2::ggplotGrob(p_lim))
  }
})

test_that("geom_hline, geom_vline, geom_rect, geom_linerange and geom_errorbar work natively in geom_track", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  d <- data.frame(
    Chr = c("A", "A", "B"),
    Pos = c(25, 75, 50),
    Val = c(2, 8, 5),
    Lo = c(1, 6, 3),
    Hi = c(3, 9, 7),
    Left = c(10, 60, 35),
    Right = c(40, 90, 65)
  )
  for (orientation in c("horizontal", "vertical", "circular")) {
    expect_no_warning({
      p <- ggideogram(k, orientation = orientation, show_names = FALSE) +
        geom_track(
          track = "qtl", data = d,
          ggplot2::aes(chr = Chr, x = Pos, y = Val),
          layers = list(
            ggplot2::geom_hline(yintercept = 5, linetype = "dashed"),
            ggplot2::geom_vline(xintercept = 50, linetype = "dotted"),
            ggplot2::geom_rect(ggplot2::aes(xmin = Left, xmax = Right, ymin = Lo, ymax = Hi), alpha = 0.2),
            ggplot2::geom_linerange(ggplot2::aes(ymin = Lo, ymax = Hi)),
            ggplot2::geom_errorbar(ggplot2::aes(ymin = Lo, ymax = Hi), width = 4),
            ggplot2::geom_point()
          ),
          axis = TRUE
        )
      ggplot2::ggplotGrob(p)
    })
  }

  # Boundary-crossing GeomRect clips cleanly inside chr_view
  v <- chr_view(as_ideogram_data(k), "A", 20, 80)
  p_clip <- ggideogram(v, show_names = FALSE) +
    geom_track(
      track = "r", data = d,
      ggplot2::aes(chr = Chr, xmin = Left, xmax = Right, ymin = Lo, ymax = Hi),
      geom = ggplot2::geom_rect()
    )
  built <- ggplot2::ggplot_build(p_clip)
  rect_df <- utils::tail(built$data, 1)[[1]]
  expect_equal(nrow(rect_df), 2L)
  expect_no_warning(ggplot2::ggplotGrob(p_clip))
})

test_that("project_chr_track handles 0-row and NA values and constant track axes build", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d <- data.frame(Chr = "A", Pos = c(20, 50, 80), Val = c(5, NA, 5))
  p <- ggideogram(k, show_names = FALSE) +
    geom_track(
      track = "const", data = d,
      ggplot2::aes(chr = Chr, x = Pos, y = Val),
      geom = ggplot2::geom_point(na.rm = TRUE), axis = TRUE
    )
  expect_no_warning(ggplot2::ggplotGrob(p))
  empty_proj <- project_chr_track(p, d[0, ], "Chr", "Pos", "Val", "const")
  expect_equal(nrow(empty_proj), 0L)
  expect_true(all(c(".x", ".y") %in% names(empty_proj)))
  na_proj <- project_chr_track(p, d, "Chr", "Pos", "Val", "const")
  expect_true(is.na(na_proj$.x[2]) && is.na(na_proj$.y[2]))
  expect_true(is.finite(na_proj$.x[1]) && is.finite(na_proj$.y[1]))
})

test_that("track rebuild preserves skeleton omissions, replace slots, global ranges and scope inheritance flags", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  d1 <- data.frame(Chr = "A", Pos = c(20, 80), Val = c(1, 5))
  d2 <- data.frame(Chr = "B", Pos = c(20, 80), Val = c(1, 50))

  # Name-only skeleton stays body-free after adding a track
  p_name <- ggplot2::ggplot(k) + geom_chr(component = "name") +
    geom_track(track = "s", data = d1, ggplot2::aes(chr = Chr, x = Pos, y = Val),
      geom = ggplot2::geom_point())
  expect_false(any(vapply(p_name$layers, function(l) inherits(l$geom, "GeomChromosome"), logical(1))))

  # Adding a new global track updates existing global statistical ranges
  p_glob <- ggideogram(k, show_names = FALSE, tracks = list(
    g1 = geom_track(value_scale = "global", data = d1,
      ggplot2::aes(chr = Chr, x = Pos, y = Val), geom = ggplot2::geom_area())
  )) +
    geom_track(track = "g2", value_scale = "global", data = d2,
      ggplot2::aes(chr = Chr, x = Pos, y = Val), geom = ggplot2::geom_point())
  expect_equal(p_glob$coordinates$layout$track_ranges[["global"]]$limits, c(0, 50))

  # replace = TRUE changing geom class and expanding content slots keeps order before trailing layers
  p_rep <- ggideogram(k, show_names = FALSE, tracks = list(
    s = geom_track(data = d1, ggplot2::aes(chr = Chr, x = Pos, y = Val),
      geom = ggplot2::geom_point(shape = 21, fill = "red"))
  )) +
    ggplot2::annotate("point", x = 1, y = 1, colour = "black")
  p_rep2 <- p_rep + geom_track(track = "s", replace = TRUE,
    layers = list(ggplot2::geom_line(), ggplot2::geom_point(colour = "blue")))
  geoms <- vapply(p_rep2$layers, function(l) class(l$geom)[1], character(1))
  line_i <- which(geoms == "GeomLine")
  pt_i <- which(geoms == "GeomPoint")
  expect_true(line_i < pt_i[1] && pt_i[1] < pt_i[2])
  expect_no_warning(ggplot2::ggplotGrob(p_rep2))

  # inherit.aes = FALSE inside geom_track prevents inheriting bad parent mappings
  pts <- data.frame(Chr = "A", Pos = 40)
  p_inh <- ggideogram(k, show_names = FALSE) +
    geom_track(
      track = "s", data = d1,
      ggplot2::aes(chr = BogusColumn, x = BogusX, y = BogusY),
      layers = list(
        geom_locus(ggplot2::aes(chr = Chr, position = Pos), data = pts, inherit.aes = FALSE,
          key_glyph = "rect")
      )
    )
  expect_no_warning(ggplot2::ggplotGrob(p_inh))

  # geom_jitter keeps PositionJitter and reverse_chr validates against source karyotype
  p_jit <- ggideogram(chr_view(as_ideogram_data(k), "A", 10, 90), reverse_chr = "B", show_names = FALSE) +
    geom_locus(ggplot2::aes(chr = Chr, position = Pos), data = d1,
      geom = ggplot2::geom_jitter(width = 0.1, height = 0))
  expect_s3_class(utils::tail(p_jit$layers, 1)[[1]]$position, "PositionChrLocus")
  expect_no_warning(ggplot2::ggplotGrob(p_jit))
})

test_that("four-genome horizontal synteny stack shares compact tracks, diamond icons and single bottom axis", {
  k4 <- data.frame(
    Genome = paste0("G", 1:4),
    Chr = "1",
    Start = 0,
    End = c(100, 95, 110, 105),
    Label = paste0("Genome ", 1:4)
  )
  sem4 <- as_ideogram_data(k4, ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome, label = Label))
  keys4 <- sem4$karyotype$.chr
  edges4 <- data.frame(
    a = keys4[1:3], b = keys4[2:4],
    s1 = c(10, 15, 20), e1 = c(40, 45, 50),
    s2 = c(12, 18, 22), e2 = c(42, 48, 52)
  )
  loci4 <- data.frame(
    Key = keys4, Pos = c(25, 30, 35, 40),
    Name = c("top_gene", "mid_gene", "", "bot_gene")
  )
  p4 <- ggideogram(sem4, orientation = "horizontal", ncol = 4, show_names = TRUE) +
    geom_chrlink(ggplot2::aes(chr1 = a, start1 = s1, end1 = e1, chr2 = b, start2 = s2, end2 = e2), edges4) +
    geom_chr(component = "axis") +
    geom_locus(ggplot2::aes(chr = Key, position = Pos), data = data.frame(Key = keys4[1], Pos = 50),
      side = "left", gap = 1.5) +
    geom_locus(ggplot2::aes(chr = Key, position = Pos, label = Name), data = loci4,
      geom = "text", position = "spread")

  expect_true(isTRUE(p4$coordinates$layout$microsynteny))
  axes <- Filter(function(l) !is.null(l$ideogram_axis_spec), p4$layers)
  expect_equal(nrow(axes[[1]]$data), 1L)
  icons <- Filter(function(l) inherits(l$geom, "GeomPoint") && identical(l$aes_params$shape, 23), p4$layers)
  expect_length(icons, 1L)
  expect_equal(nrow(icons[[1]]$data), 1L) # empty string on keys4[3] skipped cleanly
  # Untracked marker on keys4[1] left does not inflate annotation_extent on keys4[2..4]
  expect_equal(annotation_extent(p4$coordinates$layout, keys4[2], "left"), 0)
  expect_no_warning(ggplot2::ggplotGrob(p4))
})

test_that("circular bend = 0, concentric local curves, fractional insets and 0-row semantic recycling work", {
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = 100)
  pc <- ggideogram(k, orientation = "circular", show_names = FALSE) +
    geom_chrlink(ggplot2::aes(chr1 = a, position1 = p1, chr2 = b, position2 = p2),
      data.frame(a = "A", b = "B", p1 = 20, p2 = 70), bend = 0)
  expect_equal(nrow(utils::tail(pc$layers, 1)[[1]]$data), 1L)

  # Concentric intervals with identical midpoint curve cleanly without collapsing
  p_conc <- ggideogram(k[1, , drop = FALSE], orientation = "horizontal", show_names = FALSE) +
    geom_chrlink(ggplot2::aes(chr1 = a, start1 = s1, end1 = e1, chr2 = b, start2 = s2, end2 = e2),
      data.frame(a = "A", s1 = 41, e1 = 60, b = "A", s2 = 31, e2 = 70))
  poly_y <- utils::tail(p_conc$layers, 1)[[1]]$data$.pair_y
  expect_gt(diff(range(poly_y)), 0.1)

  # Fractional chr_view interval inset and one-to-one chr_data validation
  vf <- chr_view(as_ideogram_data(k), "A", 10.25, 60.75)
  child <- ggplot2::ggplot(data.frame(x = 1:2, y = 1:2), ggplot2::aes(x, y)) + ggplot2::geom_point()
  p_ins <- ggideogram(vf, show_names = FALSE) +
    geom_locus_inset(child, ggplot2::aes(chr = Chr, start = Start, end = End),
      data = data.frame(Chr = "A", Start = 5, End = 80),
      width = grid::unit(8, "mm"), height = grid::unit(8, "mm"))
  expect_no_warning(ggplot2::ggplotGrob(p_ins))
  expect_error(
    ggideogram(vf) + chr_data(data.frame(Chr = "A", Group = "x"), relationship = "one-to-one"),
    "missing"
  )
  expect_equal(length(recycle_semantic("const", 0L, "test")), 0L)
})

test_that("PositionChrLocus dodges links with points, 360-deg circular repel fits seam, and track titles clear gene labels", {
  k <- data.frame(Chr = "A", Start = 0, End = 100)
  d_dodge <- data.frame(Chr = "A", Pos = c(50, 50), Grp = c("g1", "g2"))
  pd <- ggideogram(k, orientation = "horizontal", show_names = FALSE) +
    geom_locus(ggplot2::aes(chr = Chr, position = Pos, colour = Grp, group = Grp), d_dodge,
      geom = "link", position = ggplot2::position_dodge(width = 6)) +
    geom_locus(ggplot2::aes(chr = Chr, position = Pos, colour = Grp, group = Grp), d_dodge,
      geom = "point", position = ggplot2::position_dodge(width = 6))
  bd <- ggplot2::ggplot_build(pd)
  lnk <- bd$data[[length(bd$data) - 1L]]
  pts <- bd$data[[length(bd$data)]]
  expect_equal(lnk$xend, pts$x, tolerance = 1e-8)
  expect_equal(lnk$yend, pts$y, tolerance = 1e-8)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_ring <- ggideogram(k, orientation = "circular", opening_angle = 0, radius = 25, show_names = FALSE) +
      geom_locus(ggplot2::aes(chr = Chr, position = Pos, label = Lbl),
        data.frame(Chr = "A", Pos = c(1, 50, 99), Lbl = c("SeamA", "Mid", "SeamB")),
        geom = "text", position = "repel")
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_warning(grid::grid.draw(ggplot2::ggplotGrob(p_ring)))
  }

  # Track title on a gene-model track clears gene lane labels while chromosome name stays near body
  genes <- data.frame(
    Chr = "A", Start = c(10, 55), End = c(35, 85),
    Type = c("CDS", "CDS"), Strand = c("+", "-"),
    Gene = c("LONG_GENE_ALPHA", "LONG_GENE_BETA")
  )
  pg <- ggideogram(k, orientation = "horizontal") +
    geom_track(
      track = "genes", side = "right", width = 4, label = "Annotations",
      data = genes, ggplot2::aes(chr = Chr, start = Start, end = End, type = Type, strand = Strand, gene = Gene),
      layers = list(geom_genemodel(mode = "gene", labels = TRUE))
    )
  title_layer <- Filter(function(l) !is.null(l$ideogram_track_title), pg$layers)[[1]]
  expect_gt(grid::convertUnit(title_layer$geom_params$text_clearance, "mm", valueOnly = TRUE), 10)
  expect_no_warning(ggplot2::ggplotGrob(pg))
})

test_that("chromosome axis handles na.translate, rotated labels and telocentric centromere waist clamping", {
  telo_kar <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80), CE_start = c(0, 20), CE_end = c(6, 30))
  prof <- chromosome_axis_profile(as_ideogram_data(telo_kar)$karyotype[1, , drop = FALSE], 20, 4, 32)
  expect_true(all(prof$cross >= 0 & prof$cross <= 1))
  p_ax <- ggplot2::ggplot(data.frame(Chr = c("A", NA), Val = c(4, 2)), ggplot2::aes(Chr, Val)) +
    ggplot2::geom_col() +
    scale_x_chromosome(telo_kar, limits = NULL, na.translate = TRUE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
  expect_no_warning(ggplot2::ggplotGrob(p_ax))
})

test_that("binning, feature readers, widget URLs and composite showcase script pass end to end", {
  kar <- data.frame(Chr = "A", Start = 0, End = 100)
  d_nan <- data.frame(Chr = c("A", "A", "A"), Start = c(10, 20, 30), End = c(15, 25, 35), Grp = c(1, NA_real_, NaN))
  b_nan <- bin_genome(d_nan, kar, window = 100, method = "count", group = "Grp")
  expect_equal(sum(b_nan$N), 3L)
  expect_error(bin_genome(d_nan, kar, window = 100, value = "Start", FUN = range), "single numeric")

  cb_file <- tempfile(fileext = ".txt")
  on.exit(unlink(cb_file), add = TRUE)
  writeLines(c("F\t0\t100\tp1\tgneg", "T\t0\t80\tp1\tgpos50", "NA\t0\t60\tp1\tacen"), cb_file)
  cb <- read_cytoband(cb_file)
  expect_identical(cb$Chr, c("F", "T", "NA"))

  gtf_file <- tempfile(fileext = ".gtf")
  on.exit(unlink(gtf_file), add = TRUE)
  writeLines(
    "chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id \"g1\"; gene_name \" NAC \\\"alpha\\\" 1 \"; transcript_id \"t1\";",
    gtf_file
  )
  gtf_df <- read_chr_features(gtf_file)
  expect_identical(gtf_df$gene_name, " NAC \"alpha\" 1 ")

  if (requireNamespace("ggiraph", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    pw <- ggideogram(kar) + geom_locus(
      ggplot2::aes(chr = Chr, position = Start, data_id = ID),
      data.frame(ID = "x1", Chr = "A", Start = 20),
      geom = ggiraph::geom_point_interactive
    )
    w_ok <- as_ideogram_widget(
      pw,
      data.frame(ID = "x1", Chr = "A", Start = 20, End = 30, Cat = "", URL = "HTTPS://example.org/gene"),
      category = "Cat", url = "URL"
    )
    expect_identical(w_ok$jsHooks$render[[1]]$data$records[[1]]$category, "\u672a\u5206\u7c7b")
    expect_error(
      as_ideogram_widget(
        pw,
        data.frame(ID = "x1", Chr = "A", Start = 20, End = 30, URL = "https://example.org/a b"),
        url = "URL"
      ),
      "HTTP"
    )
  }
})
