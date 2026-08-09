component_kar <- data.frame(
  Chr = c("A", "B"), Start = 0, End = c(100, 80),
  CE_start = c(40, 30), CE_end = c(55, 45)
)

test_that("chromosome metadata attaches by an explicit key", {
  metadata <- data.frame(
    Chromosome = c("B", "A"), Species = c("sp2", "sp1"), Score = c(2, 1)
  )
  base <- ggideogram(component_kar, show_names = FALSE)
  p <- base + chr_data(
    metadata, by = c(chr = "Chromosome"), name = "metadata")

  expect_false("Species" %in% names(base$coordinates$layout$chrom))
  expect_equal(p$coordinates$layout$chrom$Species, c("sp1", "sp2"))
  expect_equal(p$coordinates$layout$data$karyotype$Score, c(1, 2))
  attached <- get_ideogram_attachment(p, "metadata", type = "chr")
  expect_equal(attached$.chr, c("B", "A"))
  expect_equal(attached$Species, metadata$Species)
})

test_that("chromosome attachments make relationships and collisions explicit", {
  duplicated <- data.frame(Chr = c("A", "A"), Value = 1:2)
  p <- ggideogram(component_kar, show_names = FALSE)
  expect_error(
    p + chr_data(duplicated, name = "duplicate"),
    "duplicate key"
  )

  many <- p + chr_data(
    duplicated, relationship = "many-to-many", name = "many")
  expect_equal(nrow(get_ideogram_attachment(many, "many", "chr")), 2L)
  expect_false("Value" %in% names(many$coordinates$layout$chrom))

  expect_error(
    p + chr_data(data.frame(Chr = "A", Value = 1),
                 relationship = "one-to-one", name = "partial"),
    "missing"
  )
  expect_error(
    p + chr_data(data.frame(Chr = "Q", Value = 1), name = "unknown"),
    "unknown chromosome"
  )
  expect_error(
    p + chr_data(data.frame(Chr = "A", .row = 1), name = "collision"),
    "overwrite"
  )
})

test_that("locus attachments validate keys and genomic bounds", {
  loci <- data.frame(
    Chromosome = c("A", "B"), Pos = c(10, 20), Gene = c("g1", "g2")
  )
  p <- ggideogram(component_kar, show_names = FALSE) +
    locus_data(
      loci, by = c(chr = "Chromosome", position = "Pos"),
      name = "genes")
  attached <- get_ideogram_attachment(p, "genes", "locus")
  expect_equal(attached$.chr, loci$Chromosome)
  expect_equal(attached$.position, loci$Pos)

  duplicated <- rbind(loci[1, ], loci[1, ])
  expect_error(
    ggideogram(component_kar) + locus_data(
      duplicated, by = c(chr = "Chromosome", position = "Pos"),
      name = "duplicate"),
    "duplicate key"
  )
  many <- ggideogram(component_kar) + locus_data(
    duplicated, by = c(chr = "Chromosome", position = "Pos"),
    relationship = "many-to-many", name = "many")
  expect_equal(nrow(get_ideogram_attachment(many, "many", "locus")), 2L)
  expect_error(
    ggideogram(component_kar) + locus_data(
      data.frame(Chr = "A", Pos = 101), name = "outside"),
    "outside"
  )
})

make_component_child <- function(with_legend = FALSE) {
  p <- ggplot2::ggplot(
    data.frame(x = 1:3, y = c(1, 3, 2), group = c("a", "b", "a")),
    ggplot2::aes(x, y)
  ) + ggplot2::geom_line() + ggplot2::theme_minimal(base_size = 7)
  if (with_legend) {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = group)) +
      ggplot2::labs(colour = "Child group")
  }
  p
}

test_that("inset collisions are explicit instead of silently covering bodies", {
  child_plot <- make_component_child(with_legend = TRUE)
  child_grob <- grid::rectGrob(gp = grid::gpar(fill = "gold"))
  insets <- data.frame(Chr = c("A", "B"), Pos = c(40, 40))
  insets$Plot <- I(list(child_plot, child_grob))
  crowded <- ggideogram(component_kar, show_names = FALSE) +
    geom_chr_inset(
      data = insets,
      ggplot2::aes(chr = Chr, position = Pos, plot = Plot),
      width = grid::unit(0.25, "npc"),
      height = grid::unit(0.30, "npc"),
      gap = 0.3
    )
  expect_error(ggplot2::ggplotGrob(crowded), "Inset viewport collision")

  tracks <- track_layout(insets = track("right", width = 3, gap = 0.3))
  p <- ggideogram(
    component_kar, tracks = tracks, show_names = FALSE) +
    geom_chr_inset(
      data = insets,
      ggplot2::aes(chr = Chr, position = Pos, plot = Plot),
      width = grid::unit(0.15, "npc"),
      height = grid::unit(0.25, "npc"),
      track = "insets"
    )
  layer <- p$layers[[length(p$layers)]]
  built <- ggplot2::ggplot_build(p)

  expect_s3_class(layer$geom, "GeomChrInset")
  expect_equal(layer$geom_params$width, grid::unit(0.15, "npc"))
  expect_equal(layer$geom_params$height, grid::unit(0.25, "npc"))
  expect_equal(built$data[[length(built$data)]]$ideogram_position,
               c(40, 40))
  expect_length(built$plot$guides$aesthetics, 0L)
  expect_no_error(ggplot2::ggplotGrob(p))

  coordinates <- p$coordinates$transform(
    utils::tail(built$data, 1)[[1]], built$layout$panel_params[[1]])
  justification <- ggideogram:::inset_justification(
    p$coordinates, coordinates, "beside", NULL, NULL)
  expect_equal(justification$hjust, c(0, 0))
  expect_equal(justification$vjust, c(0.5, 0.5))
})

test_that("insets support intervals, centre placement and declared tracks", {
  child <- make_component_child()
  interval <- data.frame(Chr = "A", Start = 20, End = 60)
  interval$Plot <- I(list(child))
  centred <- ggideogram(component_kar[1, ], show_names = FALSE) +
    geom_chr_inset(
      data = interval,
      ggplot2::aes(chr = Chr, start = Start, end = End, plot = Plot),
      width = grid::unit(15, "mm"), height = grid::unit(15, "mm"),
      placement = "center"
    )
  centred_data <- utils::tail(
    ggplot2::ggplot_build(centred)$data, 1)[[1]]
  expect_equal(centred_data$ideogram_position, 40)
  expect_false("ideogram_side" %in% names(centred_data))

  tracks <- track_layout(insets = track("right", width = 1, gap = 0.2))
  tracked <- ggideogram(
    component_kar[1, ], tracks = tracks, show_names = FALSE) +
    geom_chr_inset(
      data = interval,
      ggplot2::aes(chr = Chr, start = Start, end = End, plot = Plot),
      width = grid::unit(0.15, "npc"),
      height = grid::unit(0.15, "npc"),
      track = "insets"
    )
  built <- ggplot2::ggplot_build(tracked)
  inset_data <- utils::tail(built$data, 1)[[1]]
  xy <- tracked$coordinates$transform(
    inset_data, built$layout$panel_params[[1]])
  axis <- project_chr_point(
    tracked$coordinates$layout,
    data.frame(Chr = "A", Pos = 40), "Chr", "Pos")
  expect_gt(xy$x, 0.5)
  expect_equal(axis$.position, 40)
  expect_no_error(ggplot2::ggplotGrob(tracked))
})

test_that("inset contracts reject ambiguous units and payloads", {
  child <- make_component_child()
  data <- data.frame(Chr = "A", Pos = 40)
  data$Plot <- I(list(child))
  mapping <- ggplot2::aes(chr = Chr, position = Pos, plot = Plot)

  expect_error(
    geom_chr_inset(data = data, mapping,
                   width = 0.2, height = grid::unit(1, "cm")),
    "grid::unit"
  )
  bad <- data.frame(Chr = "A", Pos = 40)
  bad$Plot <- I(list("not a plot"))
  expect_error(
    geom_chr_inset(
      data = bad, mapping,
      width = grid::unit(1, "cm"), height = grid::unit(1, "cm")),
    "list-column"
  )
  expect_error(
    geom_chr_inset(
      data = data,
      ggplot2::aes(chr = Chr, position = Pos, start = Pos, end = Pos,
                   plot = Plot),
      width = grid::unit(1, "cm"), height = grid::unit(1, "cm")),
    "either"
  )
})

test_that("ideogram grob extraction is the standard ggplotGrob boundary", {
  p <- ggideogram(component_kar)
  grob <- as_ideogram_grob(p)
  direct <- ggplot2::ggplotGrob(p)
  expect_s3_class(grob, "gtable")
  expect_equal(grob$layout, direct$layout)
  expect_equal(grob$widths, direct$widths)
  expect_equal(grob$heights, direct$heights)
  expect_error(
    as_ideogram_grob(ggplot2::ggplot()),
    "ggideogram"
  )
})
