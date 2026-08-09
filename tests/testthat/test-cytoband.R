cytoband_rows <- function(chr = "A", ...) {
  data.frame(Chr = chr, ..., stringsAsFactors = FALSE)
}
cytoband_bands <- cytoband_rows(
  Start = c(1, 301, 601), End = c(300, 600, 1000),
  Name = c("p11", "cen", "q11"),
  Stain = c("gneg", "acen", "gpos50")
)

test_that("cytoband palettes interpolate stains and accept overrides", {
  expect_equal(
    cytoband_colours(c(
      "gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen"
    )),
    c("#FFFFFF", "#C8C8C8", "#C8C8C8", "#828282", "#000000", "#D92F27")
  )
  expect_equal(
    cytoband_colours(
      c("gpos33", "gpos66"), scheme = "biovizbase"
    ),
    c("#AAAAAA", "#565656")
  )
  expect_equal(
    cytoband_colours(
      c("gpos100", "gpos50"), palette = c(gpos = "#08306B")
    ),
    c("#08306B", "#8397B5")
  )
  expect_equal(cytoband_colours("gpos100", bleach = 0.5), "#7F7F7F")
  expect_error(cytoband_colours("unknown"), "Unknown Giemsa stain")
  expect_error(cytoband_colours("gneg", bleach = 2), "between 0 and 1")
})

test_that("UCSC cytobands are read, selected and converted to a karyotype", {
  ucsc <- data.frame(
    chromosome = c("chr1", "chr2", "chr1"),
    start = c(0, 0, 100), end = c(100, 50, 200),
    name = "band", stain = c("gneg", "gneg", "acen")
  )
  bands <- read_cytoband(ucsc, chr = c("chr2", "chr1"))
  expect_named(bands, c("Chr", "Start", "End", "Name", "Stain"))
  expect_equal(bands$Chr, c("chr2", "chr1", "chr1"))
  expect_equal(bands$Start, c(1, 1, 101))

  karyotype <- cytoband_karyotype(bands)
  expect_equal(karyotype$Chr, c("chr2", "chr1"))
  expect_equal(karyotype$End, c(50, 200))
  expect_equal(karyotype$CE_start, c(0, 101))
  expect_equal(karyotype$CE_end, c(0, 200))
  expect_no_error(as_ideogram_data(karyotype))
  expect_error(read_cytoband(ucsc, chr = "chrZ"), "Not in")
})

test_that("cytoband tables enter the dimensionless renderer directly", {
  raw <- rbind(
    cytoband_bands,
    cytoband_rows(
      "B", Start = 1, End = 800, Name = "p1", Stain = "gneg"
    )
  )
  p <- ggideogram(raw, show_names = FALSE)
  built <- ggplot2::ggplot_build(p)

  expect_s3_class(p$coordinates$layout, "ideogram_layout_v2")
  expect_equal(length(built$data), 3L)
  expect_setequal(
    unique(built$data[[2]]$fill),
    cytoband_colours(raw$Stain)
  )
  expect_no_error(ggplot2::ggplotGrob(p))

  added <- p + geom_chr_cytoband(bleach = 0.4)
  expect_equal(length(added$layers), length(p$layers) + 1L)
  expect_no_error(ggplot2::ggplotGrob(added))

  plain <- ggideogram(
    data.frame(Chr = "A", Start = 0, End = 1000),
    show_names = FALSE
  )
  expect_error(plain + geom_cytoband(), "no cytoband data")
})

test_that("cytoband validation happens before drawing", {
  bad <- cytoband_bands
  bad$Stain[2] <- "gpso50"
  expect_error(
    ggideogram(cytoband_karyotype(bad), cytoband = bad),
    "gpso50"
  )
  expect_no_error(
    ggideogram(
      cytoband_karyotype(bad), cytoband = bad,
      cytoband_palette = c(gpso50 = "red")
    )
  )

  raw <- cytoband_bands
  names(raw)[names(raw) == "Stain"] <- "gieStain"
  expect_equal(cytoband_karyotype(raw), cytoband_karyotype(cytoband_bands))
  expect_no_error(ggideogram(raw))
})
