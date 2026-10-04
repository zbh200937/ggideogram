# Windowed feature counts from a GFF. The three cases the original gets wrong --
# a chromosome shorter than one window, a length that is an exact multiple of the
# window, and a chromosome with no features -- are each pinned here.

gff <- function(chr, pos, type = "gene") {
  data.frame(V1 = chr, V2 = ".", V3 = type, V4 = pos, stringsAsFactors = FALSE)
}

test_that("a chromosome shorter than one window gets a single truncated window", {
  # RIdeogram builds its breaks with seq(window, end, by = window), which is a
  # "wrong sign in 'by' argument" error the moment `end < window`. Any assembly
  # with a short scaffold in it hit this.
  out <- GFFex(gff("A", c(10, 20, 400000)),
               data.frame(Chr = "A", End = 500000), window = 1e6)
  expect_equal(nrow(out), 1L)
  expect_equal(out$Start, 1)
  expect_equal(out$End, 500000)
  expect_equal(out$Value, 3L)
})

test_that("a length that is an exact multiple of the window is not duplicated", {
  out <- GFFex(gff("A", c(1, 1000001)),
               data.frame(Chr = "A", End = 2e6), window = 1e6)
  expect_equal(out$End, c(1e6, 2e6))
  expect_equal(out$Value, c(1L, 1L))
})

test_that("a chromosome absent from the GFF still gets its zero windows", {
  out <- GFFex(gff("A", 10),
               data.frame(Chr = c("A", "B"), End = c(1e6, 2.5e6)), window = 1e6)
  expect_equal(out$Chr, c("A", rep("B", 3)))
  expect_equal(out$Value, c(1L, 0L, 0L, 0L))
  # The short last window of B is truncated to B's length, not rounded up.
  expect_equal(out$End[4], 2.5e6)
})

test_that("windows are half-open on the left, so a feature on an edge counts once", {
  out <- GFFex(gff("A", c(100, 1e6, 1000001)),
               data.frame(Chr = "A", End = 2e6), window = 1e6)
  expect_equal(out$Value, c(2L, 1L))
  expect_equal(sum(out$Value), 3L)
})

test_that("only the requested feature type is counted and unknown chromosomes are errors", {
  g <- rbind(gff("A", c(10, 20)), gff("A", 30, type = "exon"))
  out <- GFFex(g, data.frame(Chr = "A", End = 1e6))
  expect_equal(out$Value, 2L)
  expect_equal(GFFex(g, data.frame(Chr = "A", End = 1e6), feature = "exon")$Value,
               1L)
  expect_error(GFFex(rbind(g, gff("Q", 40)), data.frame(Chr = "A", End = 1e6)),
    "unknown")
})

test_that("GFF counts use source starts and validate complete selected intervals", {
  k <- data.frame(Chr = "A", Start = 100, End = 225)
  g <- gff("A", c(150, 151, 225))
  g$V5 <- c(175, 175, 225)
  out <- GFFex(g, k, window = 50)
  expect_equal(out$Start, c(101, 151, 201))
  expect_equal(out$End, c(150, 200, 225))
  expect_equal(out$Value, c(1, 1, 1))
  for (bad in list(gff("A", 100), gff("A", 226), gff("A", 150.5))) {
    expect_error(GFFex(bad, k, window = 50), "closed integer")
  }
  g$V5[1] <- 226
  expect_error(GFFex(g, k, window = 50), "bounds")
  g$V5[1] <- 149
  expect_error(GFFex(g, k, window = 50), "closed integer")
  expect_error(GFFex(gff("A", NA_real_), k), "finite")
  expect_equal(GFFex(gff("A", 150, "exon"), k, window = 50)$Value, c(0, 0, 0))
})

test_that("GFF counts preserve mapped source karyotypes and composite keys", {
  k <- as_ideogram_data(data.frame(Genome = c("a", "b"), Chrom = "1",
    Lower = 100, Upper = 225), ggplot2::aes(chr = Chrom, start = Lower,
      end = Upper, genome = Genome))
  out <- GFFex(gff(chr_key("a", "1"), 150), k, window = 50)
  expect_equal(out$Chr, rep(k$karyotype$.chr, each = 3))
  expect_equal(out$Value, c(1, 0, 0, 0, 0, 0))
})

test_that("GFFex shares the start-count engine and metadata for complete intervals", {
  g <- gff("A", c(950, 100))
  g$V5 <- c(1150, 1100)
  original <- g
  k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(2000, 500))
  d <- data.frame(Chr = g$V1, Start = g$V4, End = g$V5)
  out <- GFFex(g, k, window = 1000)
  expect_equal(out, bin_genome(d, k, window = 1000, count_position = "start"))
  expect_equal(out$Value, c(2, 0, 0))
  expect_equal(bin_genome(d, k, window = 1000, method = "count")$Value, c(1, 1, 0))
  expect_identical(attr(out, "window_summary")$count_position, "start")
  expect_identical(g, original)
})

test_that("bad arguments are named rather than failing inside cut()", {
  kar1 <- data.frame(Chr = "A", End = 1e6)
  expect_error(GFFex(gff("A", 10)[, 1:2], kar1), "at least 4 GFF columns")
  expect_error(GFFex(gff("A", 10), data.frame(Chr = "A")), "columns Chr and End")
  expect_error(GFFex(gff("A", 10), kar1, window = 0), "positive integer")
  expect_error(GFFex(gff("A", 10), kar1, window = c(1e6, 2e6)),
               "positive integer")
  expect_error(GFFex(gff("A", 10), kar1, window = Inf), "positive integer")
  expect_error(GFFex(gff("A", 10), kar1, window = 1.5), "positive integer")
  expect_error(GFFex(gff("A", 10), data.frame(Chr = c("A", "B"), End = c(1e6, NA))),
               "finite")
})
