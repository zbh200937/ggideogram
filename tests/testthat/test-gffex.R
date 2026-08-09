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

test_that("only the requested feature type on a known chromosome is counted", {
  g <- rbind(gff("A", c(10, 20)), gff("A", 30, type = "exon"), gff("Q", 40))
  out <- GFFex(g, data.frame(Chr = "A", End = 1e6))
  expect_equal(out$Value, 2L)
  expect_equal(GFFex(g, data.frame(Chr = "A", End = 1e6), feature = "exon")$Value,
               1L)
})

test_that("bad arguments are named rather than failing inside cut()", {
  kar1 <- data.frame(Chr = "A", End = 1e6)
  expect_error(GFFex(gff("A", 10)[, 1:2], kar1), "at least 4 GFF columns")
  expect_error(GFFex(gff("A", 10), data.frame(Chr = "A")), "columns Chr and End")
  expect_error(GFFex(gff("A", 10), kar1, window = 0), "single positive number")
  expect_error(GFFex(gff("A", 10), kar1, window = c(1e6, 2e6)),
               "single positive number")
  expect_error(GFFex(gff("A", 10), data.frame(Chr = c("A", "B"), End = c(1e6, NA))),
               "positive and non-missing.*`B`")
})
