# bin_genome(): the general form of what GFFex() does for one GFF feature type.

kar <- data.frame(Chr = c("A", "B"), Start = 0, End = c(2500, 800))

test_that("window edges cover the chromosome without duplicating the last one", {
  expect_equal(window_edges(2500, 1000)$upper, c(1000, 2000, 2500))
  expect_equal(window_edges(2500, 1000)$lower, c(0, 1000, 2000))
  # An exact multiple does not produce a zero-width final window.
  expect_equal(window_edges(2000, 1000)$upper, c(1000, 2000))
  # A chromosome shorter than one window is one window, not an error.
  expect_equal(window_edges(800, 1000)$upper, 800)
})

test_that("counting rows is the default, and empty windows come back as 0", {
  d <- data.frame(Chr = "A", Start = c(10, 20, 1500))
  out <- bin_genome(d, kar, window = 1000)
  expect_named(out, c("Chr", "Start", "End", "Value"))
  expect_equal(out$Chr, c("A", "A", "A", "B"))
  expect_equal(out$Value, c(2, 1, 0, 0))
  # Half-open on the left, as GFFex is: 1000 belongs to the first window.
  expect_equal(bin_genome(data.frame(Chr = "A", Start = c(1000, 1001)), kar,
                          window = 1000)$Value, c(1, 1, 0, 0))
})

test_that("a value column is summarised, and empty windows come back as NA", {
  d <- data.frame(Chr = "A", Start = c(10, 20, 1500), Q = c(30, 50, 99))
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q")$Value,
               c(40, 99, NA, NA))
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q", FUN = max)$Value,
               c(50, 99, NA, NA))
  # The fill for an empty window is settable, and extra arguments reach FUN.
  d$Q[1] <- NA
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q", empty = 0)$Value,
               c(0, 99, 0, 0))
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q", FUN = mean,
                          na.rm = TRUE)$Value, c(50, 99, NA, NA))
})

test_that("an interval is binned once, by its midpoint", {
  # 950-1150 straddles the first window edge; its midpoint puts it in the second,
  # and it is counted there once rather than in both.
  d <- data.frame(Chr = "A", Start = 950, End = 1150)
  expect_equal(bin_genome(d, kar, window = 1000)$Value, c(0, 1, 0, 0))
  # Without an End the Start decides, so the same row lands in the first.
  expect_equal(bin_genome(d[c("Chr", "Start")], kar, window = 1000)$Value,
               c(1, 0, 0, 0))
  # A midpoint exactly on an edge goes to the window that edge closes, matching
  # the half-open-on-the-left rule the rest of the binning uses.
  expect_equal(bin_genome(data.frame(Chr = "A", Start = 900, End = 1100), kar,
                          window = 1000)$Value, c(1, 0, 0, 0))
})

test_that("chromosomes absent from the data still appear, and strays are dropped", {
  d <- data.frame(Chr = c("A", "Z"), Start = c(10, 10))
  out <- bin_genome(d, kar, window = 1000)
  expect_true("B" %in% out$Chr)
  expect_false("Z" %in% out$Chr)
  expect_equal(sum(out$Value), 1)
  # Positions off the end of their own chromosome are dropped too.
  expect_equal(sum(bin_genome(data.frame(Chr = "B", Start = 5000), kar,
                              window = 1000)$Value), 0)
})

test_that("bad arguments are named", {
  d <- data.frame(Chr = "A", Start = 10)
  expect_error(bin_genome(list(), kar), "must be a data frame")
  expect_error(bin_genome(data.frame(Chr = "A"), kar), "missing column `Start`")
  expect_error(bin_genome(d, kar, value = "Nope"), "not in `data`")
  expect_error(bin_genome(d, kar, window = 0), "single positive number")
  expect_error(bin_genome(d, data.frame(Chr = "A")), "needs columns Chr and End")
  expect_error(bin_genome(d, data.frame(Chr = "A", End = -1)), "must be positive")
})

test_that("bin_genome reproduces GFFex for the case GFFex covers", {
  gff <- data.frame(V1 = c("A", "A", "B", "A"), V2 = ".",
                    V3 = c("gene", "gene", "gene", "exon"),
                    V4 = c(10, 1500, 700, 20), stringsAsFactors = FALSE)
  a <- GFFex(gff, kar, window = 1000)
  b <- bin_genome(gff[gff$V3 == "gene", c("V1", "V4")] |>
                    stats::setNames(c("Chr", "Start")), kar, window = 1000)
  expect_equal(a, b)
})
