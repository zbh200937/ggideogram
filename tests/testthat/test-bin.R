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
               c(NA, 99, 0, 0))
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q", FUN = mean,
                          na.rm = TRUE)$Value, c(50, 99, NA, NA))
})

test_that("empty fill does not replace a missing summary from observed rows", {
  d <- data.frame(Chr = "A", Start = c(10, 1500), Q = c(NA_real_, 0))
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q", empty = -1)$Value,
    c(NA, 0, -1, -1))
  expect_equal(bin_genome(d, kar, window = 1000, value = "Q", empty = 0,
    na.rm = TRUE)$Value, c(NaN, 0, 0, 0))
})

test_that("summary metadata preserves lazy unused arguments and exact na.rm names", {
  d <- data.frame(Chr = "A", Start = c(10, 1500), Q = c(NA_real_, 0))
  lazy_sum <- function(x, ...) sum(x)
  out <- bin_genome(d, kar, window = 1000, value = "Q", FUN = lazy_sum,
    ignored = stop("unused argument was evaluated"))
  expect_equal(out$Value, c(NA, 0, NA, NA))
  expect_null(attr(out, "window_summary")$na.rm)
  similar <- bin_genome(d, kar, window = 1000, value = "Q", FUN = lazy_sum,
    na.rm_extra = TRUE)
  expect_equal(similar$Value, out$Value)
  expect_null(attr(similar, "window_summary")$na.rm)
  means <- bin_genome(d, kar, window = 1000, value = "Q", na.rm = TRUE,
    ignored = stop("unused argument was evaluated"))
  expect_equal(means$Value, c(NaN, 0, NA, NA))
  expect_true(attr(means, "window_summary")$na.rm)
  empty <- bin_genome(d[FALSE, ], kar, window = 1000, value = "Q", FUN = lazy_sum,
    ignored = stop("unused argument was evaluated"))
  expect_true(all(is.na(empty$Value)))
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

test_that("count and FUN summaries use the same selectable representative positions", {
  d <- data.frame(Chr = "A", Start = c(950, 100), End = c(1150, 1100), Q = c(2, 8))
  expected <- list(midpoint = c(1, 1, 0, 0), start = c(2, 0, 0, 0), end = c(0, 2, 0, 0))
  for (position in names(expected)) {
    out <- bin_genome(d, kar, window = 1000, count_position = position)
    explicit <- bin_genome(d, kar, window = 1000, method = "count", count_position = position)
    expect_equal(out$Value, expected[[position]])
    expect_equal(explicit$Value, out$Value)
    expect_equal(attr(out, "window_summary")$count_position, position)
    expect_equal(attr(explicit, "window_summary")$count_position, position)
  }
  means <- bin_genome(d, kar, window = 1000, value = "Q", count_position = "start")
  expect_equal(means$Value, c(5, NA, NA, NA))
  expect_equal(attr(means, "window_summary")$summary, "mean")
  sums <- bin_genome(d, kar, window = 1000, value = "Q", FUN = sum, count_position = "end")
  expect_equal(sums$Value, c(NA, 10, NA, NA))
  expect_equal(attr(sums, "window_summary")$summary, "sum")
  # Preserve the mean of base indices, including an even-length interval.
  even <- data.frame(Chr = "A", Start = 1, End = 2000)
  expect_equal(bin_genome(even, kar, window = 1000)$Value, c(0, 1, 0, 0))
  for (method in c("coverage", "weighted_mean")) {
    expect_error(bin_genome(d, kar, method = method, count_position = "start"),
      "applies only to count or FUN")
  }
})

test_that("all window results record the same metadata fields without guessing value units", {
  d <- data.frame(Chr = "A", Start = c(10, 1500), End = c(20, 1600), Q = c(NA, 0))
  results <- list(
    bin_genome(d, kar, window = 1000),
    bin_genome(d, kar, window = 1000, value = "Q", empty = -1, na.rm = TRUE),
    bin_genome(d, kar, window = 1000, method = "count"),
    bin_genome(d, kar, window = 1000, method = "coverage"),
    bin_genome(d, kar, window = 1000, method = "weighted_mean", value = "Q"))
  metadata <- lapply(results, attr, "window_summary")
  for (m in metadata) {
    expect_named(m, c("method", "group", "coordinates", "unit", "value", "window",
      "count_position", "na.rm", "overlap", "empty", "summary"))
    expect_identical(m$coordinates, "1-based closed")
    expect_equal(m$window, 1000)
  }
  expect_equal(vapply(metadata, `[[`, character(1), "method"),
    c("count", "summary", "count", "coverage", "weighted_mean"))
  expect_identical(metadata[[1]]$unit, "features")
  expect_identical(metadata[[4]]$unit, "fraction")
  expect_null(metadata[[2]]$unit)
  expect_null(metadata[[5]]$unit)
  expect_identical(metadata[[2]]$value, "Q")
  expect_identical(metadata[[5]]$value, "Q")
  expect_true(metadata[[2]]$na.rm)
  expect_equal(metadata[[2]]$empty, -1)
  expect_equal(results[[2]]$Value, c(NaN, 0, -1, -1))
  points <- bin_genome(d[c("Chr", "Start")], kar, window = 1000, count_position = "end")
  expect_identical(attr(points, "window_summary")$count_position, "start")
})

test_that("empty chromosomes remain and invalid source intervals are rejected", {
  d <- data.frame(Chr = "A", Start = 10)
  out <- bin_genome(d, kar, window = 1000)
  expect_true("B" %in% out$Chr)
  expect_equal(sum(out$Value), 1)
  expect_error(bin_genome(data.frame(Chr = c("A", "Z"), Start = 10), kar), "unknown")
  expect_error(bin_genome(data.frame(Chr = "B", Start = 700, End = 900), kar), "bounds")
  expect_error(bin_genome(data.frame(Chr = "A", Start = NA_real_), kar), "finite")
  expect_error(bin_genome(data.frame(Chr = "A", Start = 1.5), kar), "integer")
})

test_that("all window methods respect nonzero chromosome starts", {
  k <- data.frame(Chr = "A", Start = 100, End = 225)
  d <- data.frame(Chr = "A", Start = c(150, 151, 225), End = c(150, 175, 225), Score = 2)
  for (method in list(NULL, "count", "coverage", "weighted_mean")) {
    out <- bin_genome(d, k, window = 50, method = method,
      value = if (identical(method, "weighted_mean")) "Score" else NULL)
    expect_equal(out$Start, c(101, 151, 201))
    expect_equal(out$End, c(150, 200, 225))
    if (is.null(method) || identical(method, "count")) expect_equal(out$Value, c(1, 1, 1))
  }
})

test_that("bad arguments are named", {
  d <- data.frame(Chr = "A", Start = 10)
  expect_error(bin_genome(list(), kar), "must be a data frame")
  expect_error(bin_genome(data.frame(Chr = "A"), kar), "missing column `Start`")
  expect_error(bin_genome(d, kar, value = "Nope"), "not in `data`")
  expect_error(bin_genome(d, kar, window = 0), "positive integer")
  expect_error(bin_genome(d, kar, window = Inf), "positive integer")
  expect_error(bin_genome(d, data.frame(Chr = "A")), "needs columns Chr and End")
  expect_error(bin_genome(d, data.frame(Chr = "A", End = -1)), "End|end")
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
