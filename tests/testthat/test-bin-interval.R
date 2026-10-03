test_that('coverage counts a union and uses short terminal windows', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = c(25, 5))
  d <- data.frame(Chr = 'A', Start = c(1, 5, 9, 23), End = c(6, 8, 12, 25))
  out <- bin_genome(d, k, window = 10, method = 'coverage')
  expect_equal(out$Value, c(1, 0.2, 0.6, 0))
  expect_equal(out$Width, c(10, 10, 5, 5))
  expect_equal(out$N, c(3, 1, 1, 0))
  expect_equal(attr(out, 'window_summary')$overlap, 'union')
  counts <- bin_genome(d, k, window = 10, method = 'count')
  expect_equal(counts$Value, c(2, 1, 1, 0))
  expect_equal(counts$Rate[3], 2e5)
  expect_equal(counts$Value, bin_genome(d, k, window = 10)$Value)
  p <- ggideogram(k, tracks = track_layout(c = track('right', limits = c(0, 1)))) +
    geom_track_point(data = transform(out, Pos = (Start + End) / 2), track = 'c',
      ggplot2::aes(chr = Chr, position = Pos, value = Value))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that('weighted intervals retain independent overlapping observations', {
  k <- data.frame(Chr = 'A', Start = 0, End = 25)
  d <- data.frame(Chr = 'A', Start = c(1, 6, 21), End = c(10, 15, 25), Score = c(2, 8, 0))
  out <- bin_genome(d, k, window = 10, method = 'weighted_mean', value = 'Score')
  expect_equal(out$Value, c(4, 8, 0))
  expect_equal(out$N, c(2, 1, 1))
  d$Score[1] <- NA_real_
  missing <- bin_genome(d, k, window = 10, method = 'weighted_mean', value = 'Score')
  expect_true(is.na(missing$Value[1]))
  expect_equal(missing$N_valid, c(1, 1, 1))
  removed <- bin_genome(d, k, window = 10, method = 'weighted_mean', value = 'Score', na.rm = TRUE)
  expect_equal(removed$Value, c(8, 8, 0))
  d$Score[] <- NA_real_
  all_missing <- bin_genome(d, k, window = 10, method = 'weighted_mean', value = 'Score', na.rm = TRUE)
  expect_true(all(is.na(all_missing$Value)))
  expect_equal(all_missing$N_valid, c(0, 0, 0))
})

test_that('multiple group fields are independent, including missing groups', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = 12)
  d <- data.frame(Chr = 'A', Start = c(1, 3, 11), End = c(2, 5, 12),
    Sample = c('s1', 's1', NA), Type = c('x', 'y', 'x'))
  out <- bin_genome(d, k, window = 10, method = 'coverage', group = c('Sample', 'Type'))
  expect_equal(nrow(out), 12)
  expect_equal(out$Value, c(0.2, 0, 0, 0, 0.3, 0, 0, 0, 0, 1, 0, 0))
  expect_true(all(is.na(out$Sample[9:12])))
  expect_equal(out$N_valid, out$N)
  shifted <- bin_genome(data.frame(Chr = 'A', Start = 11, End = 12),
    data.frame(Chr = 'A', Start = 10, End = 22), window = 10, method = 'coverage')
  expect_equal(shifted$Start, c(11, 21))
  expect_equal(shifted$Value, c(0.2, 0))
})

test_that('explicit window methods reject invalid biological ranges', {
  k <- data.frame(Chr = 'A', Start = 0, End = 10)
  d <- data.frame(Chr = 'A', Start = 0, End = 2)
  expect_error(bin_genome(d, k, method = 'coverage'), 'bounds')
  d$Start <- 1; d$Chr <- 'B'
  expect_error(bin_genome(d, k, method = 'coverage'), 'unknown')
  d$Chr <- 'A'; d$End <- 11
  expect_error(bin_genome(d, k, method = 'count'), 'bounds')
  expect_error(bin_genome(d, k, method = 'weighted_mean'), 'bounds')
  d$End <- 2
  expect_error(bin_genome(d, k, method = 'weighted_mean'), 'numeric')
  expect_error(bin_genome(d, k, method = 'coverage', group = 'Chr'), 'conflict')
  expect_error(bin_genome(d, k, window = 1.5, method = 'coverage'), 'integer')
})
