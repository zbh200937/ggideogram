draw_circular_gap_test <- function(plot, mm = 120) {
  file <- tempfile(fileext = '.pdf')
  grDevices::cairo_pdf(file, width = mm / 25.4, height = mm / 25.4)
  on.exit({ grDevices::dev.off(); unlink(file) })
  built <- ggplot2::ggplot_build(plot)
  warnings <- character()
  withCallingHandlers(grid::grid.draw(ggplot2::ggplot_gtable(built)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart('muffleWarning')
    })
  state <- get('.circular_gap_text', built$layout$panel_params[[1]]$ideogram_labels)
  list(built = built, state = state, positions = state$cache$positions,
    boxes = state$cache$boxes, warnings = warnings)
}

circular_gap_test_plot <- function(clockwise = TRUE, reverse = FALSE, track_width = 4, track_gap = 1.5, ...) {
  k <- data.frame(Chr = 'A', Start = 1000, End = 100000)
  d <- data.frame(Chr = 'A', Pos = c(1000, 50500, 100000),
    Value = c(0, .5, 1), ID = c('source-a', 'source-b', 'source-c'))
  scopes <- list(
    count = geom_track(side = 'outer', width = track_width, gap = track_gap, limits = c(0, 1),
      data = d, mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value),
      geom = ggplot2::geom_line(), label = 'Count', axis = TRUE),
    coverage = geom_track(side = 'outer', width = track_width, gap = track_gap, limits = c(0, 1),
      data = d, mapping = ggplot2::aes(chr = Chr, x = Pos, y = Value),
      geom = ggplot2::geom_point(), label = 'Cover', axis = TRUE))
  ggideogram(k, orientation = 'circular', show_names = FALSE,
    tracks = scopes, clockwise = clockwise,
    reverse_chr = if (reverse) 'A' else character(), ...)
}

circular_gap_test_overlap <- function(a, b) {
  edges <- rbind(a[c(2, 3, 4, 1), ] - a, b[c(2, 3, 4, 1), ] - b)
  normals <- cbind(-edges[, 2], edges[, 1])
  all(apply(normals, 1, function(normal) {
    pa <- range(a %*% normal); pb <- range(b %*% normal)
    min(pa[2], pb[2]) > max(pa[1], pb[1]) + 1e-8
  }))
}

test_that('gap labels keep their track radii and physical text dimensions', {
  skip_if_not(capabilities('cairo'))
  for (clockwise in c(FALSE, TRUE)) for (reverse in c(FALSE, TRUE)) {
    p <- circular_gap_test_plot(clockwise, reverse, opening_angle = 70)
    results <- lapply(c(120, 200), function(mm) draw_circular_gap_test(p, mm))
    for (r in results) {
      z <- r$positions
      expect_length(r$warnings, 0)
      expect_true(r$state$cache$fits)
      expect_equal(z$label[!z$title], c('0', '1', '0', '1'))
      expect_equal(z$label[z$title], c('Count', 'Cover'))
      expect_equal((z$x - z$original_x) * -z$ny +
        (z$y - z$original_y) * z$nx, rep(0, nrow(z)), tolerance = 1e-10)
      expect_true(all(z$shift <= 2 * z$size))
      for (i in seq_len(length(r$boxes) - 1L)) for (j in (i + 1L):length(r$boxes))
        expect_false(circular_gap_test_overlap(r$boxes[[i]], r$boxes[[j]]))
      for (id in seq_along(r$state$entries)) {
        entry <- r$state$entries[[id]]
        text <- r$state$cache$texts[[id]]
        expect_s3_class(text, 'text')
        expect_equal(text$label, entry$text$label)
        expect_equal(text$gp$fontsize, entry$text$gp$fontsize)
        expect_equal(text$rot, entry$text$rot)
        if (entry$title) expect_null(r$state$cache$leaders[[id]])
        lead <- r$state$cache$leaders[[id]]
        if (!is.null(lead)) {
          expect_s3_class(r$state$cache$segments[[id]], 'segments')
          expect_false(anyDuplicated(lead$row) > 0)
        }
      }
    }
    expect_equal(results[[1]]$positions$width, results[[2]]$positions$width)
    expect_equal(results[[1]]$positions$height, results[[2]]$positions$height)
  }
})

test_that('small openings report capacity without moving titles across tracks', {
  skip_if_not(capabilities('cairo'))
  p <- circular_gap_test_plot(track_width = 2, track_gap = .2, opening_angle = 2)
  r <- draw_circular_gap_test(p, 120)
  z <- r$positions
  expect_length(r$warnings, 1)
  expect_match(r$warnings, 'exceed the closing gap')
  expect_false(r$state$cache$fits)
  expect_equal(z$label, c('Count', '0', '1', 'Cover', '0', '1'))
  expect_equal(z$size, c(2.8, 2.2, 2.2, 2.8, 2.2, 2.2))
  expect_equal((z$x - z$original_x) * -z$ny + (z$y - z$original_y) * z$nx,
    rep(0, nrow(z)), tolerance = 1e-10)
  expect_true(all(z$shift <= 2 * z$size))
  for (id in which(vapply(r$state$entries, function(e) e$title, logical(1))))
    expect_null(r$state$cache$leaders[[id]])
})

test_that('rotated openings and track replay preserve original anchors', {
  skip_if_not(capabilities('cairo'))
  p <- circular_gap_test_plot(start_angle = 73, opening_angle = 70)
  before <- draw_circular_gap_test(p, 200)
  q <- p + geom_track(track = 'count', label = 'Genes')
  after <- draw_circular_gap_test(q, 200)
  expect_length(after$warnings, 0)
  expect_false(identical(before$state, after$state))
  expect_equal(length(after$state$entries), 4L)
  expect_equal(after$positions$label[after$positions$title], c('Genes', 'Cover'))
  expect_equal(after$positions$anchor_x, before$positions$anchor_x)
  expect_equal(after$positions$anchor_y, before$positions$anchor_y)
  z <- after$positions
  expect_equal((z$x - z$original_x) * -z$ny + (z$y - z$original_y) * z$nx,
    rep(0, nrow(z)), tolerance = 1e-10)
})
