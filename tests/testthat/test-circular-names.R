circular_name_text_boxes <- function(plot, mm) {
  file <- tempfile(fileext = '.pdf')
  grDevices::cairo_pdf(file, width = mm / 25.4, height = mm / 25.4)
  on.exit({ grDevices::dev.off(); unlink(file) })
  grid::pushViewport(grid::viewport(width = grid::unit(mm, 'mm'),
    height = grid::unit(mm, 'mm'), xscale = c(0, 1), yscale = c(0, 1)))
  b <- ggplot2::ggplot_build(plot)
  draw <- function(geom_class) {
    i <- which(vapply(b$plot$layers, function(layer)
      inherits(layer$geom, geom_class), logical(1)))
    layer <- b$plot$layers[[i]]
    do.call(layer$geom$draw_panel, c(list(data = b$data[[i]],
      panel_params = b$layout$panel_params[[1]], coord = b$plot$coordinates),
      layer$computed_geom_params[intersect(names(layer$computed_geom_params),
        layer$geom$parameters())]))
  }
  boxes <- function(grob) lapply(seq_along(grob$label), function(i) {
    gp <- do.call(grid::gpar, lapply(grob$gp, function(value)
      value[(i - 1L) %% length(value) + 1L]))
    text <- grid::textGrob(grob$label[i], gp = gp)
    w <- grid::convertWidth(grid::grobWidth(text), 'mm', valueOnly = TRUE)
    h <- grid::convertHeight(grid::grobHeight(text), 'mm', valueOnly = TRUE)
    anchor <- c(grid::convertX(grob$x[i], 'mm', valueOnly = TRUE),
      grid::convertY(grob$y[i], 'mm', valueOnly = TRUE))
    angle <- grob$rot[(i - 1L) %% length(grob$rot) + 1L] * pi / 180
    hjust <- grob$hjust[(i - 1L) %% length(grob$hjust) + 1L]
    vjust <- grob$vjust[(i - 1L) %% length(grob$vjust) + 1L]
    corners <- rbind(c(-hjust * w, -vjust * h),
      c((1 - hjust) * w, -vjust * h),
      c((1 - hjust) * w, (1 - vjust) * h),
      c(-hjust * w, (1 - vjust) * h))
    rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
    list(label = grob$label[i], corners = sweep(corners %*% rotation, 2, anchor, '+'),
      height = h, fontsize = gp$fontsize)
  })
  name_index <- which(vapply(b$plot$layers, function(layer)
    inherits(layer$geom, 'GeomCircularName'), logical(1)))
  anchors <- b$plot$coordinates$transform(b$data[[name_index]],
    b$layout$panel_params[[1]])
  axis_index <- which(vapply(b$plot$layers, function(layer)
    inherits(layer$geom, 'GeomIdeogramAxisText'), logical(1)))
  list(names = boxes(draw('GeomCircularName')),
    axes = if (length(axis_index)) boxes(draw('GeomIdeogramAxisText')) else list(),
    anchors = data.frame(x = anchors$x * mm, y = anchors$y * mm,
      nx = anchors$nx, ny = anchors$ny))
}

circular_text_boxes_overlap <- function(a, b) {
  edges <- rbind(a[c(2, 3, 4, 1), ] - a, b[c(2, 3, 4, 1), ] - b)
  normals <- cbind(-edges[, 2], edges[, 1])
  all(apply(normals, 1, function(normal) {
    pa <- range(a %*% normal); pb <- range(b %*% normal)
    min(pa[2], pb[2]) > max(pa[1], pb[1]) + 1e-8
  }))
}

test_that('multiline circular genome names clear bp labels at both device sizes', {
  skip_if_not(capabilities('cairo'))
  k <- data.frame(Genome = c('Audit-A', 'Audit-B'), Assembly = 'TAIR10',
    Chr = '1', Start = 0, End = 100000)
  keys <- with(k, chr_key(Genome, Chr, Assembly))
  d <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start, end = End,
    genome = Genome, assembly = Assembly))
  p <- ggideogram(d, orientation = 'circular', radius = 20,
    gap_angle = setNames(c(12, 20), keys), clockwise = FALSE,
    start_angle = 31, reverse_chr = keys[2], axis = TRUE,
    axis_units = 'kb', axis_breaks = seq(0, 75000, 25000),
    name_size = 2.9, axis_size = 2.7)
  drawn <- lapply(c(120, 200), function(mm) circular_name_text_boxes(p, mm))
  for (result in drawn) {
    expect_equal(vapply(result$names, `[[`, character(1), 'label'),
      c('Audit-A\n1', 'Audit-B\n1'))
    for (name in result$names) for (axis in result$axes)
      expect_false(circular_text_boxes_overlap(name$corners, axis$corners),
        info = paste(name$label, axis$label))
  }
  expect_equal(vapply(drawn[[1]]$names, `[[`, numeric(1), 'fontsize'),
    vapply(drawn[[2]]$names, `[[`, numeric(1), 'fontsize'))
  expect_equal(vapply(drawn[[1]]$names, `[[`, numeric(1), 'height'),
    vapply(drawn[[2]]$names, `[[`, numeric(1), 'height'))
})

test_that('different line counts share a circle and keep at least their physical gap', {
  skip_if_not(capabilities('cairo'))
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = 100,
    Label = c('Short', 'Genome B\nAssembly 1\nChr B'))
  d <- as_ideogram_data(k, ggplot2::aes(chr = Chr, start = Start,
    end = End, label = Label))
  p <- ggideogram(d, orientation = 'circular', gap_angle = 20,
    start_angle = 31, name_size = 2.9, name_gap = 0.5)
  for (mm in c(120, 200)) {
    result <- circular_name_text_boxes(p, mm)
    centres <- t(vapply(result$names, function(name) colMeans(name$corners), numeric(2)))
    radius <- sqrt(rowSums(sweep(centres, 2, c(mm, mm) / 2, '-')^2))
    expect_equal(radius, rep(radius[1], length(radius)), tolerance = 1e-8)
    heights <- vapply(result$names, `[[`, numeric(1), 'height')
    for (i in seq_along(result$names)) {
      anchor <- result$anchors[i, ]
      corners <- sweep(result$names[[i]]$corners, 2,
        c(anchor$x, anchor$y), '-')
      clearance <- min(corners %*% c(anchor$nx, anchor$ny))
      expect_equal(clearance, (max(heights) - heights[i]) / 2 + 0.5 * 2.9,
        tolerance = 1e-8)
    }
  }
})

test_that('name positions use display endpoints and the final track group', {
  k <- data.frame(Chr = c('A', 'B'), Start = 0, End = c(100, 80))
  tracks <- track_layout(l = track('left', width = 2), r = track('right', width = 6))
  for (orientation in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(k, orientation = orientation, tracks = tracks,
      reverse_chr = 'B', show_names = FALSE)
    for (position in c('auto', 'start', 'end', 'middle')) {
      p <- base + geom_chr_name(position = position)
      d <- utils::tail(ggplot2::ggplot_build(p)$data, 1)[[1]]
      g <- p$coordinates$layout$chrom
      if (orientation == 'circular') {
        expected <- switch(position, start = g$.arc_start, end = g$.arc_end,
          (g$.arc_start + g$.arc_end) / 2)
        expect_equal(d$x, expected)
      } else {
        axis <- if (orientation == 'vertical') 'y' else 'x'
        start <- g[[paste0('.axis_start_', axis)]]
        end <- g[[paste0('.axis_end_', axis)]]
        expected <- if (position == 'middle') (start + end) / 2 else {
          first <- position %in% c('auto', 'start')
          if (xor(first, orientation == 'horizontal')) pmax(start, end) else pmin(start, end)
        }
        expect_equal(d[[axis]], expected)
        expect_s3_class(utils::tail(p$layers, 1)[[1]]$geom, 'GeomText')
      }
    }
  }
  expect_error(geom_chr_name(position = 'center'), 'arg')
})

test_that('linear endpoint names clear same-end value axes in physical units', {
  skip_if_not(capabilities('cairo'))
  k <- data.frame(Chr = c('1', '14'), Start = 0, End = c(100, 80))
  for (orientation in c('vertical', 'horizontal')) for (names_first in c(FALSE, TRUE)) {
    base <- ggideogram(k, orientation = orientation, reverse_chr = '14', show_names = FALSE,
      tracks = list(density = geom_track(side = 'right', width = 3.5, limits = c(0, 150))))
    names <- list(geom_chr_name(position = 'start'), geom_chr_name(position = 'end'))
    axes <- list(geom_track_axis('density', position = 'start', breaks = c(0, 75, 150)),
      geom_track_axis('density', position = 'end', breaks = c(0, 75, 150)))
    p <- if (names_first) base + names + axes else base + axes + names
    fonts <- list()
    for (mm in c(120, 200)) {
      file <- tempfile(fileext = '.pdf')
      grDevices::cairo_pdf(file, width = mm / 25.4, height = mm / 25.4)
      grid::pushViewport(grid::viewport(width = grid::unit(mm, 'mm'),
        height = grid::unit(mm, 'mm'), xscale = c(0, 1), yscale = c(0, 1)))
      b <- ggplot2::ggplot_build(p)
      draw <- function(i) {
        layer <- b$plot$layers[[i]]
        do.call(layer$geom$draw_panel, c(list(data = b$data[[i]],
          panel_params = b$layout$panel_params[[1]], coord = b$plot$coordinates),
          layer$computed_geom_params[intersect(names(layer$computed_geom_params), layer$geom$parameters())]))
      }
      boxes <- function(grob) lapply(seq_along(grob$label), function(i) {
        at <- function(x) x[(i - 1L) %% length(x) + 1L]
        gp <- do.call(grid::gpar, lapply(grob$gp, at))
        text <- grid::textGrob(grob$label[i], gp = gp)
        w <- grid::convertWidth(grid::grobWidth(text), 'mm', TRUE)
        h <- grid::convertHeight(grid::grobHeight(text), 'mm', TRUE)
        x <- grid::convertX(grob$x[i], 'mm', TRUE)
        y <- grid::convertY(grob$y[i], 'mm', TRUE)
        c(x - at(grob$hjust) * w, x + (1 - at(grob$hjust)) * w,
          y - at(grob$vjust) * h, y + (1 - at(grob$vjust)) * h)
      })
      name_index <- which(vapply(b$plot$layers, function(layer)
        inherits(layer$geom, 'GeomChromosomeName'), logical(1)))
      axis_index <- which(vapply(b$plot$layers, function(layer)
        inherits(layer$geom, 'GeomIdeogramAxisText'), logical(1)))
      name_grobs <- lapply(name_index, draw)
      name_boxes <- unlist(lapply(name_grobs, boxes), recursive = FALSE)
      axis_boxes <- unlist(lapply(axis_index, function(i) boxes(draw(i))), recursive = FALSE)
      overlaps <- vapply(name_boxes, function(a) any(vapply(axis_boxes, function(z)
        min(a[2], z[2]) > max(a[1], z[1]) + 1e-8 &&
          min(a[4], z[4]) > max(a[3], z[3]) + 1e-8, logical(1))), logical(1))
      expect_false(any(overlaps), info = paste(orientation, names_first, mm))
      fonts[[as.character(mm)]] <- unlist(lapply(name_grobs, function(grob) grob$gp$fontsize))
      anchors <- linear_chromosome_name_data(p$coordinates$layout, 'start')
      expect_equal(b$data[[name_index[1]]]$x, anchors$x)
      expect_equal(b$data[[name_index[1]]]$y, anchors$y)
      grid::popViewport()
      grDevices::dev.off()
      unlink(file)
    }
    expect_equal(fonts[['120']], fonts[['200']])
  }
})

test_that('selected chromosome components preserve other name-axis clearance', {
  k <- data.frame(Chr = c('1', '14'), Start = 0, End = c(100, 80))
  p <- ggideogram(k, orientation = 'horizontal', reverse_chr = '14',
    tracks = list(density = geom_track(side = 'right', width = 3.5, limits = c(0, 150)))) +
    geom_track_axis('density', chr = '14', position = 'end', breaks = c(0, 75, 150))
  name <- Filter(function(layer) inherits(layer$geom, 'GeomChromosomeName'), p$layers)[[1]]
  before <- grid::convertWidth(name$geom_params$name_axis_clearance, 'mm', TRUE)
  q <- p + geom_chr(chr = '1', component = 'axis')
  name <- Filter(function(layer) inherits(layer$geom, 'GeomChromosomeName'), q$layers)[[1]]
  after <- grid::convertWidth(name$geom_params$name_axis_clearance, 'mm', TRUE)
  expect_equal(q$coordinates$layout$chrom$.chr, k$Chr)
  expect_equal(after, before)
  expect_gt(after[2], 0)
  expect_no_error(ggplot2::ggplotGrob(q))
})

test_that('track titles contribute their native text extent to endpoint name clearance', {
  k <- data.frame(Chr = c('1', '14'), Start = 0, End = c(100, 80))
  p <- ggideogram(k, orientation = 'vertical', reverse_chr = '14') +
    geom_track(track = 'density', side = 'left', width = 3, label = 'Genes / Mb',
      limits = c(0, 200))
  name <- Filter(function(layer) inherits(layer$geom, 'GeomChromosomeName'), p$layers)[[1]]
  title <- Filter(function(layer) !is.null(layer$ideogram_track_title), p$layers)[[1]]
  for (mm in c(120, 200)) {
    file <- tempfile(fileext = '.pdf')
    grDevices::cairo_pdf(file, width = mm / 25.4, height = mm / 25.4)
    grid::pushViewport(grid::viewport(width = grid::unit(mm, 'mm'),
      height = grid::unit(mm, 'mm'), xscale = c(0, 1), yscale = c(0, 1)))
    b <- ggplot2::ggplot_build(p)
    draw <- function(geom_class, track_title = FALSE) {
      index <- which(vapply(b$plot$layers, function(layer)
        if (track_title) !is.null(layer$ideogram_track_title) else inherits(layer$geom, geom_class), logical(1)))
      layer <- b$plot$layers[[index]]
      do.call(layer$geom$draw_panel, c(list(data = b$data[[index]],
        panel_params = b$layout$panel_params[[1]], coord = b$plot$coordinates),
        layer$computed_geom_params[intersect(names(layer$computed_geom_params), layer$geom$parameters())]))
    }
    text <- draw('GeomChromosomeName')
    title_text <- draw(NULL, TRUE)
    name_height <- grid::convertHeight(grid::grobHeight(grid::textGrob(text$label[1],
      gp = text$gp)), 'mm', TRUE)
    title_height <- grid::convertHeight(grid::grobHeight(grid::textGrob(title_text$label,
      gp = title_text$gp)), 'mm', TRUE)
    name_bottom <- grid::convertY(text$y[1], 'mm', TRUE) - name_height / 2
    title_top <- grid::convertY(title_text$y, 'mm', TRUE) + title_height
    expect_gte(name_bottom - title_top, name$geom_params$name_gap * name$aes_params$size - 1e-8)
    expect_equal(text$gp$fontsize, rep(ideogram_text_size('chromosome') * ggplot2::.pt, 2))
    grid::popViewport(); grDevices::dev.off(); unlink(file)
  }
})
