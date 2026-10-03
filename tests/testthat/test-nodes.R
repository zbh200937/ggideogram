node_kar <- data.frame(Chr = c('A', 'B'), Start = 0, End = c(100, 80))

test_that('nodes retain source fields, unique ids and exact closed-interval anchors', {
  d <- data.frame(Gene = c('one', 'two', 'coincident'), Chr = 'A',
    Start = c(10L, 21L, 10L), End = c(10L, 30L, 10L), Strand = c('+', '-', '.'))
  n <- chr_nodes(d, ggplot2::aes(id = Gene, chr = Chr, start = Start, end = End))
  expect_s3_class(n, 'chr_nodes')
  expect_identical(lapply(n[names(d)], identity), as.list(d))
  expect_identical(n$.node_id, d$Gene)
  expect_equal(n$.node_position, c(9.5, 25, 9.5))
  expect_identical(n$.node_start, as.numeric(d$Start))
  expect_identical(n$.node_end, as.numeric(d$End))
  expect_s3_class(n[1:2, ], 'chr_nodes')
  p <- chr_nodes(data.frame(ID = c('a', 'b'), Chr = 'A', Pos = c(0, 0)),
    ggplot2::aes(id = ID, chr = Chr, position = Pos))
  expect_equal(p$.node_position, c(0, 0))
  d$Gene[3] <- 'one'
  expect_error(chr_nodes(d, ggplot2::aes(id = Gene, chr = Chr, start = Start, end = End)), 'unique')
  expect_error(chr_nodes(d, ggplot2::aes(id = Gene, chr = Chr, start = Start)), 'end')
  expect_error(chr_nodes(d, ggplot2::aes(id = Gene, chr = Chr, position = Start, start = Start, end = End)), 'either')
  d$Gene[3] <- 'three'; d$Start[1] <- 0
  expect_error(chr_nodes(d, ggplot2::aes(id = Gene, chr = Chr, start = Start, end = End)), 'one-based')
  expect_error(chr_nodes(data.frame(ID = '', Chr = 'A', Pos = 1),
    ggplot2::aes(id = ID, chr = Chr, position = Pos)), 'node identifiers')
})

test_that('id connections share manual source projections and retain edge and node metadata', {
  source <- data.frame(Gene = c('a', 'b', 'same'), Chr = c('A', 'B', 'A'),
    Pos = c(10, 30, 10), Strand = c('+', '-', '.'))
  nodes <- chr_nodes(source, ggplot2::aes(id = Gene, chr = Chr, position = Pos))
  edges <- data.frame(From = c('a', 'same', 'a'), To = c('b', 'a', 'b'),
    Score = c(0.3, 0.7, 0.9), Class = c('cross', 'local', 'cross'))
  manual <- edges
  for (j in 1:2) {
    i <- match(edges[[c('From', 'To')[j]]], source$Gene)
    manual[[paste0('Chr', j)]] <- source$Chr[i]
    manual[[paste0('Pos', j)]] <- source$Pos[i]
  }
  for (direction in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(node_kar, orientation = direction, reverse_chr = 'A')
    by_id <- base + geom_chr_connection(ggplot2::aes(from = From, to = To,
      colour = Class, alpha = Score), edges, nodes = nodes)
    by_bp <- base + geom_chr_connection(ggplot2::aes(chr1 = Chr1, position1 = Pos1,
      chr2 = Chr2, position2 = Pos2, colour = Class, alpha = Score), manual)
    id <- utils::tail(by_id$layers, 1)[[1]]; bp <- utils::tail(by_bp$layers, 1)[[1]]
    expect_identical(class(id$geom)[1], 'GeomSegment')
    cols <- c('.pair_id', '.pair_chr1', '.pair_start1', '.pair_chr2', '.pair_start2',
      '.pair_x', '.pair_y', '.pair_xend', '.pair_yend')
    expect_equal(id$data[cols], bp$data[cols])
    first <- id$data[!duplicated(id$data$.pair_id), , drop = FALSE]
    expect_equal(first$.node_id1, edges$From)
    expect_equal(first$.node_id2, edges$To)
    expect_equal(first$Score, edges$Score)
    expect_equal(first$.node1_Strand, source$Strand[match(edges$From, source$Gene)])
    built_id <- utils::tail(ggplot2::ggplot_build(by_id)$data, 1)[[1]]
    built_bp <- utils::tail(ggplot2::ggplot_build(by_bp)$data, 1)[[1]]
    expect_equal(built_id[c('x', 'y', 'xend', 'yend', 'colour', 'alpha', 'group')],
      built_bp[c('x', 'y', 'xend', 'yend', 'colour', 'alpha', 'group')])
    expect_no_error(ggplot2::ggplotGrob(by_id))
  }
})

test_that('full-source interval nodes use existing joint clipping and reversed correspondence', {
  source <- data.frame(ID = c('left', 'right', 'other'), Chr = c('A', 'A', 'B'),
    Start = c(11, 51, 5), End = c(50, 90, 10))
  nodes <- chr_nodes(source, ggplot2::aes(id = ID, chr = Chr, start = Start, end = End))
  edges <- data.frame(From = c('left', 'left'), To = c('right', 'other'), Direction = '-')
  for (direction in c('vertical', 'horizontal', 'circular')) {
    base <- ggideogram(chr_view(node_kar, 'A', 20, 80), orientation = direction, reverse_chr = 'A')
    p <- base + geom_chr_synteny(ggplot2::aes(from = From, to = To, orientation = Direction),
      edges, nodes = nodes)
    layer <- utils::tail(p$layers, 1)[[1]]
    expect_identical(class(layer$geom)[1], 'GeomPolygon')
    d <- layer$data
    expect_equal(unique(d$.pair_id), 1)
    expect_equal(unique(d$.pair_from1), 20)
    expect_equal(unique(d$.pair_to1), 50)
    expect_equal(unique(d$.pair_from2), 80)
    expect_equal(unique(d$.pair_to2), 50)
    expect_equal(unique(d$.node_start1), 11)
    expect_equal(unique(d$.node_end2), 90)
    expect_equal(nrow(nodes), 3)
    expect_no_error(ggplot2::ggplotGrob(p))
  }
})

test_that('node validation uses full source bounds and explicit multi-genome keys', {
  kar <- data.frame(Genome = c('rice', 'wild'), Assembly = c('IRGSP', 'OR'),
    Chr = '1', Start = 0, End = c(100, 80))
  k <- as_ideogram_data(kar, ggplot2::aes(genome = Genome, assembly = Assembly,
    chr = Chr, start = Start, end = End))
  source <- data.frame(ID = c('rice:gene', 'wild:gene'), Genome = kar$Genome,
    Assembly = kar$Assembly, Chr = '1', Pos = c(30, 40), Gene = 'same-name')
  n <- chr_nodes(source, ggplot2::aes(id = ID, chr = chr_key(Genome, Chr, Assembly), position = Pos))
  e <- data.frame(From = 'rice:gene', To = 'wild:gene')
  p <- ggideogram(k, orientation = 'horizontal', reverse_chr = chr_key('wild', '1', 'OR')) +
    geom_chr_connection(ggplot2::aes(from = From, to = To), e, nodes = n)
  d <- utils::tail(p$layers, 1)[[1]]$data
  expect_equal(d$.pair_chr1, chr_key('rice', '1', 'IRGSP'))
  expect_equal(d$.pair_chr2, chr_key('wild', '1', 'OR'))
  expect_equal(d$.node1_Gene, 'same-name')
  expect_equal(d$.node2_Gene, 'same-name')
  expect_equal(unproject_chr_point(p$coordinates$layout, d,
    '.pair_chr2', '.pair_xend', '.pair_yend')$.position, 40)
  expect_no_error(ggplot2::ggplotGrob(p))
  e$To <- 'missing'
  expect_error(ggideogram(k) + geom_chr_connection(ggplot2::aes(from = From, to = To),
    e, nodes = n), 'Unknown node.*missing')
  e$To <- 'rice:gene'
  bad <- source; bad$Pos[2] <- 81
  bad_nodes <- chr_nodes(bad, ggplot2::aes(id = ID, chr = chr_key(Genome, Chr, Assembly), position = Pos))
  local <- ggideogram(chr_view(k, chr_key('rice', '1', 'IRGSP'), 20, 50))
  expect_error(local + geom_chr_connection(ggplot2::aes(from = From, to = To),
    e, nodes = bad_nodes), 'source bounds')
  bad <- chr_nodes(source, ggplot2::aes(id = ID, chr = Chr, position = Pos))
  expect_error(ggideogram(k) + geom_chr_connection(ggplot2::aes(from = From, to = To),
    e, nodes = bad), 'unknown source chromosomes')
  expect_error(ggideogram(k) + geom_chr_synteny(ggplot2::aes(from = From, to = To),
    e, nodes = n), 'interval nodes')
})
