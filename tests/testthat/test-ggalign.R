ggalign_karyotype <- data.frame(
  Chr = c("A", "B", "C"), Start = 0, End = c(100, 80, 50),
  Label = c("Alpha", "Beta", "Gamma")
)

test_that("ggalign ordering and groups preserve chromosome keys and labels", {
  skip_if_not_installed("ggalign", "1.2.0")
  matrix <- matrix(seq_len(9), nrow = 3,
    dimnames = list(ggalign_karyotype$Chr, c("s1", "s2", "s3")))
  semantic <- as_ideogram_data(ggalign_karyotype,
    ggplot2::aes(chr = Chr, start = Start, end = End, label = Label))
  for (axis in c("x", "y")) {
    for (grouped in c(FALSE, TRUE)) {
      for (labels in list(ggplot2::waiver(), c(A = "a", B = "b", C = "c"))) {
        keys <- list()
        parent <- guide_chromosome_axis(semantic)
        guide <- ggplot2::ggproto(NULL, parent,
          extract_key = function(self, scale, aesthetic, ...) {
            key <- parent$extract_key(scale, aesthetic, ...)
            keys[[length(keys) + 1L]] <<- key
            key
          })
        layout <- if (axis == "x") ggalign::stack_discretev(matrix) else
          ggalign::stack_discreteh(matrix)
        layout <- layout + ggalign::align_order(c(3, 1, 2))
        if (grouped) layout <- layout + ggalign::align_group(c("one", "two", "one"))
        plot <- if (axis == "x") {
          layout + ggalign::ggalign() +
            ggplot2::geom_point(ggplot2::aes(x = .discrete_x, y = value)) +
            scale_x_chromosome(semantic, limits = NULL, labels = labels, guide = guide)
        } else {
          layout + ggalign::ggalign() +
            ggplot2::geom_point(ggplot2::aes(x = value, y = .discrete_y)) +
            scale_y_chromosome(semantic, limits = NULL, labels = labels, guide = guide)
        }
        expect_s3_class(ggalign::ggalignGrob(plot), "gtable")
        expect_gt(length(keys), 0L)
        expected <- if (inherits(labels, "waiver")) {
          c(A = "Alpha", B = "Beta", C = "Gamma")
        } else labels
        expect_setequal(unlist(lapply(keys, `[[`, ".chr")), ggalign_karyotype$Chr)
        for (key in keys) {
          if (!grouped) expect_identical(key$.chr, c("C", "A", "B"))
          expect_equal(as.character(key$.label), unname(expected[key$.chr]))
        }
      }
    }
  }
})

ggalign_text <- function(grob) {
  if (inherits(grob, "text")) return(as.character(grob$label))
  unlist(lapply(c(grob$grobs, as.list(grob$children)), ggalign_text),
    use.names = FALSE)
}

test_that("ggalign grobs retain their panels and nested insets at locus anchors", {
  skip_if_not_installed("ggalign", "1.2.0")
  values <- data.frame(x = 1:3, y = c(3, 1, 2))
  first <- ggplot2::ggplot(values, ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::labs(title = "left-child")
  second <- first + ggplot2::labs(title = "right-child")
  children <- list(
    ggalign::align_plots(first, second),
    ggalign::stack_discretev(matrix(1:6, 3)) +
      ggalign::ggfree(first) + ggalign::ggfree(second),
    first + ggalign::inset(second,
      vp = grid::viewport(.7, .7, width = .5, height = .5))
  )
  host <- ggideogram(ggalign_karyotype[1, ],
    tracks = list(child = geom_track(width = 35)))
  for (child in children) {
    child <- ggalign::ggalignGrob(child)
    anchors <- data.frame(Chr = "A", Pos = 50)
    anchors$Plot <- list(child)
    direct <- geom_locus_inset(child, data = anchors,
      mapping = ggplot2::aes(chr = Chr, position = Pos), track = "child",
      width = grid::unit(.65, "npc"), height = grid::unit(.7, "npc"))
    mapped <- geom_locus_inset(data = anchors,
      mapping = ggplot2::aes(chr = Chr, position = Pos, plot = Plot), track = "child",
      width = grid::unit(.65, "npc"), height = grid::unit(.7, "npc"))
    for (component in list(direct, mapped)) {
      plot <- host + component
      grob <- ggplot2::ggplotGrob(plot)
      expect_true(all(c("left-child", "right-child") %in% ggalign_text(grob)))
      built <- ggplot2::ggplot_build(plot)
      expect_equal(tail(built$data, 1)[[1]]$ideogram_position, 50)
    }
  }
  heatmap <- suppressMessages(ggalign::ggheatmap(matrix(1:9, 3))) +
    ggplot2::scale_fill_continuous(name = "child-signal")
  plot <- host + geom_locus_inset(suppressMessages(ggalign::ggalignGrob(heatmap)),
    data = data.frame(Chr = "A", Pos = 50),
    mapping = ggplot2::aes(chr = Chr, position = Pos), track = "child",
    width = grid::unit(.65, "npc"), height = grid::unit(.7, "npc"))
  expect_true("child-signal" %in% ggalign_text(ggplot2::ggplotGrob(plot)))
  expect_length(ggplot2::ggplot_build(plot)$plot$guides$aesthetics, 0L)
})

test_that("ggalign compositions render through the standard save boundary", {
  skip_if_not_installed("ggalign", "1.2.0")
  skip_if_not_installed("ragg")
  marks <- data.frame(Chr = c("A", "B"), Pos = c(30, 40), Type = c("one", "two"))
  plot <- ggideogram(ggalign_karyotype) +
    geom_locus(data = marks, ggplot2::aes(chr = Chr, position = Pos, colour = Type), size = 2)
  circle <- ggideogram(ggalign_karyotype, orientation = "circular")
  children <- list(
    ggalign::align_plots(plot, circle, guides = "r"),
    ggalign::stack_discretev(matrix(1:6, 3)) + ggalign::ggfree(plot)
  )
  files <- vapply(children, function(x) tempfile(fileext = ".png"), character(1))
  on.exit(unlink(files), add = TRUE)
  for (i in seq_along(children)) {
    expect_no_warning(ggplot2::ggsave(files[i], children[[i]],
      width = 180, height = 110, units = "mm", device = ragg::agg_png, dpi = 72))
    expect_gt(file.info(files[i])$size, 0)
  }
  outward <- plot + ggalign::inset(circle,
    vp = grid::viewport(.75, .75, width = .4, height = .4))
  expect_s3_class(ggalign::ggalignGrob(outward), "gtable")
})
