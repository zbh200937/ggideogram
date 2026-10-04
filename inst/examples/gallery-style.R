# Shared physical text and legend settings for the example galleries.
gallery_axis <- function(...) list(...)
gallery_legend_spacing <- ggideogram::theme_ideogram()$legend.box.spacing

gallery_grid <- function(karyotype, tracks, width_mm, height_mm,
    orientations = 'vertical', ...) {
  semantic <- ggideogram::as_ideogram_data(karyotype)
  n <- nrow(semantic$karyotype)
  columns <- seq_len(n)
  candidates <- expand.grid(ncol = columns[n %% columns == 0],
    orientation = orientations, stringsAsFactors = FALSE)
  aspects <- vapply(seq_len(nrow(candidates)), function(i) {
    axes <- semantic$karyotype$.chr[seq.int(1L, n, by = candidates$ncol[i])]
    plot <- ggideogram::ggideogram(semantic, ncol = candidates$ncol[i],
      orientation = candidates$orientation[i], tracks = tracks, axis = axes, ...)
    diff(plot$coordinates$limits$x) / diff(plot$coordinates$limits$y)
  }, numeric(1))
  best <- which.min(abs(log(aspects / (width_mm / height_mm))))
  list(ncol = candidates$ncol[best], orientation = candidates$orientation[best],
    aspect = aspects[best])
}

gallery_theme <- function(...) {
  ggplot2::theme(
    text = ggplot2::element_text(size = 8.25, colour = 'black'),
    plot.title = ggplot2::element_text(size = 10, face = 'plain', hjust = 0),
    plot.subtitle = ggplot2::element_text(size = 8.25),
    plot.caption = ggplot2::element_text(size = 7.5, hjust = 0),
    legend.position = 'bottom',
    legend.title = ggplot2::element_text(size = 8.25),
    legend.text = ggplot2::element_text(size = 8.25),
    legend.key.height = grid::unit(2.2, 'mm'),
    legend.key.width = grid::unit(4.2, 'mm'),
    legend.spacing.x = grid::unit(2.5, 'mm'),
    legend.box.spacing = gallery_legend_spacing,
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    plot.margin = ggplot2::margin(8, 16, 6, 24)
  ) + ggplot2::theme(...)
}

gallery_key <- function(data, params, size) {
  key <- ggplot2::draw_key_polygon(data, params, size)
  grid::editGrob(key, width = grid::unit(3.5, 'mm'), height = grid::unit(1.8, 'mm'))
}

gallery_save <- function(plot, directory, name, width = 185, height) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(directory, paste0(name, '.png')), plot,
    width = width, height = height, units = 'mm', dpi = 220,
    device = ragg::agg_png, bg = 'white')
  ggplot2::ggsave(file.path(directory, paste0(name, '.pdf')), plot,
    width = width, height = height, units = 'mm',
    device = grDevices::cairo_pdf, bg = 'white')
}
