# Shared physical text and legend settings for the example galleries.
gallery_theme <- function(...) {
  ggplot2::theme(
    text = ggplot2::element_text(size = 8.25, colour = 'black'),
    legend.position = 'bottom',
    legend.title = ggplot2::element_text(size = 8.25),
    legend.text = ggplot2::element_text(size = 8.25),
    legend.key.height = grid::unit(2.2, 'mm'),
    legend.key.width = grid::unit(4.2, 'mm'),
    legend.spacing.x = grid::unit(2.5, 'mm'),
    legend.box.spacing = grid::unit(1.3, 'mm'),
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
