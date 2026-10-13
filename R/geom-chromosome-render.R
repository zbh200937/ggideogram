# Chromosome caps are resolved in the actual panel viewport, including insets.
chromosome_cap_ratio <- function(layout, panel_params) {
  if (is_circular_layout(layout)) return(1)
  mm <- c(grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE),
    grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE))
  scale <- mm / c(diff(panel_params$x.range), diff(panel_params$y.range))
  if (layout$orientation == "horizontal") scale[2] / scale[1] else scale[1] / scale[2]
}

chromosome_render_data <- function(data, layout, cap_ratio, curve_points = layout$curve_points) {
  complete <- !(layout$chrom$.cut_start %||% rep(FALSE, nrow(layout$chrom))) &
    !(layout$chrom$.cut_end %||% rep(FALSE, nrow(layout$chrom)))
  width <- layout$chromosome_width
  if (any(complete)) width <- min(width, min(layout$chrom$.display_length[complete]) / cap_ratio)
  pieces <- lapply(split(seq_len(nrow(data)), data$group), function(rows) {
    style <- data[rows[1], , drop = FALSE]
    g <- layout$chrom[layout$index[[as.character(style$ideogram_chr)]], , drop = FALSE]
    start <- style$ideogram_start
    end <- style$ideogram_end
    if (is.null(start) || is.na(start)) start <- g$.start
    if (is.null(end) || is.na(end)) end <- g$.end
    polygon <- chromosome_section_polygon(g, width,
      curve_points, start, end, cap_ratio = cap_ratio)
    result <- style[rep(1L, nrow(polygon)), , drop = FALSE]
    result$x <- polygon$x
    result$y <- polygon$y
    result
  })
  do.call(rbind, pieces)
}

chromosome_render_grob <- function(data, panel_params, coord, geom, params,
                                   curve_points = coord$layout$curve_points) {
  if (!nrow(data)) return(grid::nullGrob())
  if (is_circular_layout(coord$layout)) {
    return(do.call(geom$draw_panel, c(list(data, panel_params, coord), params)))
  }
  grid::gTree(data = data, panel_params = panel_params, coord = coord,
    geom = geom, params = params, curve_points = curve_points, cl = "chr_body_render")
}

GeomChromosomePolygon <- ggplot2::ggproto("GeomChromosomePolygon", ggplot2::GeomPolygon,
  required_aes = c("x", "y", "ideogram_chr"),
  default_aes = utils::modifyList(ggplot2::GeomPolygon$default_aes,
    ggplot2::aes(ideogram_start = NA_real_, ideogram_end = NA_real_)),
  draw_panel = function(data, panel_params, coord, rule = "evenodd",
      lineend = "butt", linejoin = "round", linemitre = 10, na.rm = FALSE,
      curve_points = 32) {
    chromosome_render_grob(data, panel_params, coord, ggplot2::GeomPolygon,
      list(rule = rule, lineend = lineend, linejoin = linejoin, linemitre = linemitre), curve_points)
  })

GeomChromosomePath <- ggplot2::ggproto("GeomChromosomePath", ggplot2::GeomPath,
  required_aes = c("x", "y", "ideogram_chr"),
  draw_panel = function(data, panel_params, coord, lineend = "butt",
      linejoin = "round", linemitre = 10, na.rm = FALSE, curve_points = 32) {
    chromosome_render_grob(data, panel_params, coord, ggplot2::GeomPath,
      list(lineend = lineend, linejoin = linejoin, linemitre = linemitre), curve_points)
  })

#' @method makeContent chr_body_render
#' @importFrom grid makeContent
#' @export
makeContent.chr_body_render <- function(x) {
  ratio <- chromosome_cap_ratio(x$coord$layout, x$panel_params)
  data <- chromosome_render_data(x$data, x$coord$layout, ratio, x$curve_points)
  body <- do.call(x$geom$draw_panel, c(list(data, x$panel_params, x$coord), x$params))
  grid::setChildren(x, grid::gList(body))
}
