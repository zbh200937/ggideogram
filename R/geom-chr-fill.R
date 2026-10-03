#' @noRd
geom_chr_fill <- function(mapping = NULL, data = NULL, track = NULL,
    stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = FALSE) {
  if (!identical(stat, "identity")) stopf("Chromosome fills require stat = 'identity'.")
  structure(list(mapping = mapping, data = data, track = track,
    position = position, params = list(...), na.rm = na.rm,
    show.legend = show.legend, inherit.aes = inherit.aes),
    class = "ggideogram_fill_component")
}

GeomIdeogramRect <- ggplot2::ggproto("GeomIdeogramRect", ggplot2::GeomRect,
  required_aes = c("xmin", "xmax", "ymin", "ymax", "ideogram_chr"),
  draw_panel = function(data, panel_params, coord, linejoin = "mitre",
      lineend = "butt", na.rm = FALSE) {
    if (!nrow(data)) return(grid::nullGrob())
    mask_chr_geometry(data, panel_params, coord, function(rows) {
      ggplot2::GeomRect$draw_panel(rows, panel_params, coord,
        linejoin = linejoin, lineend = lineend)
    })
  })

#' @method ggplot_add ggideogram_fill_component
#' @export
ggplot_add.ggideogram_fill_component <- function(object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  spec <- track_table_row(layout, object$track)
  if (spec$side != "overlay") stopf("`geom_chr_fill()` requires an overlay track.")
  data <- object$data
  if (!is.data.frame(data)) stopf("Chromosome fills require data-frame `data`.")
  fields <- mapped_fields(data, object$mapping, c("chr", "start", "end"), what = "mapping")
  chr <- validate_chr(fields$chr, "chr")
  first <- validate_coordinate(fields$start, "start")
  last <- validate_coordinate(fields$end, "end")
  if (any(first < 1 | first > last | first != floor(first) | last != floor(last))) {
    stopf("Region annotations require positive integer start <= end.")
  }
  source <- layout$data$karyotype %||% layout$chrom
  index <- match(chr, source$.chr)
  if (anyNA(index) || any(first - 1 < source$.start[index] | last > source$.end[index])) {
    stopf("Region annotations are outside the source chromosome bounds.")
  }
  data$.source_start <- data$.source_start %||% first
  data$.source_end <- data$.source_end %||% last
  visible_index <- match(chr, layout$chrom$.chr)
  lo <- pmax(first - 1, layout$chrom$.start[visible_index])
  hi <- pmin(last, layout$chrom$.end[visible_index])
  keep <- !is.na(visible_index) & lo < hi
  if (!any(keep)) return(plot)
  data <- data[keep, , drop = FALSE]
  chr <- chr[keep]; lo <- lo[keep]; hi <- hi[keep]
  a <- offset_chr_points_signed(layout, chr,
    project_positions_checked(layout, chr, lo), rep(spec$low_offset, length(chr)))
  b <- offset_chr_points_signed(layout, chr,
    project_positions_checked(layout, chr, hi), rep(spec$high_offset, length(chr)))
  data$.xmin <- pmin(a$x, b$x); data$.xmax <- pmax(a$x, b$x)
  data$.ymin <- pmin(a$y, b$y); data$.ymax <- pmax(a$y, b$y)
  data$ideogram_chr <- chr
  mapping <- track_mapping_without(object$mapping, c("chr", "start", "end"))
  mapping <- combine_track_mapping(mapping, ggplot2::aes(
    xmin = .data$.xmin, xmax = .data$.xmax, ymin = .data$.ymin, ymax = .data$.ymax,
    ideogram_chr = .data$ideogram_chr))
  params <- utils::modifyList(list(colour = NA), object$params)
  if ("colour" %in% names(mapping) && !"colour" %in% names(object$params)) params$colour <- NULL
  key <- params$key_glyph
  params$key_glyph <- NULL
  plot + do.call(ggplot2::layer, c(list(data = data, mapping = mapping,
    stat = "identity", geom = GeomIdeogramRect, position = object$position,
    show.legend = object$show.legend, inherit.aes = object$inherit.aes, key_glyph = key),
    list(params = c(list(na.rm = object$na.rm), params))))
}
