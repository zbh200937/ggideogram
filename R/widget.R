#' View an ideogram as an interactive HTML widget
#'
#' Interactive layers can use `ggiraph::geom_point_interactive()` through
#' [geom_track()], or other compatible geoms with the public projections.
#' Map `tooltip` and `data_id` explicitly in those layers. The original ggplot
#' remains suitable for [ggplot2::ggsave()].
#'
#' Optional `data` connects each unique ID to its source chromosome interval.
#' Clicking a matching mark selects the record; the widget can filter by
#' category, open record links and download selected rows as UTF-8 TSV.
#' All supplied columns are retained in the download. Selection uses source
#' coordinates, including the full original interval in a local view.
#' The declared physical SVG size is retained in narrow containers by scrolling.
#'
#' @param plot A [ggideogram()] plot with compatible interactive layers.
#' @param data Optional data frame of original records. `id` matches the
#'   interactive layer's `data_id`; identifiers must be unique.
#' @param id,chr,start,end Column names in `data`. For points, `end = NULL`
#'   uses the start position for both endpoints. Multi-genome `chr` must contain
#'   the same [chr_key()] values mapped in the plot.
#' @param category Optional category column for a display filter.
#' @param url Optional column of HTTP(S) record links.
#' @param width,height SVG output dimensions in millimetres.
#'
#' @return A `girafe` htmlwidget. Save with [htmlwidgets::saveWidget()].
#' @export
as_ideogram_widget <- function(plot, data = NULL, id = "ID", chr = "Chr",
                               start = "Start", end = "End", category = NULL,
                               url = NULL, width = 185, height = 90) {
  for (package in c("ggiraph", "htmlwidgets")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stopf("Install the optional package `%s` to create an interactive ideogram.", package)
    }
  }
  check_positive_layout(width, "width")
  check_positive_layout(height, "height")
  layout <- ideogram_plot_layout(plot)
  payload <- list(width = width, height = height, records = list(), fields = character(),
                  filter = !is.null(category))
  if (!is.null(data)) {
    check_projection_data(data)
    ids <- validate_chr(projection_column(data, id, "id"), "id")
    if (anyDuplicated(ids)) stopf("Widget record IDs must be unique.")
    chromosomes <- validate_chr(projection_column(data, chr, "chr"), "chr")
    starts <- validate_coordinate(projection_column(data, start, "start"), "start")
    ends <- if (is.null(end)) starts else
      validate_coordinate(projection_column(data, end, "end"), "end")
    source <- layout$data$karyotype
    index <- match(chromosomes, source$.chr)
    if (anyNA(index)) stopf("Widget records contain unknown source chromosomes.")
    if (any(starts < source$.start[index] | ends > source$.end[index] | starts > ends)) {
      stopf("Widget records require start <= end within the source chromosome bounds.")
    }
    categories <- if (is.null(category)) rep("", nrow(data)) else
      as.character(projection_column(data, category, "category"))
    categories[is.na(categories)] <- "\u672a\u5206\u7c7b"
    urls <- if (is.null(url)) rep("", nrow(data)) else
      as.character(projection_column(data, url, "url"))
    urls[is.na(urls)] <- ""
    if (any(nzchar(urls) & !grepl("^https?://", urls))) {
      stopf("Record links must be HTTP(S) URLs or empty.")
    }
    payload$fields <- names(data)
    payload$records <- lapply(seq_len(nrow(data)), function(i) {
      list(id = ids[i], chr = (source$.display_name %||% source$.chr)[index[i]],
           start = starts[i], end = ends[i],
           category = categories[i], url = urls[i],
           values = unname(lapply(data, function(column) column[i])))
    })
  }
  widget <- ggiraph::girafe(ggobj = plot, width_svg = width / 25.4,
    height_svg = height / 25.4, options = list(
      ggiraph::opts_sizing(rescale = FALSE),
      ggiraph::opts_selection(type = "none"),
      ggiraph::opts_hover(css = "opacity:0.75;"),
      ggiraph::opts_tooltip(css = paste0("background-color:white;color:#202020;",
        "border:1px solid #bac4cd;padding:5px 7px;font-size:11px;line-height:1.45;")),
      ggiraph::opts_toolbar(saveaspng = FALSE, hidden = "fullscreen")))
  script <- system.file("htmlwidgets", "ideogram-controls.js", package = "ggideogram")
  htmlwidgets::onRender(widget, paste(readLines(script, warn = FALSE), collapse = "\n"),
                        data = payload)
}
