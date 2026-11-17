# Native Positions work on displayed locus coordinates, while the semantic
# fields remain the source of the genomic anchor and later layout rebuilds.
locus_position <- function(position) {
  if (inherits(position, "PositionIdentity") ||
      inherits(position, "PositionChrRepel") ||
      inherits(position, "PositionChrLocus")) return(position)
  ggplot2::ggproto("PositionChrLocus", position,
    setup_params = function(data) list(),
    setup_data = function(data, params) data,
    compute_layer = function(self, data, params, layout) {
      if (!nrow(data)) return(data)
      chr_layout <- layout$coord$layout
      is_link <- inherits(chr_layout, "ideogram_layout_v2") &&
        "ideogram_link" %in% names(data) && any(data$ideogram_link %in% TRUE)
      if (inherits(chr_layout, "ideogram_layout_v2")) {
        data <- transform_ideogram_semantics(chr_layout, data)
        data$ideogram_projected <- TRUE
      }
      if (is_link && all(c("xend", "yend") %in% names(data))) {
        data$.ideogram_anchor_x <- data$x
        data$.ideogram_anchor_y <- data$y
        data$x <- data$xend
        data$y <- data$yend
        data$xend <- NULL
        data$yend <- NULL
      }
      settings <- position$setup_params(data)
      data <- position$setup_data(data, settings)
      data <- position$compute_layer(data, settings, layout)
      if (is_link && nrow(data)) {
        if (all(c(".ideogram_anchor_x", ".ideogram_anchor_y") %in% names(data))) {
          data$xend <- data$x
          data$yend <- data$y
          data$x <- data$.ideogram_anchor_x
          data$y <- data$.ideogram_anchor_y
          data[c(".ideogram_anchor_x", ".ideogram_anchor_y")] <- NULL
        } else {
          anchor <- transform_ideogram_semantics(chr_layout, data)
          data$x <- anchor$x
          data$y <- anchor$y
        }
      }
      data
    })
}
