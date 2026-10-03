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
      if (inherits(chr_layout, "ideogram_layout_v2")) {
        data <- transform_ideogram_semantics(chr_layout, data)
        data$ideogram_projected <- TRUE
      }
      settings <- position$setup_params(data)
      data <- position$setup_data(data, settings)
      data <- position$compute_layer(data, settings, layout)
      if (inherits(chr_layout, "ideogram_layout_v2") &&
          "ideogram_link" %in% names(data) && any(data$ideogram_link %in% TRUE)) {
        anchor <- transform_ideogram_semantics(chr_layout, data)
        data$x <- anchor$x
        data$y <- anchor$y
      }
      data
    })
}
