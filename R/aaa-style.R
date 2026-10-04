# Shared visual defaults for chromosome components.

ideogram_text_size <- function(role) {
  unname(c(chromosome = 3.2, label = 3, track = 2.8, bp = 2.8, value = 2.2)[role])
}

ideogram_colour <- function(role) {
  unname(c(body = "#F7F7F7", outline = "#4D4D4D", text = "#202020",
    axis = "#666666", annotation = "#626A73", gene = "#4477AA")[role])
}

ideogram_linewidth <- function(role) {
  unname(c(body = 0.4, axis = 0.3, annotation = 0.22)[role])
}
