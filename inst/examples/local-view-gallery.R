# Local genomic views from the bundled RIdeogram human data.
# Run from a checkout:
# Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/local-view-gallery.R")'
# Outputs: work/api-optimization/examples/local-view-gallery/{01-local,02-overview,03-windows}.{png,pdf}
# The density tables contain counts in 1 Mb windows. Their observed midpoints
# are selected without interpolation; source intervals and values are unchanged.

library(ggideogram)
library(ggplot2)
library(patchwork)
source(system.file("examples", "gallery-style.R", package = "ggideogram"))

data(human_karyotype, package = "ggideogram")
data(gene_density, package = "ggideogram")
data(LTR_density, package = "ggideogram")
gene_density$Position <- (gene_density$Start - 1 + gene_density$End) / 2
LTR_density$Position <- (LTR_density$Start - 1 + LTR_density$End) / 2
colours <- c(Genes = "#4477AA", `LTR elements` = "#EE6677")
output <- file.path("work", "api-optimization", "examples", "local-view-gallery")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

local_plot <- function(chr, from, to) {
  view <- chr_view(human_karyotype, chr, from, to)
  genes <- view_chr_data(gene_density, view, position = "Position")
  ltr <- view_chr_data(LTR_density, view, position = "Position")
  p <- ggideogram(view, orientation = "horizontal",
    max_chr_length = 36, chromosome_width = 0.9,
    tracks = track_layout(genes = geom_track(side = "right", width = 3.5, gap = 0.7, limits = c(0,
        150), reverse = TRUE), ltr = geom_track(side = "right", width = 3.5, gap = 0.9, limits = c(0,
        650), reverse = TRUE)),
    axis = TRUE, axis_side = "left", axis_breaks = seq(from, to, length.out = 7),
    axis_units = "Mb", axis_colour = "#000000", axis_linewidth = 0.22,
    axis_tick_length = 1.2, axis_size = 2.9,
    name_colour = "#000000", name_size = 2.9, padding = 1.5,
    fill = "#F0F2F4", colour = "#626A73", linewidth = 0.22) +
    geom_track(data = genes, track = "genes", mapping = aes(chr = Chr, x = Position, y = Value,
        colour = "Genes"), geom = ggplot2::geom_line(linewidth = 0.35)) +
    geom_track(data = ltr, track = "ltr", mapping = aes(chr = Chr, x = Position, y = Value,
        colour = "LTR elements"), geom = ggplot2::geom_line(linewidth = 0.35)) +
    geom_track(track = "genes", axis = list(position = "start", breaks = c(0, 75, 150), size = 2.9,
        colour = "#000000", linewidth = 0.22, tick_length = 1.2)) +
    geom_track(track = "ltr", axis = list(position = "start", breaks = c(0, 300, 600), size = 2.9,
        colour = "#000000", linewidth = 0.22, tick_length = 1.2)) +
    scale_colour_manual(values = colours, name = "Count / 1 Mb window") +
    gallery_theme(legend.direction = "horizontal",
      plot.margin = margin(8, 16, 6, 24))
  p
}

# 1. A local region spanning the annotated Chr1 centromere.
local <- local_plot("1", 110e6, 140e6)
gallery_save(local, output, "01-local", height = 82)

# 2. The same region located in a chromosome overview.
selected <- human_karyotype[1:6, ]
overview <- ggideogram(selected, ncol = 6,
  max_chr_length = 18, chromosome_gap = 2.7,
  axis = "1", axis_units = "Mb", axis_colour = "#000000",
  axis_size = 2.9, axis_linewidth = 0.22, name_colour = "#000000",
  name_size = 2.9, padding = 1.5, fill = "#F0F2F4", colour = "#626A73", linewidth = 0.22) +
  geom_locus(geom = "interval", data = data.frame(Chr = "1", Start = 1.1e+08, End = 1.4e+08),
      aes(chr = Chr, start = Start, end = End, colour = "110–140 Mb"), linewidth = 2) +
  scale_colour_manual(values = c(`110–140 Mb` = "#228833"), name = "Chr1 window") +
  gallery_theme(legend.position = "right", plot.margin = margin(8, 16, 6, 24))
combined <- (overview / local) + plot_layout(heights = c(2, 1))
gallery_save(combined, output, "02-overview", height = 160)

# 3. Two Chr17 windows with the same bp span and identical value limits.
left <- local_plot("17", 0, 18e6)
right <- local_plot("17", 60e6, 78e6)
windows <- (left / right) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
gallery_save(windows, output, "03-windows", height = 152)

writeLines(c(
  "01: Chr1 110–140 Mb; both ends are truncated, centromere retained.",
  "02: Chr1–6 overview locates the same Chr1 window highlighted in green.",
  "03: Chr17 0–18 and 60–78 Mb; equal genomic spans and value limits.",
  "Data: human_karyotype, gene_density, LTR_density bundled with ggideogram,",
  "inherited from RIdeogram: https://github.com/TickingClock1992/RIdeogram .",
  "Density values are original counts per 1 Mb window; no interpolation.",
  paste("R", getRversion(), "; ggplot2", packageVersion("ggplot2"))
), file.path(output, "captions.txt"))
