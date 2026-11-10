# Run from the source checkout:
# Rscript -e 'pkgload::load_all("."); source("inst/examples/component-spacing.R")'
library(ggideogram)
library(ggplot2)

# The bundled human gene counts are measured in their original 1 Mb windows.
k <- human_karyotype[1:3, ]
d <- gene_density[gene_density$Chr %in% k$Chr, ]
d$Mid <- (d$Start - 1 + d$End) / 2
make_plot <- function(orientation) ggideogram(k, orientation = orientation,
  opening_angle = if (orientation == "circular") 60 else NULL,
  tracks = list(density = geom_track(data = d,
    aes(chr = Chr, x = Mid, y = Value), geom = geom_line(linewidth = 0.3),
    label = "Genes / bin", axis = list(
      position = if (orientation == "circular") "auto" else "start", n = 2))))

out <- Sys.getenv("GGIDEOGRAM_EXAMPLE_OUTPUT", "work/component-spacing")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, width, height) {
  ggsave(file.path(out, paste0(name, ".png")), p, width = width, height = height,
    units = "mm", dpi = 180, bg = "white")
}

for (orientation in c("horizontal", "vertical", "circular")) {
  p <- make_plot(orientation)
  save_plot(p, paste0(orientation, "-180"), 180, 140)
  save_plot(p, paste0(orientation, "-120"), 120, 100)
}

# Removing the title and axis releases their automatic text space.
clean <- make_plot("horizontal") +
  geom_track(track = "density", label = NULL, axis = FALSE)
save_plot(clean, "track-text-removed", 180, 100)

# Original randomly sampled RNA annotations illustrate point placement only.
rna <- Random_RNAs_500[Random_RNAs_500$Chr %in% k$Chr, ]
rna$Mid <- (rna$Start - 1 + rna$End) / 2
markers <- ggideogram(k) + geom_locus(data = rna,
  aes(chr = Chr, position = Mid, colour = Type), gap = 3, size = 1.5)
save_plot(markers, "separate-marker-lanes", 180, 130)
