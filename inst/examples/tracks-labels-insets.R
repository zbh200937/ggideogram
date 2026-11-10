# Run from the source checkout:
# Rscript -e 'pkgload::load_all("."); source("inst/examples/tracks-labels-insets.R")'
library(ggideogram)
library(ggplot2)

k <- human_karyotype[1:2, ]
d <- gene_density[gene_density$Chr %in% k$Chr, ]
d$Mid <- (d$Start - 1 + d$End) / 2
rna <- Random_RNAs_500[Random_RNAs_500$Chr %in% k$Chr, ]
rna <- do.call(rbind, lapply(split(rna, rna$Chr), function(x) head(x, 3)))
rna$Mid <- (rna$Start - 1 + rna$End) / 2
children <- data.frame(Chr = k$Chr, Pos = k$End / 2)
children$Plot <- lapply(k$Chr, function(chr) {
  ggplot(d[d$Chr == chr, ], aes(Value)) +
    geom_histogram(bins = 12, fill = '#4477AA', colour = 'white', linewidth = .2) +
    labs(x = 'Genes / bin', y = 'Bins') + theme_minimal(base_size = 6)
})

make_plot <- function(orientation) {
  circular <- orientation == 'circular'
  ggideogram(k, orientation = orientation, ncol = if (circular) NULL else 2,
    radius = if (circular) 60 else NULL,
    opening_angle = if (circular) 60 else NULL,
    tracks = list(
      density = geom_track(side = 'left', width = 3, data = d,
        aes(chr = Chr, x = Mid, y = Value),
        layers = list(geom_line(linewidth = .25, colour = '#4477AA'),
          geom_point(size = .35, colour = '#4477AA')),
        label = 'Genes / bin', axis = TRUE),
      labels = geom_track(width = if (circular) 16 else 12),
      distribution = geom_track(side = if (circular) 'left' else 'right',
        width = if (circular) 40 else 20))) +
    geom_locus(data = rna, aes(chr = Chr, position = Mid, label = Type),
      geom = 'text', position = 'spread', track = 'labels', size = 2.5) +
    geom_locus_inset(data = children, aes(chr = Chr, position = Pos, plot = Plot),
      track = 'distribution', width = grid::unit(28, 'mm'), height = grid::unit(22, 'mm'))
}

out <- Sys.getenv('GGIDEOGRAM_EXAMPLE_OUTPUT', 'work/tracks-labels-insets')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (orientation in c('horizontal', 'vertical', 'circular')) {
  p <- make_plot(orientation)
  for (size in c(160, 220)) {
    ggsave(file.path(out, paste0(orientation, '-', size, '.png')),
      p, width = size, height = size, units = 'mm', dpi = 180, bg = 'white')
  }
}
