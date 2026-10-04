# Run: micromamba run -n multiomics Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/label-layout-gallery.R")'
# Arabidopsis TAIR10 / Araport11 gene spans, Ensembl Plants release 62.
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))
g <- read_chr_features(system.file('extdata', 'arabidopsis-gene-cluster.gff3', package = 'ggideogram'))
g$Gene <- sub('gene:', '', g$ID)
g$Position <- (g$Start - 1 + g$End) / 2
g$Label <- ifelse(!is.na(g$Name) & g$Name != g$Gene,
  paste(g$Gene, g$Name, sep = ' · '), g$Gene)
k <- data.frame(Chr = '1', Start = 0, End = 30427671)
output <- Sys.getenv('GGIDEOGRAM_EXAMPLE_OUTPUT', 'work/label-layout-gallery')

label_base <- function(view, tracks, orientation, breaks) {
  ggideogram(view, orientation = orientation, max_chr_length = 40,
    chromosome_width = 0.9, tracks = tracks, axis = TRUE, axis_side = 'left',
    axis_breaks = breaks, axis_units = 'kb', axis_colour = 'black', axis_size = 2.7,
    axis_linewidth = 0.22, axis_tick_length = 1.2,
    name_colour = 'black', name_size = 3, fill = '#F0F2F4',
    colour = '#626A73', linewidth = 0.22, padding = 1.2) + gallery_theme()
}

# 1. Sparse gene names, with true midpoint loci and fixed-size native points.
sparse <- g[1:3, ]
sparse$Label <- sparse$Name
view <- chr_view(k, '1', 2500, 15000)
p <- label_base(view, track_layout(labels = geom_track(side = "right", width = 9, gap = 0.7)),
  'horizontal', seq(3000, 15000, 3000)) +
  geom_locus(geom = "text", position = "spread", data = sparse, track = "labels", aes(chr = Chr, position = Position,
      label = Label)) +
  geom_locus(geom = "point", data = sparse, aes(chr = Chr, position = Position), gap = 0,
      side = "right", size = 1.3, colour = "#4477AA")
p <- p + labs(title = 'Arabidopsis · Chr1 2.5–15 kb',
  caption = 'TAIR10 / Araport11. Three gene names at their source midpoints.')
gallery_save(p, output, '01-sparse-names', height = 55)

# 2. A real 100 kb gene cluster, with all label edges on one straight baseline.
view <- chr_view(k, '1', 0, 100000)
tracks <- track_layout(labels = geom_track(side = "right", width = 28, gap = 0.7))
dense <- label_base(view, tracks, 'vertical', seq(0, 100000, 25000)) +
  geom_locus(geom = "text", position = "spread", data = g, track = "labels", aes(chr = Chr,
      position = Position, label = Label)) +
  geom_locus(geom = "point", data = g, aes(chr = Chr, position = Position), gap = 0, side = "right",
      size = 1.3, colour = "#4477AA")
dense <- dense + labs(title = 'Arabidopsis · Chr1 0–100 kb',
  caption = 'TAIR10 / Araport11; 24 source genes.')
gallery_save(dense, output, '02-gene-cluster-wide', width = 125, height = 135)

# 3. The identical input and plot in a narrower, shorter publication slot.
gallery_save(dense, output, '03-gene-cluster-narrow', width = 85, height = 100)
writeLines(c(
  '01: Three neighbouring genes, labels anchored at annotated gene midpoints.',
  '02–03: First 24 annotated genes of Arabidopsis Chr1, shown in 0–100 kb.',
  'Gene IDs and distinct annotated symbols share one label line.',
  'All labels start on the same straight baseline; leaders retain true source loci.',
  'Wide: 125 × 135 mm; narrow: 85 × 100 mm; text size 3 mm in both.',
  'Source: TAIR10 / Araport11, Ensembl Plants release 62; GFF subset packaged with ggideogram.',
  'Reference: https://jokergoo.github.io/circlize/reference/circos.genomicLabels.html'
), file.path(output, 'captions.txt'))
