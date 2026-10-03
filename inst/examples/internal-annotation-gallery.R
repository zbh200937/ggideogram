# Run: Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/internal-annotation-gallery.R")'
# References: karyoploteR r0/r1 lanes and kpHeatmap interval semantics.
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))
data(human_karyotype, package = 'ggideogram')
data(gene_density, package = 'ggideogram')
data(LTR_density, package = 'ggideogram')
output <- 'work/api-optimization/examples/internal-annotation-gallery'

k <- human_karyotype[1:4, ]
genes <- gene_density[gene_density$Chr %in% k$Chr, ]
ltr <- LTR_density[LTR_density$Chr %in% k$Chr, ]
# Normalize the real final partial window as well as the full 1 Mb windows.
genes$Rate <- genes$Value * 1e6 / (genes$End - genes$Start + 1)
ltr$Rate <- ltr$Value * 1e6 / (ltr$End - ltr$Start + 1)
high <- ceiling(max(c(genes$Rate, ltr$Rate)) / 100) * 100
fill_scale <- function() scale_fill_gradient(low = '#F2F5F8', high = '#4477AA',
  limits = c(0, high), name = 'Features / Mb',
  guide = guide_colourbar(barwidth = grid::unit(32, 'mm'), barheight = grid::unit(2, 'mm')))
inside <- track_layout(genes = geom_track(side = "overlay", width = 1.16, offset = -0.62, limits = c(0,
    1)), ltr = geom_track(side = "overlay", width = 1.16, offset = 0.62, limits = c(0,
    1)))

# 1. Two internal variables with shared count units and a visible center gutter.
heatmap <- ggideogram(k, ncol = 4, max_chr_length = 22, chromosome_width = 2.4,
  chromosome_gap = 3.4, tracks = inside,
  axis = '1', axis_units = 'Mb', axis_colour = 'black', axis_size = 2.7,
  axis_linewidth = 0.22, axis_tick_length = 1.2,
  name_colour = 'black', name_size = 3, fill = '#F0F2F4',
  colour = '#626A73', linewidth = 0.22, padding = 1.2) +
  geom_chr(component = "fill", data = genes, track = "genes", aes(chr = Chr, start = Start,
      end = End, fill = Rate)) +
  geom_chr(component = "fill", data = ltr, track = "ltr", aes(chr = Chr, start = Start,
      end = End, fill = Rate)) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  geom_track(data = transform(k, Pos = End, Lane = 0.5), track = "genes", mapping = aes(chr = Chr,
      x = Pos, y = Lane), geom = ggplot2::geom_text(label = "G", size = 2.5, vjust = 1.7)) +
  geom_track(data = transform(k, Pos = End, Lane = 0.5), track = "ltr", mapping = aes(chr = Chr,
      x = Pos, y = Lane), geom = ggplot2::geom_text(label = "L", size = 2.5, vjust = 1.7)) +
  fill_scale() + gallery_theme()
gallery_save(heatmap, output, '01-internal-heatmap', height = 130)

# 2. TAIR10 / Araport11 gene spans and explicitly annotated mRNA starts.
f <- read_chr_features(system.file('extdata', 'arabidopsis-first-genes.gff3',
  package = 'ggideogram'))
regions <- f[f$Type == 'gene', ]
regions$Label <- sub('gene:', '', regions$ID)
regions$Midpoint <- (regions$Start - 1 + regions$End) / 2
starts <- f[f$Type == 'mRNA', ]
starts$Position <- ifelse(starts$Strand == '+', starts$Start, starts$End)
starts <- unique(starts[c('Chr', 'Position')])
regions$Lane <- starts$Lane <- 0.5
ar <- data.frame(Chr = '1', Start = 0, End = 30427671)
bands <- ggideogram(chr_view(ar, '1', 2500, 15000), orientation = 'horizontal',
  max_chr_length = 40, chromosome_width = 1.4,
  tracks = track_layout(body = geom_track(side = "overlay", width = 1.4, limits = c(0, 1))),
  axis = TRUE, axis_side = 'left', axis_units = 'kb', axis_breaks = seq(3000, 15000, 3000),
  axis_colour = 'black', axis_size = 2.7, axis_linewidth = 0.22, axis_tick_length = 1.2,
  name_colour = 'black', name_size = 3, fill = '#F0F2F4',
  colour = '#626A73', linewidth = 0.22, padding = 0.8) +
  geom_chr(component = "fill", data = regions, track = "body", aes(chr = Chr, start = Start,
      end = End, fill = Strand), key_glyph = gallery_key) +
  geom_track(data = starts, track = "body", mapping = aes(chr = Chr, x = Position, y = Lane,
      shape = "Transcript start"), geom = ggplot2::geom_point(size = 1.4, fill = "white",
      colour = "#202020", stroke = 0.25)) +
  geom_track(data = regions, track = "body", mapping = aes(chr = Chr, x = Midpoint, y = Lane,
      label = Label), geom = ggplot2::geom_text(size = 2.7)) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  scale_fill_manual(values = c('+' = '#AACCEE', '-' = '#F4CBD2'), name = 'Strand') +
  scale_shape_manual(values = c('Transcript start' = 21), name = NULL) +
  gallery_theme(plot.margin = margin(8, 14, 6, 28))
gallery_save(bands, output, '02-regions-and-starts', height = 48)

# 3. The two internal variables plus a native external gene-density line track.
view <- chr_view(human_karyotype, '1', 110e6, 140e6)
line <- genes[genes$Chr == '1', ]
line$Position <- (line$Start - 1 + line$End) / 2
line <- view_chr_data(line, view, position = 'Position')
tracks <- track_layout(genes = inside$genes, ltr = inside$ltr, density = geom_track(side = "right",
    width = 3.5, gap = 0.8, limits = c(0, 150), reverse = TRUE))
combined <- ggideogram(view, orientation = 'horizontal', max_chr_length = 40,
  chromosome_width = 2.4, tracks = tracks,
  axis = TRUE, axis_side = 'left', axis_breaks = seq(110e6, 140e6, 5e6), axis_units = 'Mb',
  axis_colour = 'black', axis_size = 2.7, axis_linewidth = 0.22, axis_tick_length = 1.2,
  name_colour = 'black', name_size = 3, fill = '#F0F2F4',
  colour = '#626A73', linewidth = 0.22, padding = 1.2) +
  geom_chr(component = "fill", data = genes, track = "genes", aes(chr = Chr, start = Start,
      end = End, fill = Rate)) +
  geom_chr(component = "fill", data = ltr, track = "ltr", aes(chr = Chr, start = Start,
      end = End, fill = Rate)) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  geom_track(data = line, track = "density", mapping = aes(chr = Chr, x = Position, y = Rate),
      geom = ggplot2::geom_line(colour = "#4477AA", linewidth = 0.35)) +
  geom_track(track = "density", axis = list(position = "start", breaks = c(0, 75, 150),
      size = 2.5, colour = "black", linewidth = 0.22, tick_length = 1.2)) +
  geom_track(data = data.frame(Chr = "1", Pos = 1.11e+08, Rate = 140), track = "density",
      mapping = aes(chr = Chr, x = Pos, y = Rate), geom = ggplot2::geom_text(label = "Genes / Mb",
          size = 2.5, hjust = 0)) +
  geom_track(data = data.frame(Chr = "1", Pos = 1.4e+08, Lane = 0.5), track = "genes", mapping = aes(chr = Chr,
      x = Pos, y = Lane), geom = ggplot2::geom_text(label = "G", size = 2.5, hjust = -0.3)) +
  geom_track(data = data.frame(Chr = "1", Pos = 1.4e+08, Lane = 0.5), track = "ltr", mapping = aes(chr = Chr,
      x = Pos, y = Lane), geom = ggplot2::geom_text(label = "L", size = 2.5, hjust = -0.3)) +
  fill_scale() + gallery_theme(plot.margin = margin(8, 18, 6, 28))
gallery_save(combined, output, '03-internal-and-external', height = 75)
writeLines(c(
  '01: Human Chr1–4; G = genes, L = LTR elements. Shared count-per-Mb fill scale.',
  '02: Arabidopsis Chr1 2.5–15 kb; gene intervals coloured by strand.',
  'Circles: distinct annotated mRNA starts (Start on +, End on -).',
  '03: Human Chr1 110–140 Mb; G above L, external line = genes per Mb.',
  'Density values: bundled RIdeogram gene_density and LTR_density;',
  'Rate = original count divided by the actual window width, multiplied by 1e6.',
  'Plant annotation: TAIR10 / Araport11, Ensembl Plants release 62.',
  'Reference: https://bernatgel.github.io/karyoploter_tutorial/Tutorial/DataPositioning/DataPositioning.html'
), file.path(output, 'captions.txt'))
