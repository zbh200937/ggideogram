# Run: micromamba run -n multiomics Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/interactive-gallery.R")'
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))
output <- 'work/api-optimization/examples/interactive-gallery'
k <- data.frame(Chr = '1', Start = 0, End = 30427671)
genes <- function(file) {
  d <- read_chr_features(system.file('extdata', file, package = 'ggideogram'))
  d <- d[d$Type == 'gene', ]
  d$ID <- sub('^gene:', '', d$ID)
  d$Mid <- (d$Start - 1 + d$End) / 2
  d$URL <- paste0('https://plants.ensembl.org/Arabidopsis_thaliana/Gene/Summary?g=', d$ID)
  d$Tip <- sprintf('%s (%s)\nChr%s: %s–%s bp; strand %s',
    d$Name, d$ID, d$Chr, d$Start, d$End, d$Strand)
  d
}
gene_colours <- function() scale_fill_manual(values = c('+' = '#AACCEE', '-' = '#F4CBD2'),
  breaks = c('+', '-'), name = 'Strand')
common <- list(orientation = 'horizontal', chromosome_width = 0.75,
  max_chr_length = 55, fill = '#F0F2F4', colour = '#626A73', linewidth = 0.22,
  name_size = 2.9, name_colour = 'black', axis = TRUE, axis_units = 'kb',
  axis_colour = 'black', axis_size = 2.9, axis_linewidth = 0.22, axis_tick_length = 1.2)
save_interactive <- function(p, metadata, name, height, category, url = NULL) {
  gallery_save(p, output, name, height = height)
  widget <- as_ideogram_widget(p, metadata, category = category, url = url,
    width = 185, height = height)
  htmlwidgets::saveWidget(widget, file.path(output, paste0(name, '.html')),
    selfcontained = FALSE, libdir = 'lib', title = 'ggideogram')
}

# 1. Real gene midpoint markers, with physical labels and gene-page links.
d <- genes('arabidopsis-gene-cluster.gff3')
tracks <- track_layout(body = geom_track(side = "overlay", width = 0.75, limits = c(0, 1)), names = geom_track(side = "right",
    width = 3.5, gap = 0.8))
p <- do.call(ggideogram, c(list(data = chr_view(k, '1', 0, 100000),
  tracks = tracks, axis_breaks = seq(0, 100000, 20000)), common)) +
  geom_track(data = d, track = "body", mapping = aes(chr = Chr, x = Mid, y = 0.5, fill = Strand,
      data_id = ID, tooltip = Tip), geom = ggiraph::geom_point_interactive(shape = 21, size = 1.5,
      stroke = 0.18, colour = "#626A73")) +
  geom_track(data = d[d$Name %in% c("NAC001", "DCL1", "LHY", "KCS1"), ],
    track = "names", aes(chr = Chr, x = Mid, y = .5, label = Name, data_id = ID, tooltip = Tip),
    geom = ggiraph::geom_text_interactive(size = 2.9, show.legend = FALSE)) +
  gene_colours() + gallery_theme()
save_interactive(p, d, '01-gene-markers', height = 45, category = 'Strand', url = 'URL')

# 2. A cropped gene-span view; metadata still contains the complete source genes.
d <- genes('arabidopsis-first-genes.gff3')
view <- chr_view(k, '1', 4500, 13000)
visible <- view_chr_data(d, view)
visible$EdgeStart <- pmax(visible$Start - 1, view$view$start)
visible$Mid <- (visible$EdgeStart + visible$End) / 2
tracks <- track_layout(names = geom_track(side = "right", width = 1.8, gap = 0.4, limits = c(0,
    1)))
p <- do.call(ggideogram, c(list(data = view, tracks = tracks,
  padding = 2, axis_side = 'right', axis_breaks = seq(5000, 13000, 2000)), common))
rects <- project_chr_interval(p, visible, 'Chr', 'EdgeStart', 'End')
p <- p + ggiraph::geom_rect_interactive(data = rects,
    aes(xmin = .x_start, xmax = .x_end, ymin = .y_start - 0.375,
      ymax = .y_start + 0.375, fill = Strand, data_id = ID, tooltip = Tip),
    colour = '#626A73', linewidth = 0.22, inherit.aes = FALSE, key_glyph = gallery_key) +
  geom_track(data = visible, track = "names", mapping = aes(chr = Chr, x = Mid, y = 0.5,
      label = Name, data_id = ID, tooltip = Tip), geom = ggiraph::geom_text_interactive(size = 2.9,
      show.legend = FALSE)) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  gene_colours() + gallery_theme()
save_interactive(p, d, '02-gene-intervals', height = 45, category = 'Strand', url = 'URL')

# 3. Native interactive columns in two bp-aligned tracks, real 1 Mb counts.
data(human_karyotype, package = 'ggideogram')
data(gene_density, package = 'ggideogram')
data(LTR_density, package = 'ggideogram')
view <- chr_view(human_karyotype, '1', 110000000, 140000000)
g <- gene_density; g$Category <- 'Genes'
l <- LTR_density; l$Category <- 'LTR'
d <- rbind(g, l)
d$Rate <- d$Value * 1e6 / (d$End - d$Start + 1)
d$Mid <- (d$Start - 1 + d$End) / 2
d <- view_chr_data(d, view, position = 'Mid')
d$ID <- paste(d$Category, d$Chr, d$Start, sep = ':')
d$Tip <- sprintf('%s\nChr%s: %s–%s bp\nCount: %s; %.1f / Mb',
  d$Category, d$Chr, d$Start, d$End, d$Value, d$Rate)
tracks <- track_layout(genes = geom_track(side = "right", width = 2.6, gap = 0.7, limits = c(0,
    150), reverse = TRUE), ltr = geom_track(side = "right", width = 2.6, gap = 1, limits = c(0,
    650), reverse = TRUE))
p <- ggideogram(view, orientation = 'horizontal', chromosome_width = 0.75,
  max_chr_length = 55, tracks = tracks, axis = TRUE, axis_side = 'left',
  axis_breaks = seq(110000000, 140000000, 5000000), axis_units = 'Mb',
  name_size = 2.9, name_colour = 'black', colour = '#626A73', linewidth = 0.22,
  axis_size = 2.9, axis_colour = 'black', axis_linewidth = 0.22, axis_tick_length = 1.2)
for (category in c('Genes', 'LTR')) {
  p <- p + geom_track(data = d[d$Category == category, ], track = if (category == "Genes") "genes" else "ltr",
      mapping = aes(chr = Chr, x = Mid, y = Rate, fill = Category, data_id = ID, tooltip = Tip),
      geom = ggiraph::geom_col_interactive(width = 7e+05, colour = NA, key_glyph = gallery_key))
}
p <- p + geom_track(track = "genes", axis = gallery_axis(breaks = c(0, 75, 150), colour = "black",
    linewidth = 0.22, tick_length = 1.2)) +
  geom_track(track = "ltr", axis = gallery_axis(breaks = c(0, 325, 650), colour = "black",
      linewidth = 0.22, tick_length = 1.2)) +
  scale_fill_manual(values = c(Genes = '#4477AA', LTR = '#CC6677'), name = 'Feature') +
  gallery_theme(plot.margin = margin(8, 34, 6, 24))
save_interactive(p, d, '03-window-tracks', height = 65, category = 'Category')
writeLines(c(
  '01: Arabidopsis TAIR10, first 24 gene midpoints in Chr1 0–100 kb; hover, select, gene links and strand filter.',
  '02: Chr1 4.5–13 kb, three source gene spans; cropped display, original complete intervals in TSV export.',
  '03: Human Chr1 110–140 Mb, original RIdeogram gene/LTR counts normalized by actual window width.',
  'Native point, rect and column geoms use ggiraph; static PNG/PDF retain the same plotted data.',
  'Source GFF3: Ensembl Plants release 62, TAIR10 / Araport11.',
  'Open an HTML file alongside its lib directory; controls select marks and export original rows.'
), file.path(output, 'captions.txt'))
