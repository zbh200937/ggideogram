# Shared chromosome coordinates, native tracks and source-node links.
# After installing ggideogram:
# Rscript -e 'source(system.file("examples", "combined-gallery.R", package = "ggideogram"))'

library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))

output <- Sys.getenv('GGIDEOGRAM_EXAMPLE_OUTPUT', 'work/combined-review')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
extdata <- function(file) system.file('extdata', file, package = 'ggideogram')
save_pair <- function(plot, name, width, height) {
  ggsave(file.path(output, paste0(name, '.png')), plot, device = ragg::agg_png,
    width = width, height = height, units = 'mm', dpi = 150, bg = 'white')
  ggsave(file.path(output, paste0(name, '.pdf')), plot, device = grDevices::cairo_pdf,
    width = width, height = height, units = 'mm', bg = 'white')
}
style <- function() theme(text = element_text(size = 9),
  legend.position = 'bottom', legend.key = element_blank(),
  legend.text = element_text(size = 8), legend.title = element_text(size = 8),
  legend.key.size = grid::unit(4, 'mm'), legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = gallery_legend_spacing,
  plot.title = element_text(size = 11, margin = margin(b = 5, unit = 'mm')))
palette <- c(rice = '#4477AA', wild = '#CC6677')
genome_colour <- function() scale_colour_manual(values = palette, limits = names(palette),
  drop = FALSE, name = 'Genome')
genome_fill <- function() scale_fill_manual(values = palette, limits = names(palette),
  drop = FALSE, name = 'Genome')
guide_style <- function() guides(colour = guide_legend(override.aes =
  list(size = 1.6, shape = 16, stroke = .22, linewidth = .24, fill = NA)))

# 1. Rice and wild rice: distinguish identical chromosome names by genome and
# assembly, then summarize their complete gene annotations in 1 Mb windows.
k <- read.delim(extdata('rice-comparison-karyotype.tsv'),
  colClasses = c(Chr = 'character'))
k$Key <- with(k, chr_key(Genome, Chr, Assembly))
rice_model <- as_ideogram_data(k, aes(chr = Chr, start = Start, end = End,
  genome = Genome, assembly = Assembly, homolog = Homolog, label = Label))
rice_genes <- read.delim(extdata('rice-comparison-genes.tsv'),
  colClasses = c(Chr = 'character'))
rice_genes$Chr <- with(rice_genes, chr_key(Genome, Chr, Assembly))
rice_genes$Pos <- (rice_genes$Start - 1 + rice_genes$End) / 2
rice_genes$Node <- paste(rice_genes$Genome, rice_genes$Gene, sep = ':')
rice_bins <- bin_genome(rice_genes, rice_model, window = 1e6, method = 'count')
rice_bins$Pos <- (rice_bins$Start - 1 + rice_bins$End) / 2
rice_bins$Genome <- k$Genome[match(rice_bins$Chr, k$Key)]
rice_coverage <- bin_genome(rice_genes, rice_model, window = 1e6, method = 'coverage')
rice_coverage$Pos <- (rice_coverage$Start - 1 + rice_coverage$End) / 2
rice_coverage$Genome <- k$Genome[match(rice_coverage$Chr, k$Key)]

# Interval links retain the selected alignment boundaries. Markers identify
# the nearest annotated gene to each block midpoint in each assembly.
blocks <- read.delim(extdata('rice-comparison-blocks.tsv'),
  colClasses = c(Chr1 = 'character', Chr2 = 'character'))
blocks$Key1 <- with(blocks, chr_key(Genome1, Chr1, Assembly1))
blocks$Key2 <- with(blocks, chr_key(Genome2, Chr2, Assembly2))
blocks$Mid1 <- (blocks$Start1 - 1 + blocks$End1) / 2
blocks$Mid2 <- (blocks$Start2 - 1 + blocks$End2) / 2
anchors <- do.call(rbind, lapply(seq_len(nrow(blocks)), function(i) {
  a <- rice_genes[rice_genes$Chr == blocks$Key1[i], ]
  b <- rice_genes[rice_genes$Chr == blocks$Key2[i], ]
  rbind(a[which.min(abs(a$Pos - blocks$Mid1[i])), ],
    b[which.min(abs(b$Pos - blocks$Mid2[i])), ])
}))
anchors <- anchors[!duplicated(anchors$Node), ]
rice_labels <- anchors[!duplicated(anchors$Chr), ]
rice_labels$Rate <- rice_bins$Rate[vapply(seq_len(nrow(rice_labels)), function(i)
  which(rice_bins$Chr == rice_labels$Chr[i] & rice_labels$Pos[i] >= rice_bins$Start - 1 &
    rice_labels$Pos[i] <= rice_bins$End), integer(1))]
label_chr <- match(rice_labels$Chr, k$Key)
label_fraction <- (rice_labels$Pos - k$Start[label_chr]) /
  (k$End[label_chr] - k$Start[label_chr])
label_fraction[rice_labels$Genome == 'wild'] <- 1 - label_fraction[rice_labels$Genome == 'wild']
rice_labels$Hjust <- as.numeric(label_fraction > .5)

# Each genome owns its outward tracks. Rate is genes per Mb; Count is the
# actual number per window, including the shorter final window. Coverage is
# the fraction of each window covered by the union of annotated genes.
rice_chr <- k$Key[k$Genome == 'rice']
wild_chr <- k$Key[k$Genome == 'wild']
rate_limits <- c(0, ceiling(max(rice_bins$Rate) / 50) * 50)
count_limits <- c(0, ceiling(max(rice_bins$Value) / 50) * 50)
rate_track <- function(chr, side, axis_position) geom_track(chr = chr, side = side,
  width = 4, gap = 1.6, data = rice_bins,
  mapping = aes(chr = Chr, x = Pos, y = Rate, colour = Genome),
  limits = rate_limits, label = 'Rate (genes/Mb)',
  axis = gallery_axis(position = axis_position, breaks = rate_limits),
  layers = list(geom_col(width = 5e5, fill = 'grey85', colour = NA,
    show.legend = FALSE), geom_line(linewidth = .24, show.legend = FALSE),
    geom_point(size = 1.6, stroke = .22)))
coverage_track <- function(chr, side, axis_position) geom_track(chr = chr, side = side,
  width = 4, gap = 1.6, data = rice_coverage,
  mapping = aes(chr = Chr, x = Pos, y = Value, colour = Genome),
  limits = c(0, 1), label = 'Cover (fraction)',
  axis = gallery_axis(position = axis_position, breaks = c(0, 1)),
  layers = list(geom_line(linewidth = .24, show.legend = FALSE)))
label_track <- function(chr, side) geom_track(chr = chr, side = side,
  width = 3.5, gap = 1.6, data = rice_labels, limits = rate_limits,
  mapping = aes(chr = Chr, x = Pos, y = Rate, label = Gene, colour = Genome,
    hjust = Hjust),
  layers = list(geom_text(size = 2.5, show.legend = FALSE)))
count_track <- function(chr, side, axis_position, track = NULL) geom_track(
  track = track, chr = chr, side = side,
  width = 4, gap = 1.6, data = rice_bins,
  mapping = aes(chr = Chr, x = Pos, y = Value, fill = Genome),
  limits = count_limits, label = 'Count (genes/window)',
  axis = gallery_axis(position = axis_position, breaks = count_limits),
  layers = list(geom_col(width = 5e5, alpha = .7, show.legend = FALSE)))

# Ordinary columns, line and point reuse the same raw bp/value scope; domain
# fill and locus text use their own scopes on those same chromosome axes.
rice_tracks <- list(
  wild_rate = rate_track(wild_chr, 'left', 'start'),
  rice_rate = rate_track(rice_chr, 'right', 'end'),
  wild_coverage = coverage_track(wild_chr, 'left', 'start'),
  rice_coverage = coverage_track(rice_chr, 'right', 'end'),
  blocks = geom_track(side = 'overlay', width = .8, data = anchors,
    mapping = aes(chr = Chr, start = Start, end = End, fill = Genome),
    layers = list(geom_chr(component = 'fill', show.legend = FALSE))),
  wild_labels = label_track(wild_chr, 'left'),
  rice_labels = label_track(rice_chr, 'right'))
rice_plot <- ggideogram(rice_model, orientation = 'horizontal',
  order_by = 'homolog', genome_order = c('wild', 'rice'),
  reverse_chr = wild_chr, ncol = 2, chromosome_gap = 4,
  tracks = rice_tracks, max_chr_length = 90, fill = '#9ABCB7', linewidth = .24,
  axis = FALSE,
  name_size = 3.2, name_gap = 1.8) +
  geom_chr(chr = wild_chr, component = 'axis',
    side = 'left', units = 'Mb', breaks = seq(0, 4e7, 1e7)) +
  geom_chr(chr = rice_chr, component = 'axis',
    side = 'right', units = 'Mb', breaks = seq(0, 4e7, 1e7)) +
  geom_chrlink(data = blocks, type = 'interval',
    aes(chr1 = Key1, start1 = Start1, end1 = End1, chr2 = Key2,
      start2 = Start2, end2 = End2, orientation = Orientation),
    side1 = 'left', side2 = 'right', gap = 0, fill = '#AAB2BD', alpha = .4) +
  geom_chrlink(data = blocks,
    aes(chr1 = Key1, position1 = Mid1, chr2 = Key2, position2 = Mid2),
    side1 = 'left', side2 = 'right', gap = 0, colour = '#687889',
    linewidth = .24, alpha = .8) +
  geom_locus(data = anchors[anchors$Genome == 'wild', ],
    aes(chr = Chr, position = Pos, colour = Genome),
    side = 'right', gap = 0, size = 1.6, stroke = .22, show.legend = FALSE) +
  geom_locus(data = anchors[anchors$Genome == 'rice', ],
    aes(chr = Chr, position = Pos, colour = Genome),
    side = 'left', gap = 0, size = 1.6, stroke = .22, show.legend = FALSE) +
  genome_colour() + genome_fill() + guide_style() + style() +
  labs(title = 'Rice and wild rice · Chr1 annotations',
    caption = paste('Eight selected LASTZ_NET blocks: two longest from each of four query regions.',
      'Lines connect block midpoints; markers show nearest annotated genes.',
      'Nominal 1 Mb bins: counts/window, width-normalized genes/Mb and gene-union coverage.', sep = '\n'))

# A later track joins the same layout as all existing links and markers.
rice_plot <- rice_plot +
  count_track(wild_chr, 'left', 'start', track = 'wild_counts') +
  count_track(rice_chr, 'right', 'end', track = 'rice_counts')
save_pair(rice_plot, 'rice-scope-comprehensive', 185, 185)

# 2. Arabidopsis: keep the complete source GFF and chromosome extent, and let
# chr_view select the visible 7–9 kb window for the circular host and its tracks.
arab_path <- extdata('arabidopsis-first-genes.gff3')
features <- read_chr_features(arab_path)
regions <- readLines(arab_path)
regions <- regions[grepl('^##sequence-region[[:space:]]', regions)]
regions <- read.table(text = sub('^##sequence-region[[:space:]]+', '', regions),
  col.names = c('Chr', 'Start', 'End'), colClasses = c('character', 'numeric', 'numeric'))
regions$Start <- regions$Start - 1
arab_karyotype <- regions[regions$Chr == '1', ]
arab_view <- chr_view(arab_karyotype, '1', 7000, 9000)
source(system.file('examples', 'gene-model-data.R', package = 'ggideogram'))
models <- gene_model_data(features)
arab_counts <- bin_genome(features[features$Type == 'gene', ],
  arab_karyotype, window = 250, method = 'count')
arab_counts$Pos <- (arab_counts$Start - 1 + arab_counts$End) / 2
arab_coverage <- bin_genome(features[features$Type == 'exon', ],
  arab_karyotype, window = 250, method = 'coverage')
arab_coverage$Pos <- (arab_coverage$Start - 1 + arab_coverage$End) / 2

# Exon rows have no GFF ID here; derive a source key from transcript and interval.
# Link the first and last visible exons of one transcript through those node IDs.
arab_exons <- features[features$Type == 'exon', ]
arab_exons$Node <- paste(arab_exons$Parent, arab_exons$Start, arab_exons$End, sep = ':')
arab_exons <- arab_exons[!duplicated(arab_exons$Node), ]
arab_exons$Pos <- (arab_exons$Start - 1 + arab_exons$End) / 2
exon_anchors <- arab_exons[arab_exons$Pos >= 7000 & arab_exons$Pos <= 9000, ]
exon_anchors <- exon_anchors[exon_anchors$Parent == exon_anchors$Parent[1], ]
exon_anchors <- exon_anchors[c(1, nrow(exon_anchors)), ]
arab_nodes <- chr_nodes(arab_exons, aes(id = Node, chr = Chr, start = Start, end = End))
arab_edges <- data.frame(From = exon_anchors$Node[1], To = exon_anchors$Node[2])
arab_tracks <- list(
  genes = geom_track(side = 'outer', width = 3.5, gap = .6, data = models,
    mapping = aes(chr = Chr, start = Start, end = End, gene = Gene,
      type = Type, strand = Strand, fill = Strand), label = 'Gene',
    layers = list(geom_genemodel(mode = 'gene', labels = FALSE, block_height = .5))),
  exons = geom_track(side = 'overlay', width = .8, data = arab_exons,
    mapping = aes(chr = Chr, start = Start, end = End, fill = Strand),
    layers = list(geom_chr(component = 'fill'))),
  count = geom_track(side = 'inner', width = 2.2, gap = 1.7, data = arab_counts,
    mapping = aes(chr = Chr, x = Pos, y = Value), label = 'Genes',
    axis = gallery_axis(breaks = c(0, 1), labels = as.character),
    layers = list(geom_col(width = 180, fill = '#B6C7DF'),
      geom_point(size = 1.6, stroke = .22, colour = '#4477AA'))),
  coverage = geom_track(side = 'inner', width = 2.2, gap = 1.7, limits = c(0, 1), data = arab_coverage,
    mapping = aes(chr = Chr, x = Pos, y = Value), label = 'Cover',
    axis = gallery_axis(breaks = c(0, 1), labels = as.character),
    layers = list(geom_line(linewidth = .24, colour = '#438776'))))
arab_plot <- ggideogram(arab_view, orientation = 'circular',
  radius = 28, opening_angle = 60, reverse_chr = '1', tracks = arab_tracks, axis = TRUE,
  axis_breaks = c(7000, 7500, 8000, 8500), axis_units = 'kb',
  name_size = 3.2, fill = '#9ABCB7', linewidth = .24) +
  geom_chrlink(data = arab_edges, nodes = arab_nodes, aes(from = From, to = To),
    linewidth = .24, colour = '#A26995') +
  geom_locus(data = exon_anchors, aes(chr = Chr, position = Pos), track = 'coverage', track_position = 1,
    size = 1.6, stroke = .22, colour = '#A26995') +
  scale_fill_manual(values = c('+' = '#4477AA', '-' = '#CC6677'), name = 'Strand') +
  style() + labs(title = 'Arabidopsis · ARV1 locus, Chr1 7–9 kb',
    caption = paste('250 bp windows: blue, gene counts; green, exon-union coverage (0–1).',
      'Purple link: first and last visible exons of transcript AT1G01020.2.', sep = '\n'))
save_pair(arab_plot, 'arabidopsis-scope-comprehensive', 220, 220)
