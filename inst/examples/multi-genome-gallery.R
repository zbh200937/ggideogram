# Run: micromamba run -n multiomics Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/multi-genome-gallery.R")'
# Regenerate inputs: python3 inst/examples/prepare-multi-genome-data.py
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))
input <- function(name) {
  d <- read.delim(system.file('extdata', name, package = 'ggideogram'), colClasses = 'character')
  for (field in intersect(c('Start', 'End'), names(d))) d[[field]] <- as.numeric(d[[field]])
  d
}
semantic <- function(k) as_ideogram_data(k,
  aes(chr = Chr, start = Start, end = End, genome = Genome, assembly = Assembly,
      homolog = Homolog, label = Label))
output <- 'work/api-optimization/examples/multi-genome-gallery'
style <- list(name_colour = 'black', name_size = 2.9, colour = '#626A73',
  linewidth = 0.22, axis_colour = 'black', axis_size = 2.9, axis_linewidth = 0.22,
  axis_tick_length = 1.2, axis_units = 'Mb', fill = '#F0F2F4', padding = 1)

# 1. All 19 nuclear chromosomes from each Populus haplotype, common bp scale.
k <- input('poplar-haplotypes.tsv')
k$Key <- with(k, chr_key(Genome, Chr, Assembly))
p <- do.call(ggideogram, c(list(data = semantic(k), order_by = 'genome',
  ncol = 19, max_chr_length = 22, chromosome_width = 0.65, chromosome_gap = 2.2,
  row_gap = 4, tracks = track_layout(body = geom_track(side = "overlay", width = 0.65)),
  axis = k$Key[k$Chr == '01'], axis_breaks = seq(0, 4e7, 1e7)), style)) +
  geom_chr(component = "fill", data = transform(k, Start = 1), track = "body", aes(chr = Key,
      start = Start, end = End, fill = Haplotype), key_glyph = gallery_key) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  scale_fill_manual(values = c(A = '#AACCEE', B = '#F4CBD2'), name = 'Haplotype') +
  gallery_theme(plot.margin = margin(14, 16, 6, 36))
p <- p + labs(title = 'Populus trichocarpa · phased haplotypes',
  caption = 'GWHERCL00000000 / GWHERCM00000000. Nineteen chromosomes per haplotype.')
gallery_save(p, output, '01-poplar-haplotypes', height = 160)

# 2. A/B/D subgenomes, grouped horizontally; identical chromosome names remain distinct.
k <- input('wheat-subgenomes.tsv')
k$Key <- with(k, chr_key(Genome, Chr, Assembly))
p <- do.call(ggideogram, c(list(data = semantic(k), order_by = 'genome',
  orientation = 'horizontal', genome_order = c('A', 'B', 'D'), ncol = 7,
  max_chr_length = 18, chromosome_width = 0.65, chromosome_gap = 2.2, row_gap = 3.5,
  axis_side = 'right',
  tracks = track_layout(body = geom_track(side = "overlay", width = 0.65)),
  axis = k$Key[k$Chr == '7'], axis_breaks = function(range) seq(0, range[2], 2e8)), style)) +
  geom_chr(component = "fill", data = transform(k, Start = 1), track = "body", aes(chr = Key,
      start = Start, end = End, fill = Subgenome), key_glyph = gallery_key) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  scale_fill_manual(values = c(A = '#AACCEE', B = '#F4CBD2', D = '#DCEAC4'), name = 'Subgenome') +
  gallery_theme(plot.margin = margin(12, 16, 12, 27))
p <- p + labs(title = 'Wheat · A, B and D subgenomes',
  caption = 'Chinese Spring CS RefSeq v2.1. Twenty-one nuclear chromosomes; common bp scale.')
gallery_save(p, output, '02-wheat-subgenomes', height = 87)

# 3. Full Chr1 in rice and wild rice; sampled LASTZ blocks and complete gene density.
k <- input('rice-comparison-karyotype.tsv')
k$Key <- with(k, chr_key(Genome, Chr, Assembly))
genes <- input('rice-comparison-genes.tsv')
genes$Chr <- with(genes, chr_key(Genome, Chr, Assembly))
counts <- bin_genome(genes, semantic(k), window = 1e6, method = 'count')
counts$Mid <- (counts$Start - 1 + counts$End) / 2
blocks <- read.delim(system.file('extdata', 'rice-comparison-blocks.tsv', package = 'ggideogram'))
blocks$Key1 <- with(blocks, chr_key(Genome1, Chr1, Assembly1))
blocks$Key2 <- with(blocks, chr_key(Genome2, Chr2, Assembly2))
high <- ceiling(max(counts$Rate) / 50) * 50
tracks <- track_layout(body = geom_track(side = "overlay", width = 0.9), upper_density = geom_track(side = "left",
    width = 3, gap = 0.8, limits = c(0, high)), lower_density = geom_track(side = "right",
    width = 3, gap = 0.8, limits = c(0, high), reverse = TRUE))
upper_chr <- k$Key[1]
lower_chr <- k$Key[2]
p <- do.call(ggideogram, c(list(data = semantic(k), order_by = 'homolog',
  orientation = 'horizontal', max_chr_length = 55, chromosome_width = 0.9,
  chromosome_gap = 0.8, tracks = tracks), style)) +
  geom_chrlink(type = "interval", data = blocks, aes(chr1 = Key1, start1 = Start1, end1 = End1,
      chr2 = Key2, start2 = Start2, end2 = End2, orientation = Orientation), fill = "#AAB2BD",
      alpha = 0.3, colour = NA) +
  geom_chrlink(type = "point", data = transform(blocks, Pos1 = (Start1 + End1)/2, Pos2 = (Start2 +
      End2)/2), aes(chr1 = Key1, position1 = Pos1, chr2 = Key2, position2 = Pos2), colour = "#818E9C",
      linewidth = 0.22, alpha = 0.7) +
  geom_chr(component = "fill", data = counts, track = "body", aes(chr = Chr, start = Start,
      end = End, fill = Rate)) +
  geom_track(data = counts[counts$Chr == upper_chr, ], track = "upper_density", mapping = aes(chr = Chr,
      x = Mid, y = Rate), geom = ggplot2::geom_line(colour = "#4477AA", linewidth = 0.25)) +
  geom_track(data = counts[counts$Chr == lower_chr, ], track = "lower_density", mapping = aes(chr = Chr,
      x = Mid, y = Rate), geom = ggplot2::geom_line(colour = "#4477AA", linewidth = 0.25)) +
  geom_track(track = "upper_density", axis = gallery_axis(chr = upper_chr, breaks = c(0, high/2,
      high), colour = "black", linewidth = 0.22, tick_length = 1.2)) +
  geom_track(track = "lower_density", axis = gallery_axis(chr = lower_chr, breaks = c(0, high/2,
      high), colour = "black", linewidth = 0.22, tick_length = 1.2)) +
  geom_chr(component = "axis", chr = upper_chr, side = "left", units = "Mb", breaks = function(range) seq(0,
      range[2], 1e+07), size = 2.9, colour = "black", linewidth = 0.22, tick_length = 1.2) +
  geom_chr(component = "axis", chr = lower_chr, side = "right", units = "Mb", breaks = function(range) seq(0,
      range[2], 1e+07), size = 2.9, colour = "black", linewidth = 0.22, tick_length = 1.2) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  scale_fill_gradient(low = '#F2F5F8', high = '#4477AA', limits = c(0, high),
    name = 'Genes / Mb', guide = guide_colourbar(barwidth = grid::unit(32, 'mm'),
      barheight = grid::unit(2, 'mm'))) +
  gallery_theme(plot.margin = margin(12, 30, 9, 86))
p <- p + labs(title = 'Rice and wild rice · Chr1',
  caption = 'IRGSP-1.0 / OR_W1943. Genes per Mb; eight sampled LASTZ_NET alignment blocks.')
gallery_save(p, output, '03-rice-comparison', height = 80)
writeLines(c(
  '01: Populus trichocarpa, 19 nuclear chromosomes in each phased haplotype; shared bp scale.',
  '02: Chinese Spring CS RefSeq v2.1, all 21 nuclear chromosomes grouped into A/B/D subgenomes.',
  '03: Rice IRGSP-1.0 and wild rice OR_W1943, Chr1; complete protein-coding gene counts per 1 Mb.',
  '03: Bands are the two longest Chr1–Chr1 LASTZ_NET blocks from each of four rice queries:',
  '1000000–1999999, 7000000–7999999, 14000000–14199999 and 22200000–22399999 bp.',
  'Their midpoint connections improve visibility; original endpoint widths and directions are preserved.',
  'Gene counts use each feature midpoint once and normalize the actual final window width.',
  'Sources and assembly accessions: inst/extdata/multi-genome-source.txt.'
), file.path(output, 'captions.txt'))
