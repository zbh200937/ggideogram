# Run from the source repository:
# Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/optimized-gallery.R")'
library(ggideogram)
library(ggplot2)
output <- Sys.getenv('GGIDEOGRAM_EXAMPLE_OUTPUT', 'work/optimized-review')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
extdata <- function(file) system.file('extdata', file, package = 'ggideogram')
save_pair <- function(plot, name, width = 185, height = width) {
  ggsave(file.path(output, paste0(name, '.png')), plot, width = width, height = height,
    units = 'mm', dpi = 220, bg = 'white', device = ragg::agg_png)
  ggsave(file.path(output, paste0(name, '.pdf')), plot, width = width, height = height,
    units = 'mm', bg = 'white', device = grDevices::cairo_pdf)
}
style <- function() theme(legend.position = 'bottom', legend.key = element_blank(),
  legend.key.width = grid::unit(4, 'mm'), legend.key.height = grid::unit(3, 'mm'),
  plot.title = element_text(size = 11, face = 'plain', margin = margin(b = 5, unit = 'mm')),
  plot.caption = element_text(size = 8, hjust = 0))
strand <- c('+' = '#4477AA', '-' = '#CC6677')
strand_scales <- function() list(scale_fill_manual(values = strand, name = 'Strand'),
  scale_colour_manual(values = strand, name = 'Strand'), guides(colour = 'none'))

# 1. Full source models and annotation summaries share one local circular axis.
# The exon coverage is a union fraction, not read depth or expression.
f <- read_chr_features(extdata('arabidopsis-first-genes.gff3'))
genes <- f[f$Type == 'gene', ]
genes$Mid <- (genes$Start - 1 + genes$End) / 2
tx <- f[f$Type == 'mRNA', ]
f$Transcript <- ifelse(f$Type == 'mRNA', f$ID, f$Parent)
models <- f[f$Transcript %in% tx$ID, ]
models$Gene <- sub('gene:', '', tx$Parent[match(models$Transcript, tx$ID)])
k <- data.frame(Chr = '1', Start = 0, End = 30427671, Label = 'Chr1')
semantic <- as_ideogram_data(k, aes(chr = Chr, start = Start, end = End, label = Label))
count <- bin_genome(genes, k, window = 1000, method = 'count')
coverage <- bin_genome(f[f$Type == 'exon', ], k, window = 1000, method = 'coverage')
count$Mid <- (count$Start - 1 + count$End) / 2
coverage$Mid <- (coverage$Start - 1 + coverage$End) / 2
local_tracks <- list(
  models = geom_track(side = 'inner', width = 4, gap = .8, data = models,
    aes(chr = Chr, start = Start, end = End, gene = Gene, type = Type,
      strand = Strand, fill = Strand), label = 'Genes',
    layers = list(geom_genemodel(mode = 'gene', labels = FALSE, block_height = .55))),
  count = geom_track(side = 'inner', width = 2.8, gap = .8, data = count,
    aes(chr = Chr, x = Mid, y = Value), limits = c(0, 1),
    geom = geom_col(width = 800, fill = '#9BBACD', colour = NA), label = 'Count',
    axis = list(breaks = c(0, 1))),
  coverage = geom_track(side = 'inner', width = 2.3, gap = 2, data = coverage,
    aes(chr = Chr, x = Mid, y = Value), limits = c(0, 1),
    layers = list(geom_line(linewidth = .32, colour = '#526C76')), label = 'Cover',
    axis = list(breaks = c(0, 1)))
)
local_base <- function(opening_angle = 65) ggideogram(chr_view(semantic, '1', 4500, 14500),
  orientation = 'circular', radius = 22, opening_angle = opening_angle,
  tracks = local_tracks, axis = TRUE, axis_side = 'inner', axis_units = 'kb',
  axis_breaks = c(5000, 8000, 11000, 14000), axis_gap = .8,
  fill = '#EEF0F2', colour = '#626A73', linewidth = .22) +
  geom_locus(data = genes, aes(chr = Chr, position = Mid, colour = Strand),
    side = 'outer', gap = 0, size = 1.7, show.legend = FALSE) +
  strand_scales() + style() +
  labs(title = 'Arabidopsis · Chr1 4.5–14.5 kb')
local_labels <- function(width, opening_angle = 65) local_base(opening_angle) +
  geom_locus(geom = "text", position = "spread", data = genes, aes(chr = Chr, position = Mid, label = Name),
    side = 'outer', label_width = width, gap = .35)
# Declare label space for each output slot; physical text stays at 3 mm.
p_local <- local_labels(10)
save_pair(p_local, '01-local-shared-circle', 155)
save_pair(local_labels(18, opening_angle = 65), '01-local-shared-circle-120', 120)
save_pair(local_labels(8), '01-local-shared-circle-200', 200)

# 2. Two genomes with identical chromosome names, outer annotation tracks and
# source-coordinate interval relations. Native line + point share one track.
kr <- read.delim(extdata('rice-comparison-karyotype.tsv'), colClasses = 'character')
kr$Start <- as.numeric(kr$Start); kr$End <- as.numeric(kr$End)
kr$Key <- with(kr, chr_key(Genome, Chr, Assembly))
dr <- as_ideogram_data(kr, aes(chr = Chr, start = Start, end = End,
  genome = Genome, assembly = Assembly, label = Label))
gr <- read.delim(extdata('rice-comparison-genes.tsv'))
gr$Key <- with(gr, chr_key(Genome, Chr, Assembly))
counts <- bin_genome(data.frame(Chr = gr$Key, Start = gr$Start, End = gr$End),
  dr, window = 1e6, method = 'count')
counts$Mid <- (counts$Start - 1 + counts$End) / 2
blocks <- read.delim(extdata('rice-comparison-blocks.tsv'))
blocks$Key1 <- with(blocks, chr_key(Genome1, Chr1, Assembly1))
blocks$Key2 <- with(blocks, chr_key(Genome2, Chr2, Assembly2))
rice <- counts[counts$Chr == kr$Key[kr$Genome == 'rice'], ]
wild <- counts[counts$Chr == kr$Key[kr$Genome == 'wild'], ]
count_limit <- ceiling(max(counts$Value) / 50) * 50
p_rice <- ggideogram(dr, ncol = 2, chromosome_gap = 10,
  axis = FALSE, reverse_chr = kr$Key[kr$Genome == 'wild'],
  fill = '#EEF0F2', colour = '#626A73', linewidth = .22) +
  geom_chrlink(type = 'interval', data = blocks,
    aes(chr1 = Key1, start1 = Start1, end1 = End1,
      chr2 = Key2, start2 = Start2, end2 = End2,
      orientation = Orientation, fill = Orientation), alpha = .5, colour = NA) +
  geom_track(track = 'rice_count', side = 'left', width = 3, data = rice,
    aes(chr = Chr, x = Mid, y = Value), label = 'Genes / window', limits = c(0, count_limit),
    layers = list(geom_line(colour = '#4477AA', linewidth = .3),
                  geom_point(colour = '#4477AA', size = .75)),
    axis = list(chr = rice$Chr[1], breaks = c(0, count_limit), position = 'end')) +
  geom_track(track = 'wild_count', side = 'right', width = 3, data = wild,
    aes(chr = Chr, x = Mid, y = Value), label = NULL, limits = c(0, count_limit),
    geom = geom_col(width = 8e5, fill = '#CC6677', alpha = .8),
    axis = list(chr = wild$Chr[1], breaks = c(0, count_limit), position = 'start')) +
  geom_chr(chr = rice$Chr[1], component = 'axis', side = 'left',
    units = 'Mb', breaks = c(0, 2e7, 4e7)) +
  geom_chr(chr = wild$Chr[1], component = 'axis', side = 'right',
    units = 'Mb', breaks = c(0, 2e7)) +
  scale_fill_manual(values = c('+' = '#86A5BC', '-' = '#CC6677'),
    name = 'Block order', labels = c('Same', 'Reversed')) + style()
save_pair(p_rice, '02-rice-outer-tracks', 185, 155)

# 3. A real 12-gene annotation cluster, with radial text edges on one arc.
cluster <- read_chr_features(extdata('arabidopsis-gene-cluster.gff3'))
cluster <- cluster[cluster$Type == 'gene', ]
cluster <- head(cluster[order(cluster$Start), ], 12)
cluster$Mid <- (cluster$Start - 1 + cluster$End) / 2
cluster$Label <- ifelse(is.na(cluster$Name) | !nzchar(cluster$Name),
  sub('gene:', '', cluster$ID), cluster$Name)
cluster_view <- semantic
p_dense_base <- ggideogram(cluster_view, orientation = 'circular', radius = 25,
  opening_angle = 32, axis = TRUE, axis_side = 'inner', axis_units = 'Mb', axis_n = 4,
  fill = '#EEF0F2', colour = '#626A73', linewidth = .22) +
  geom_locus(data = cluster, aes(chr = Chr, position = Mid, colour = Strand),
    side = 'outer', gap = 0, size = 1.5, show.legend = FALSE) +
  scale_colour_manual(values = strand, name = 'Strand') + style()
dense_labels <- function(width) p_dense_base +
  geom_locus(geom = "text", position = "spread", data = cluster, aes(chr = Chr, position = Mid, label = Label),
    side = 'outer', label_width = width, gap = .4)
p_dense <- dense_labels(18)
save_pair(p_dense, '03-dense-label-ring', 185)
save_pair(dense_labels(34), '03-dense-label-ring-120', 120)
save_pair(dense_labels(16), '03-dense-label-ring-200', 200)

# 4. Source feature IDs and a reversed interval relation share the local axis.
pairs <- read.delim(extdata('arabidopsis-synteny-pairs.tsv'))
pairs <- pairs[pairs$Block == 15, ]
node_source <- unique(rbind(
  data.frame(ID = pairs$Gene1, Chr = pairs$Chr1, Start = pairs$Start1, End = pairs$End1),
  data.frame(ID = pairs$Gene2, Chr = pairs$Chr2, Start = pairs$Start2, End = pairs$End2)))
nodes <- chr_nodes(node_source, aes(id = ID, chr = Chr, start = Start, end = End))
blocks_a <- read.delim(extdata('arabidopsis-synteny-blocks.tsv'))
blocks_a <- blocks_a[blocks_a$Block == 15, ]
node_source$Mid <- (node_source$Start - 1 + node_source$End) / 2
p_nodes <- ggideogram(chr_view(semantic, '1', 23150000, 23540000),
  orientation = 'horizontal', axis = TRUE, axis_side = 'left', axis_units = 'Mb',
  axis_breaks = c(23150000, 23250000, 23350000, 23450000, 23540000),
  fill = '#EEF0F2', colour = '#626A73', linewidth = .22) +
  geom_chrlink(type = 'interval', data = blocks_a,
    aes(chr1 = Chr1, start1 = Start1, end1 = End1,
      chr2 = Chr2, start2 = Start2, end2 = End2, orientation = Orientation),
    fill = '#CC6677', alpha = .22, colour = NA) +
  geom_chrlink(data = pairs, nodes = nodes, aes(from = Gene1, to = Gene2),
    colour = '#4477AA', linewidth = .23, alpha = .9) +
  geom_locus(data = node_source, aes(chr = Chr, position = Mid),
    side = 'right', gap = 0, size = 1.4, colour = '#4477AA') + style()
save_pair(p_nodes, '04-local-node-arcs', 185, 55)

# 5. Ordinary independent panels use patchwork; the ideogram remains a native
# ggplot component; the summary preserves its own categorical coordinates.
if (requireNamespace('patchwork', quietly = TRUE)) {
  summary <- data.frame(Group = c('First genes', 'Cluster genes'),
    N = c(sum(count$N), nrow(cluster)))
  summary_plot <- ggplot(summary, aes(Group, N)) + geom_col(fill = '#9BBACD', width = .65) +
    scale_y_continuous(breaks = function(limits) unique(round(pretty(limits, n = 4)))) +
    labs(x = NULL, y = 'Bundled genes') + theme_classic(base_size = 11)
  composed <- patchwork::wrap_plots(p_local, summary_plot, widths = c(2, 1)) +
    patchwork::plot_layout(guides = 'keep')
  save_pair(composed, '05-external-composition', 230, 155)
}

writeLines(c(
  '01: Arabidopsis TAIR10 / Araport11 Chr1 4.5–14.5 kb. From the chromosome inward: gene-level exon/CDS/UTR structures, midpoint gene counts per 1 kb bin, exon-union coverage fraction and original-bp axis. Outer radial labels start on one arc; leaders and points keep their source gene anchors.',
  '02: O. sativa IRGSP-1.0 and O. rufipogon OR_W1943 Chr1, complete bundled gene annotations summarized per 1 Mb. Wild-rice display is reversed. The eight LASTZ_NET alignment blocks are the two longest from each of four queried regions; the source block boundaries are retained.',
  '03: First 12 genes ordered by source position in the bundled Arabidopsis gene-cluster annotation, displayed on the full Chr1 coordinate range. Near text edges and leader ends share one circular baseline; radial text faces upright. Native points retain the source anchors. Label track widths are declared for each output slot; physical label size stays at 3 mm.',
  '04: Arabidopsis MCScanX block 15, six real reversed gene pairs in Chr1. ID-based node connections and the interval band use the same original-bp projection and default arcs.',
  '05: Independent annotation-count summary combined by patchwork; source counts are of the bundled examples, not whole-genome gene totals.',
  'Data sources: extdata/arabidopsis-first-genes.gff3, arabidopsis-gene-cluster.gff3, arabidopsis-synteny-source.txt and multi-genome-source.txt.'
), file.path(output, 'captions.txt'))
