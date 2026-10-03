# Run: Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/gene-structure-gallery.R")'
# Arabidopsis thaliana TAIR10 / Araport11, Ensembl Plants release 62.
# Source (Chr1, features ending <= 20 kb):
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))

features <- read_chr_features(system.file('extdata', 'arabidopsis-first-genes.gff3',
  package = 'ggideogram'))
transcripts <- features[features$Type == 'mRNA', ]
# Resolve the explicit GFF Parent links; no identifier is inferred from position.
features$Transcript <- ifelse(features$Type == 'mRNA', features$ID, features$Parent)
models <- features[features$Transcript %in% transcripts$ID, ]
models$Gene <- sub('gene:', '', transcripts$Parent[match(models$Transcript, transcripts$ID)])
models$Transcript <- sub('transcript:', '', models$Transcript)
models$Part <- ifelse(models$Type == 'CDS', 'CDS',
  ifelse(models$Type %in% c('UTR', 'five_prime_UTR', 'three_prime_UTR'), 'UTR', NA))
models$Outline <- ifelse(models$Type == 'exon', 'Exon', NA)
k <- data.frame(Chr = '1', Start = 0, End = 30427671)
output <- 'work/api-optimization/examples/gene-structure-gallery'
dir.create(output, recursive = TRUE, showWarnings = FALSE)

structure_key <- gallery_key

base_plot <- function(from, to, tracks, breaks) {
  ggideogram(chr_view(k, '1', from, to), orientation = 'horizontal',
    max_chr_length = 40, chromosome_width = 1.5 * 0.35, tracks = tracks,
    axis = TRUE, axis_side = 'left', axis_breaks = breaks,
    axis_units = 'kb', axis_colour = 'black', axis_size = 2.7,
    axis_tick_length = 1.2, axis_linewidth = 0.22,
    name_colour = 'black', name_size = 3, padding = 0.5,
    fill = '#F0F2F4', colour = '#626A73', linewidth = 0.22) +
    gallery_theme(plot.margin = margin(8, 10, 6, 105))
}
part_scales <- function() list(
  scale_fill_manual(values = c('CDS' = '#4477AA', 'UTR' = '#AACCEE'),
    breaks = c('CDS', 'UTR'), na.value = NA, name = NULL),
  scale_colour_manual(values = c('Exon' = '#626A73'), breaks = 'Exon',
    na.value = NA, name = NULL),
  guides(fill = guide_legend(order = 1, override.aes = list(colour = NA)),
    colour = guide_legend(order = 2, override.aes = list(fill = NA, linewidth = 0.22))))
save_figure <- function(plot, name, height) {
  gallery_save(plot, output, name, height = height)
}

# 1. Three neighbouring genes: gene-level unions, explicit strand legend.
genes <- base_plot(2500, 15000,
  track_layout(genes = geom_track(side = "right", width = 4.5, gap = 0.7)), seq(3000, 15000, 3000)) +
  geom_genemodel(mode = "gene", data = models, track = "genes", aes(chr = Chr, start = Start,
      end = End, type = Type, strand = Strand, gene = Gene, fill = Strand), arrow_spacing_bp = 550,
      arrow_margin_bp = 180, arrow_min_bp = 360, label_size = 2.9, key_glyph = structure_key,
      lane_order = c("AT1G01010", "AT1G01020", "AT1G01030")) +
  scale_fill_manual(values = c('+' = '#4477AA', '-' = '#EE6677'), name = 'Strand')
save_figure(genes, '01-neighbouring-genes', 64)

# 2. The six annotated ARV1 transcript structures, in stable identifier order.
arv1 <- models[models$Gene == 'AT1G01020', ]
tx_order <- sort(unique(arv1$Transcript))
isoforms <- base_plot(6500, 9500,
  track_layout(genes = geom_track(side = "right", width = 9, gap = 0.7)), seq(6500, 9500, 500)) +
  geom_genemodel(mode = "transcript", data = arv1, track = "genes", aes(chr = Chr, start = Start,
      end = End, type = Type, strand = Strand, transcript = Transcript, fill = Part, colour = Outline),
      lane_order = tx_order, arrow_spacing_bp = 125, arrow_margin_bp = 45, arrow_min_bp = 90,
      label_size = 2.9, key_glyph = structure_key) +
  part_scales()
save_figure(isoforms, '02-arv1-transcripts', 68)

# 3. A clipped window, plus the number of transcript exon unions overlapping
# each 100 bp bin. This is an annotation summary, not expression/coverage data.
starts <- seq(7000, 8800, 100)
exons <- arv1[arv1$Type == 'exon', ]
counts <- data.frame(Chr = '1', Position = starts + 50,
  Value = vapply(starts, function(s) length(unique(exons$Transcript[
    exons$End > s & exons$Start <= s + 100])), integer(1)))
window <- base_plot(7000, 8900,
  track_layout(genes = geom_track(side = "right", width = 9, gap = 0.7), count = geom_track(side = "right",
      width = 4, gap = 1.2, limits = c(0, 6), reverse = TRUE)),
  seq(7000, 8800, 300)) +
  geom_genemodel(mode = "transcript", data = arv1, track = "genes", aes(chr = Chr, start = Start,
      end = End, type = Type, strand = Strand, transcript = Transcript, fill = Part, colour = Outline),
      lane_order = tx_order, arrow_spacing_bp = 80, arrow_margin_bp = 28, arrow_min_bp = 56,
      label_size = 2.9, key_glyph = structure_key) +
  geom_track(data = counts, track = "count", mapping = aes(chr = Chr, x = Position, y = Value),
      geom = ggplot2::geom_col(width = 90, fill = "#4477AA")) +
  geom_track(data = data.frame(Chr = "1", Pos = 7000, N = 3), track = "count", mapping = aes(chr = Chr,
      x = Pos, y = N), geom = ggplot2::geom_text(label = "Transcripts\nwith exon", hjust = 1.2,
      size = 2.7, colour = "black")) +
  geom_track(track = "count", axis = list(breaks = c(0, 3, 6), size = 2.7, colour = "black",
      linewidth = 0.25)) +
  part_scales()
save_figure(window, '03-window-and-exon-count', 88)
writeLines(c(
  '01: Chr1 2.5–15 kb; gene-level union of annotated exon/CDS/UTR blocks.',
  '02: Six ARV1 transcripts, AT1G01020.1–.6; arrows indicate negative strand.',
  '03: Chr1 7–8.9 kb, full models clipped to the displayed window.',
  'Bars: number of distinct transcripts with an exon overlapping each 100 bp bin (0–6).',
  'Source: Arabidopsis TAIR10 / Araport11, Ensembl Plants release 62.',
  'Blue fill: CDS; light blue fill: explicitly annotated UTR; outlines: full exon extent.',
  paste('R', getRversion(), '; ggplot2', packageVersion('ggplot2'))
), file.path(output, 'captions.txt'))
