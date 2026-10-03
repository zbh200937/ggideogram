# Run: micromamba run -n multiomics Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/connection-gallery.R")'
# MCScanX's official Arabidopsis example, default analysis; see extdata provenance.
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))
blocks <- read.delim(system.file('extdata', 'arabidopsis-synteny-blocks.tsv', package = 'ggideogram'))
pairs <- read.delim(system.file('extdata', 'arabidopsis-synteny-pairs.tsv', package = 'ggideogram'))
pairs$Position1 <- (pairs$Start1 + pairs$End1) / 2
pairs$Position2 <- (pairs$Start2 + pairs$End2) / 2
k <- data.frame(Chr = c('1', '5'), Start = 0, End = c(30427671, 26975502))
output <- 'work/api-optimization/examples/connection-gallery'
colors <- c('+' = '#4477AA', '-' = '#CC6677')
direction_scale <- function(aesthetic = 'colour') {
  fun <- if (aesthetic == 'fill') scale_fill_manual else scale_colour_manual
  fun(values = colors, name = 'Block order', breaks = c('+', '-'),
    labels = c('Same', 'Reversed'))
}
pair_base <- function() ggideogram(k, orientation = 'horizontal',
  max_chr_length = 55, chromosome_width = 0.6, chromosome_gap = 11,
  tracks = track_layout(body = geom_track(side = "overlay", width = 0.6)),
  name_colour = 'black', name_size = 3,
  fill = '#F0F2F4', colour = '#626A73', linewidth = 0.22, padding = 1) +
  geom_chr(component = "axis", chr = "1", side = "left", units = "Mb", breaks = function(range) seq(0,
      range[2], 1e+07), size = 2.7, colour = "black", linewidth = 0.22, tick_length = 1.2) +
  geom_chr(component = "axis", chr = "5", side = "right", units = "Mb", breaks = function(range) seq(0,
      range[2], 1e+07), size = 2.7, colour = "black", linewidth = 0.22, tick_length = 1.2) +
  gallery_theme(legend.box.spacing = grid::unit(4.5, 'mm'))

# 1. One actual gene pair at the middle rank of each Chr1–Chr5 block.
cross <- pairs[pairs$Chr1 != pairs$Chr2, ]
representatives <- do.call(rbind, lapply(split(cross, cross$Block), function(d) d[ceiling(nrow(d) / 2), ]))
a <- data.frame(Chr = representatives$Chr1, Pos = representatives$Position1,
  Direction = representatives$Orientation)
b <- data.frame(Chr = representatives$Chr2, Pos = representatives$Position2,
  Direction = representatives$Orientation)
p <- pair_base() + geom_chrlink(type = "point", data = representatives, aes(chr1 = Chr1, position1 = Position1,
    chr2 = Chr2, position2 = Position2, colour = Orientation), linewidth = 0.3, alpha = 0.8) +
  geom_locus(geom = "point", data = a, aes(chr = Chr, position = Pos, colour = Direction),
      side = "right", gap = 0, size = 1.3) +
  geom_locus(geom = "point", data = b, aes(chr = Chr, position = Pos, colour = Direction),
      side = "left", gap = 0, size = 1.3) + direction_scale()
gallery_save(p, output, '01-paired-genes', height = 58)

# 2. All nine actual blocks; coloured intervals show their original endpoints.
cross_blocks <- blocks[blocks$Chr1 != blocks$Chr2, ]
spans <- rbind(data.frame(Chr = cross_blocks$Chr1, Start = cross_blocks$Start1,
  End = cross_blocks$End1, Direction = cross_blocks$Orientation),
  data.frame(Chr = cross_blocks$Chr2, Start = cross_blocks$Start2,
    End = cross_blocks$End2, Direction = cross_blocks$Orientation))
p <- pair_base() + geom_chrlink(type = "interval", data = cross_blocks, aes(chr1 = Chr1, start1 = Start1,
    end1 = End1, chr2 = Chr2, start2 = Start2, end2 = End2, orientation = Orientation,
    fill = Orientation), alpha = 0.3, colour = NA, key_glyph = gallery_key) +
  geom_chr(component = "fill", data = spans, track = "body", aes(chr = Chr, start = Start,
      end = End, fill = Direction), key_glyph = gallery_key) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  direction_scale('fill')
gallery_save(p, output, '02-collinear-blocks', height = 58)

# 3. Actual reverse block 15 within Chr1, with all six constituent gene pairs.
local <- pairs[pairs$Block == 15, ]
block <- blocks[blocks$Block == 15, ]
genes <- rbind(data.frame(Chr = local$Chr1, Start = local$Start1, End = local$End1,
  Pos = local$Position1, Gene = local$Gene1, Side = 'left'),
  data.frame(Chr = local$Chr2, Start = local$Start2, End = local$End2,
    Pos = local$Position2, Gene = local$Gene2, Side = 'right'))
view <- chr_view(k, '1', 23150000, 23550000)
tracks <- track_layout(body = geom_track(side = "overlay", width = 0.6), upper = geom_track(side = "left",
    width = 5, gap = 0.8), lower = geom_track(side = "right", width = 5, gap = 0.8))
p <- ggideogram(view, orientation = 'horizontal', max_chr_length = 55,
  chromosome_width = 0.6, tracks = tracks, axis = TRUE,
  axis_breaks = seq(23200000, 23500000, 100000), axis_units = 'Mb',
  axis_colour = 'black', axis_size = 2.7, axis_linewidth = 0.22, axis_tick_length = 1.2,
  name_colour = 'black', name_size = 3, fill = '#F0F2F4', colour = '#626A73',
  linewidth = 0.22, padding = 1) +
  geom_chrlink(type = "interval", data = block, aes(chr1 = Chr1, start1 = Start1, end1 = End1,
      chr2 = Chr2, start2 = Start2, end2 = End2, orientation = Orientation), side1 = "left",
      side2 = "right", fill = "#CC6677", alpha = 0.2, colour = NA) +
  geom_chrlink(type = "point", data = local, aes(chr1 = Chr1, position1 = Position1, chr2 = Chr2,
      position2 = Position2), side1 = "left", side2 = "right", colour = "#CC6677", linewidth = 0.22,
      alpha = 0.8) +
  geom_chr(component = "fill", data = genes, track = "body", aes(chr = Chr, start = Start,
      end = End), fill = "#4477AA") +
  geom_locus(geom = "text", position = "repel", data = genes, track = c(left = "upper", right = "lower"),
      columns = 2, aes(chr = Chr, position = Pos, label = Gene, side = Side), size = 2.9,
      seed = 9) +
  geom_locus(geom = "point", data = genes[genes$Side == "left", ], aes(chr = Chr, position = Pos),
      side = "left", gap = 0, colour = "#4477AA", size = 1.3) +
  geom_locus(geom = "point", data = genes[genes$Side == "right", ], aes(chr = Chr, position = Pos),
      side = "right", gap = 0, colour = "#4477AA", size = 1.3) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) + gallery_theme()
gallery_save(p, output, '03-local-reverse-block', height = 58)
writeLines(c(
  '01: Representative gene pairs, one middle-ranked pair per Chr1–Chr5 MCScanX block.',
  '02: All nine Chr1–Chr5 blocks, endpoint spans and explicit relative order.',
  '03: Chr1 23.15–23.55 Mb, block 15 and all its six gene pairs; blue intervals are gene spans.',
  'Source: MCScanX official at.gff / at.blast; default MCScanX parameters.',
  'Block direction is the MCScanX plus/minus relationship, not gene transcription strand.',
  'Chromosome lengths: Ensembl Plants TAIR10.',
  'Input subsets and generation instructions: inst/extdata/arabidopsis-synteny-source.txt.',
  'Reference: https://thackl.github.io/gggenomes/reference/geom_link.html'
), file.path(output, 'captions.txt'))
