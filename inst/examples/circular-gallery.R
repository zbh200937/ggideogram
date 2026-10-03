# Run from the repository: micromamba run -n multiomics Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/circular-gallery.R")'
library(ggideogram)
library(ggplot2)
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))
output <- 'work/api-optimization/examples/circular-gallery'
circle_theme <- function(plot) gallery_theme(plot.margin = plot$theme$plot.margin,
  legend.box.spacing = plot$theme$legend.box.spacing)
circle_style <- list(orientation = 'circular', name_size = 2.9, name_colour = 'black',
  fill = '#F0F2F4', colour = '#626A73', linewidth = 0.22, chromosome_width = 0.85,
  axis_size = 2.7, axis_colour = 'black', axis_linewidth = 0.22,
  axis_tick_length = 1.2, axis_units = 'Mb', padding = 0.8)

# 1. Complete GRCh38 karyotype and the package's original 1 Mb gene counts.
data(human_karyotype, package = 'ggideogram')
data(gene_density, package = 'ggideogram')
k <- human_karyotype
counts <- gene_density
counts$Width <- counts$End - counts$Start + 1
counts$Mid <- (counts$Start - 1 + counts$End) / 2
counts$Rate <- counts$Value / (counts$Width / 1e6)
count_max <- ceiling(max(counts$Value) / 50) * 50
rate_max <- ceiling(max(counts$Rate) / 50) * 50
gaps <- setNames(rep(2, nrow(k)), k$Chr); gaps[nrow(k)] <- 40
tracks <- track_layout(count = geom_track(side = "inner", width = 2.8, gap = 0.6, limits = c(0,
    count_max)), rate = geom_track(side = "outer", width = 1.8, gap = 0.6, limits = c(0,
    rate_max)), body = geom_track(side = "overlay", width = 0.85, limits = c(0, 1)))
centromeres <- transform(k, Mid = (CE_start + CE_end) / 2)
p_human <- do.call(ggideogram, c(list(data = k, radius = 25, gap_angle = gaps,
  tracks = tracks, axis = '1', axis_breaks = seq(0, max(k$End), 1e8)), circle_style)) +
  geom_track(data = counts, track = "count", mapping = aes(chr = Chr, x = Mid, y = Value,
      width = Width), geom = ggplot2::geom_col(position = "identity", fill = "#9BBACD",
      colour = NA)) +
  geom_track(data = counts, track = "rate", mapping = aes(chr = Chr, x = Mid, y = Rate),
      geom = ggplot2::geom_line(colour = "#4477AA", linewidth = 0.22)) +
  geom_track(data = centromeres, track = "body", mapping = aes(chr = Chr, x = Mid, y = 0.5),
      geom = ggplot2::geom_point(colour = "#A85469", size = 1.2)) +
  geom_track(track = "count", axis = list(chr = "1", position = "gap", breaks = c(0, count_max),
      size = 2.7, colour = "black", linewidth = 0.22, tick_length = 1)) +
  geom_track(track = "rate", axis = list(chr = "1", position = "gap", breaks = c(0, rate_max),
      size = 2.7, colour = "black", linewidth = 0.22, tick_length = 1)) +
  annotate('text', x = 0, y = 0, label = 'Human · GRCh38\n24 chromosomes\n\nOuter: genes / Mb\nInner: genes per window\nPink: centromere midpoint',
    size = 2.9, lineheight = 1.25, colour = 'black')
p_human <- p_human + circle_theme(p_human)
gallery_save(p_human, output, '01-human-multitrack', height = 185)

# 2. All ten deposited example blocks: nine Chr1–Chr5 blocks and Chr1 block 15.
blocks <- read.delim(system.file('extdata', 'arabidopsis-synteny-blocks.tsv', package = 'ggideogram'))
pairs <- read.delim(system.file('extdata', 'arabidopsis-synteny-pairs.tsv', package = 'ggideogram'))
k <- data.frame(Chr = c('1', '5'), Start = 0, End = c(30427671, 26975502))
ends <- rbind(data.frame(Chr = pairs$Chr1, Start = pairs$Start1, End = pairs$End1),
  data.frame(Chr = pairs$Chr2, Start = pairs$Start2, End = pairs$End2))
count <- bin_genome(ends, k, window = 1e6, method = 'count')
count$Mid <- (count$Start - 1 + count$End) / 2
count_max <- ceiling(max(count$Value) / 10) * 10
spans <- rbind(data.frame(Chr = blocks$Chr1, Start = blocks$Start1, End = blocks$End1,
  Direction = blocks$Orientation), data.frame(Chr = blocks$Chr2, Start = blocks$Start2,
    End = blocks$End2, Direction = blocks$Orientation))
representatives <- do.call(rbind, lapply(split(pairs, pairs$Block),
  function(d) d[ceiling(nrow(d) / 2), ]))
representatives$Pos1 <- (representatives$Start1 + representatives$End1) / 2
representatives$Pos2 <- (representatives$Start2 + representatives$End2) / 2
p_synteny <- do.call(ggideogram, c(list(data = k, radius = 23,
  gap_angle = c('1' = 8, '5' = 40),
  tracks = track_layout(count = geom_track(side = "inner", width = 2.8, gap = 0.8, limits = c(0,
      count_max)), body = geom_track(side = "overlay", width = 0.85)),
  axis = TRUE, axis_breaks = seq(0, max(k$End), 1e7)), circle_style)) +
  geom_chrlink(type = "interval", data = blocks, aes(chr1 = Chr1, start1 = Start1, end1 = End1,
      chr2 = Chr2, start2 = Start2, end2 = End2, orientation = Orientation, fill = Orientation),
      gap = 0.3, alpha = 0.35, colour = NA, key_glyph = gallery_key) +
  geom_chrlink(type = "point", data = representatives, aes(chr1 = Chr1, position1 = Pos1,
      chr2 = Chr2, position2 = Pos2), gap = 0.3, colour = "#65727E", linewidth = 0.18, alpha = 0.65) +
  geom_track(data = count, track = "count", mapping = aes(chr = Chr, x = Mid, y = Value,
      width = Width), geom = ggplot2::geom_col(position = "identity", fill = "#B5C4CE",
      colour = NA)) +
  geom_chr(component = "fill", data = spans, track = "body", aes(chr = Chr, start = Start,
      end = End, fill = Direction), key_glyph = gallery_key) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  geom_track(track = "count", axis = list(chr = "1", position = "gap", breaks = c(0, count_max), size = 2.7,
      colour = "black", linewidth = 0.22, tick_length = 1)) +
  scale_fill_manual(values = c('+' = '#4477AA', '-' = '#CC6677'), name = 'Block order',
    breaks = c('+', '-'), labels = c('Same', 'Reversed'))
p_synteny <- p_synteny + circle_theme(p_synteny)
gallery_save(p_synteny, output, '02-arabidopsis-synteny', height = 185)

# 3. All 21 Chinese Spring chromosomes, ordered A/B/D with explicit group gaps.
k <- read.delim(system.file('extdata', 'wheat-subgenomes.tsv', package = 'ggideogram'),
  colClasses = 'character')
k$Start <- as.numeric(k$Start); k$End <- as.numeric(k$End)
k$Key <- with(k, chr_key(Genome, Chr, Assembly))
semantic <- as_ideogram_data(k,
  aes(chr = Chr, start = Start, end = End, genome = Genome, assembly = Assembly,
    homolog = Homolog, label = Label))
gaps <- setNames(ifelse(k$Chr == '7', 8, 1.5), k$Key)
gaps[nrow(k)] <- 40
p_wheat <- do.call(ggideogram, c(list(data = semantic, radius = 23, gap_angle = gaps,
  order_by = 'genome', genome_order = c('A', 'B', 'D'),
  tracks = track_layout(body = geom_track(side = "overlay", width = 0.85)),
  axis = k$Key[k$Chr == '1'], axis_breaks = c(0, 4e8)), circle_style)) +
  geom_chr(component = "fill", data = transform(k, Start = 1), track = "body", aes(chr = Key,
      start = Start, end = End, fill = Subgenome), key_glyph = gallery_key) +
  geom_chr(component = "body", fill = NA, colour = "#626A73", linewidth = 0.22, cytoband = FALSE) +
  scale_fill_manual(values = c(A = '#AACCEE', B = '#F4CBD2', D = '#DCEAC4'), name = 'Subgenome') +
  annotate('text', x = 0, y = 0, label = 'Triticum aestivum\nChinese Spring · CS v2.1\n21 chromosomes',
    size = 2.9, lineheight = 1.25, colour = 'black')
p_wheat <- p_wheat + circle_theme(p_wheat)
gallery_save(p_wheat, output, '03-wheat-subgenomes', height = 185)

writeLines(c(
  '01: Human GRCh38; all 24 chromosomes, original 1 Mb gene counts, and supplied centromeres.',
  'Outer line: genes per Mb; inner columns: genes per input window. The final window uses its actual width.',
  '02: Arabidopsis TAIR10 Chr1 and Chr5; all ten bundled MCScanX blocks, with one middle-ranked gene pair per block.',
  'Inner columns count the displayed gene-pair endpoints per 1 Mb window; these are not genome-wide gene counts.',
  'Same/reversed is the explicit MCScanX block relationship. Curve endpoints sit inside the count track.',
  '03: Triticum aestivum Chinese Spring CS RefSeq v2.1; all 21 nuclear chromosomes ordered A/B/D.',
  'All figures allocate sectors by bp length with a shared centerline scale; declared gaps are excluded from that scale.',
  'Inputs: package data human_karyotype/gene_density; extdata arabidopsis-synteny-blocks.tsv, arabidopsis-synteny-pairs.tsv, wheat-subgenomes.tsv.',
  'Sources: the package data help pages, extdata/arabidopsis-synteny-source.txt and extdata/multi-genome-source.txt.'
), file.path(output, 'captions.txt'))
