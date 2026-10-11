# ggideogram

[![R-CMD-check](https://github.com/zbh200937/ggideogram/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zbh200937/ggideogram/actions/workflows/R-CMD-check.yaml)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](LICENSE)

**Draw chromosomes, gene models, quantitative tracks and genomic links with ggplot2.**

`ggideogram()` returns an ordinary **ggplot object**. Chromosomes and tracks share source bp coordinates, while colours, points, lines, text and legends follow ggplot2. Add native geom layers, anchor complete ggplot/grob insets, or combine ideograms with patchwork and cowplot.

The host plot owns chromosome coordinates and layout. Named tracks own space and value ranges; their layers own data mappings and visual parameters. Complete insets own their coordinates, themes and guides. Horizontal, vertical and circular layouts use the same source data and projection contract. Data summaries record their scientific definitions, while style defaults remain independently adjustable.

![Local Arabidopsis linear view with gene models, gene counts and exon coverage](inst/examples/optimized-local-linear.png)

Arabidopsis TAIR10 / Araport11, Chr1 4.5–14.5 kb. Protein-coding gene structures
retain source strand and feature boundaries; the lower tracks show their
midpoint counts per 1 kb window and the union coverage of their exons. [Figure captions](inst/examples/captions.txt)
describe the data and encodings of all gallery figures.

[Quick start](#quick-start) · [Public functions](#public-functions) · [Example gallery](#example-gallery)

## Installation

Requires R ≥ 4.1 and ggplot2 ≥ 3.5. This installation block also installs patchwork for the composition examples below.

```r
install.packages(c("remotes", "patchwork"))
remotes::install_github("zbh200937/ggideogram")
```

After cloning the repository, you can also run `R CMD INSTALL .` from the source directory. Optional features require ggrepel ≥ 0.9.6 for repelled labels, ggiraph and htmlwidgets for interactive plots, and GenomicRanges for GRanges input.

## Quick start

Run the R blocks below in order. All example data are bundled with the package.

```r
library(ggideogram)
library(ggplot2)

data(human_karyotype, package = "ggideogram")
data(gene_density, package = "ggideogram")

# Closed intervals have plotting boundaries at Start - 1 and End.
density <- transform(gene_density, Mid = (Start - 1 + End) / 2,
  Rate = Value / ((End - Start + 1) / 1e6))

p <- ggideogram(human_karyotype, chr = c("1", "2"),
  orientation = "horizontal", axis = TRUE) +
  geom_track(track = "density", side = "right", width = 3,
    data = density, mapping = aes(chr = Chr, x = Mid, y = Rate),
    label = "Genes / Mb", axis = list(breaks = c(0, 50)),
    layers = list(
      geom_line(colour = "#4477AA", linewidth = 0.4),
      geom_point(colour = "#4477AA", size = 0.7)
    )) +
  labs(title = "Human gene density") +
  theme(plot.margin = margin(5, 12, 5, 30, unit = "mm"))
p
```

The bundled `gene_density$Value` contains counts in nominal 1 Mb windows, including shorter terminal windows. `Rate` above divides each count by its actual window width in Mb. The bundled RNA markers are a random sample of 500 GENCODE annotations; their counts describe sampled records.

In track mappings, **`x` always means source bp and `y` always means the raw observation value**, regardless of the display orientation. Layers in the same named track share space and value ranges; individual layers can override data, mappings and styles.

### Horizontal, vertical and circular layouts

Use the same data and track definitions in all three layouts. `reverse_chr` changes the display direction of selected chromosomes while preserving source coordinates.

```r
tracks <- list(
  density = geom_track(side = "right", width = 3,
    data = density, mapping = aes(chr = Chr, x = Mid, y = Rate),
    geom = geom_line(colour = "#4477AA", linewidth = 0.4))
)

p_vertical <- ggideogram(human_karyotype, chr = c("1", "2"),
  orientation = "vertical", tracks = tracks)
p_horizontal <- ggideogram(human_karyotype, chr = c("1", "2"),
  orientation = "horizontal", reverse_chr = "2", tracks = tracks)
p_circular <- ggideogram(human_karyotype,
  orientation = "circular", opening_angle = 25, tracks = tracks)
p_circular
```

`side = "left"` means left in a vertical layout, above in a horizontal layout, and inward in a circular layout; `"right"` selects the opposite side. Circular layouts also accept `"inner"` / `"outer"`. Use `"overlay"` for tracks inside chromosome bodies. Track `width` and `gap` are measured in chromosome body widths; native point sizes, text sizes and line widths retain ggplot2's physical sizing.

### Adjusting a track

Reuse the same `track` name to update a track. Omitted settings stay unchanged; new layers are appended unless `replace = TRUE` is set.

Ordinary beside tracks default to a width of 3 chromosome body widths; gene-model tracks default to 6, and overlay tracks to 1. Explicit widths take priority.

For example, `p + geom_track(track = "density", width = 4, gap = 0.5)` adjusts space, `label = list(size = 2.8)` changes the title size, and `axis = list(n = 2)` adjusts the number of value-axis ticks.

Track tick text defaults to 2.2 mm, and `theme_ideogram()` provides a 6 mm
legend gap. Circular plots may reserve more space for outward text. Both
settings can be adjusted independently of track geometry:

```r
p + geom_track(track = "density", axis = list(size = 2.2, label_gap = 0.35)) +
  theme(legend.position = "bottom", legend.box.spacing = grid::unit(6, "mm"))
```

`limits`, `value_scale`, `transform` and `reverse` control the track's value scale. Set `label = NULL` to remove its title or `axis = FALSE` to hide its value axis. Native point, line, col, area, ribbon, tile, text, boxplot and violin layers enter through `geom` or `layers`. Boxplot and violin statistics run on raw values, with widths specified in bp. Compatible third-party identity geoms can use the same interface; see the [interactive track examples](inst/examples/interactive-gallery.R).

### Choosing the right setting

| Change | Setting | Meaning / unit |
| --- | --- | --- |
| Visible genomic interval | `chr_view(start, end)` | Source bp boundaries; annotations retain their original coordinates |
| Chromosome arrangement | `ggideogram(orientation, ncol, reverse_chr)` | Display direction and grid arrangement |
| Plot proportions | `max_chr_length`, `chromosome_gap`, `row_gap`; circular `radius`, `opening_angle` | Layout dimensions; opening angle in degrees |
| Track allocation | `geom_track(width, gap, side)` | Chromosome body widths |
| Quantitative comparison | `limits`, `value_scale`, `transform`, `reverse` | Raw value range, sharing, transformation and display direction |
| Track title / ticks | `label = list(...)`, `axis = list(...)` | Text size in mm; title/label gaps in em; tick length in mm |
| Marks and lines | Native `geom_*(size, linewidth, alpha, ...)` | Standard ggplot2 styling |
| Gene text | `geom_genemodel(label_size, label_family, label_fontface, label_colour)` | Text styled independently of exon and intron geometry |
| Legends / general text | `theme()`, `guides()`, `scale_*()` | Theme text in pt; continuous and discrete native guides |
| Export slot | `ggsave(width, height, units)` | Physical output dimensions; native text and mark sizes stay fixed |

For comparisons, use shared limits across the relevant tracks or chromosomes.
Per-chromosome ranges emphasize local variation and should be labelled explicitly.
Automatic ranges include native tile footprints and deterministic nudge, stack
and fill adjustments. Missing observations remain missing; a completely missing
layer can share limits with another layer or use explicit limits.

Use a linear view for precise local positions and gene structures, circular
views for chromosome relationships, and chromosome x/y scales for summaries.
Choose the output size first, then allocate tracks and labels for that slot.
The gallery uses vector PDF output and separates explanatory captions from
the plotting area; PNG files provide previews.

### Extension boundaries

| Input | Interface | Behaviour |
| --- | --- | --- |
| Native / compatible third-party identity geom | `geom_track(geom = ...)` or `layers = list(...)` | Native Geom and Position with bp/value projection |
| Area, ribbon, boxplot or violin | Same track interface | Supported native statistics run before projection |
| A different statistical summary | Summarize first, or place its complete plot with `geom_locus_inset()` | An arbitrary Stat needs an explicit adapter |
| Independent colour scales | `ggnewscale::new_scale_*()` between layers | Native scale ownership survives track updates |
| Complete ggplot or grob | `geom_locus_inset()` | Independent coordinates, theme, guides and physical viewport |
| Ideogram as a panel / inset | patchwork, cowplot, `as_ideogram_grob()` | Standard ggplot/grob composition |

The CI matrix targets ggplot2 3.5.0 with compatible versions of ggrepel, patchwork and
cowplot, and the current ggplot2 release with current optional extensions.
The ggiraph feature requires a mutually compatible ggiraph/ggplot2 pair;
ggiraph 0.9.0 requires ggplot2 ≥ 3.5.2, and 0.9.6 requires ≥ 4.0.0.

## Local views and gene models

This example uses bundled Arabidopsis TAIR10 / Araport11 annotations. It resolves GFF3 Parent relationships between transcripts and genes, expands shared features with multiple parents, then draws a local view of Chr1.

```r
features <- read_chr_features(system.file("extdata",
  "arabidopsis-first-genes.gff3", package = "ggideogram"))
tx <- features[features$Type == "mRNA", ]
parents <- strsplit(ifelse(features$Type == "mRNA",
  features$ID, features$Parent), ",", fixed = TRUE)
models <- features[rep(seq_len(nrow(features)), lengths(parents)), ]
models$Transcript <- unlist(parents, use.names = FALSE)
models <- models[models$Transcript %in% tx$ID, ]
parents <- strsplit(tx$Parent[match(models$Transcript, tx$ID)], ",", fixed = TRUE)
models <- models[rep(seq_len(nrow(models)), lengths(parents)), ]
models$Gene <- sub("gene:", "", unlist(parents, use.names = FALSE))

arabidopsis <- data.frame(Chr = c("1", "5"), Start = 0,
  End = c(30427671, 26975502))
view <- chr_view(arabidopsis, chr = "1", start = 2500, end = 15000)

p_genes <- ggideogram(view, orientation = "horizontal", axis = TRUE) +
  geom_track(track = "genes", side = "right", width = 5,
    data = models,
    mapping = aes(chr = Chr, start = Start, end = End,
      gene = Gene, type = Type, strand = Strand, fill = Strand),
    layers = list(geom_genemodel(mode = "gene", label_size = 2.8))) +
  scale_fill_manual(values = c("+" = "#4477AA", "-" = "#CC6677")) +
  theme(plot.margin = margin(5, 5, 5, 25, unit = "mm"))
p_genes
```

Exon outlines show complete exon intervals, CDS and explicitly annotated UTRs use fills, and intron arrows indicate strand direction. `mode = "gene"` shows the union of each feature type within a gene; `mode = "transcript"` with a `transcript` mapping draws individual transcripts. Local tracks accept complete source annotations, clip structures to the display window, and retain original bp labels.

Use `label_fontface = "italic"` for gene symbols when appropriate for the
organism's nomenclature. `label_family = NULL` inherits the plot's
`base_family`; `label_colour` and `label_size` control the text independently
of the feature fills and intron lines.

`read_karyotype()` reads chrom.sizes / FAI files. `read_chr_features()` reads BED / GFF3 / GTF, and `as_chr_features()` accepts data frames or GRanges. **The BED reader converts coordinates to 1-based closed intervals; GFF3, GTF and GRanges retain their annotation coordinates.** Chromosome names and assembly versions should match the karyotype.

For window summaries, `bin_genome()` provides feature counts (`count`), union coverage fractions (`coverage`), and means weighted by overlap length (`weighted_mean`). Map the results directly to quantitative tracks.

Counts and `FUN` summaries accept `count_position = "midpoint"`, `"start"` or `"end"`. The default midpoint is `(Start + End) / 2`; a row without `End` uses `Start`. `GFFex()` selects one GFF feature type and uses the same engine with `count_position = "start"`, retaining its feature-start definition. Coverage and weighted means use complete interval overlaps.

All summary results carry `attr(result, "window_summary")` with the method, source value column, coordinate convention, nominal window width, representative position, missing-value argument, overlap rule and empty fill. Count units are features and coverage units are fractions; numerical summaries record their source value column with unspecified units. The `FUN` entry retains the four columns `Chr`, `Start`, `End`, `Value`. Explicit methods also return actual `Width`, `N` and `N_valid`; counts include `Rate` per Mb using each window's actual width.

## Locus annotations and links

`geom_locus()` adds points, text, intervals and leaders. Text accepts `position = "identity"`, `"spread"` or `"repel"`. Ordered spreading uses actual text dimensions while leaders retain the source loci. `track` and `track_position` place annotations at a declared track's near edge, middle or far edge. See the [label examples](inst/examples/label-layout-gallery.R) for complete code.

The following plot shows Chr1–Chr5 collinear blocks from the official MCScanX Arabidopsis example.

```r
blocks <- read.delim(system.file("extdata",
  "arabidopsis-synteny-blocks.tsv", package = "ggideogram"))
blocks <- blocks[blocks$Chr1 != blocks$Chr2, ]

p_links <- ggideogram(arabidopsis, orientation = "circular") +
  geom_chrlink(type = "interval", data = blocks,
    mapping = aes(chr1 = Chr1, start1 = Start1, end1 = End1,
      chr2 = Chr2, start2 = Start2, end2 = End2,
      orientation = Orientation, fill = Orientation),
    alpha = 0.4, colour = NA) +
  scale_fill_manual(values = c("+" = "#4477AA", "-" = "#CC6677"),
    name = "Block order")
p_links
```

`type = "point"` connects two loci; `type = "interval"` preserves interval widths at both ends. You can also build a node table with unique IDs using `chr_nodes()`, then connect nodes through `aes(from, to)`. For multiple genomes, `chr_key()` distinguishes genome, chromosome and assembly identities; see the [multi-genome examples](inst/examples/multi-genome-gallery.R).

## Insets and plot composition

A complete inset retains its own coordinates, theme and guides. `geom_locus_inset()` anchors it to the host plot's bp coordinates.

```r
child <- ggplot(subset(density, Chr == "1"), aes(Value)) +
  geom_histogram(bins = 15, fill = "#4477AA", colour = "white") +
  labs(x = "Genes / window", y = "Windows") +
  theme_minimal(base_size = 8)

p_inset <- ggideogram(human_karyotype, chr = "1",
  orientation = "horizontal",
  tracks = list(summary = geom_track(side = "right", width = 14))) +
  geom_locus_inset(plot = child,
    data = data.frame(Chr = "1", Pos = 1.2e8),
    mapping = aes(chr = Chr, position = Pos), track = "summary",
    width = grid::unit(40, "mm"), height = grid::unit(30, "mm"))
p_inset

# Use the ideogram as a panel or embed it in another plot.
p_combined <- patchwork::wrap_plots(p, child, widths = c(2, 1))
p_embedded <- child + patchwork::inset_element(p_circular,
  left = 0.5, bottom = 0.5, right = 1, top = 1)
p_combined

ggsave("ideogram-combined.pdf", p_combined,
  width = 180, height = 100, units = "mm")
```

Chromosomes can also form the native discrete x/y axis of an ordinary ggplot, aligning with its layers by chromosome key.

```r
p_axis <- ggplot(subset(gene_density, Chr %in% c("1", "2", "3")),
  aes(Chr, Value)) +
  geom_boxplot(fill = "#AACCEE") +
  scale_x_chromosome(human_karyotype, limits = NULL) +
  labs(x = NULL, y = "Genes / window") +
  theme_minimal()
p_axis
```

`limits = NULL` trains the axis on the plotted chromosomes. The default keeps every chromosome in the supplied karyotype; `drop` applies when limits are trained.

## Public functions

| Purpose | Functions |
| --- | --- |
| Build plots; draw chromosome bodies, names, ticks, bands and internal fills | `ggideogram()`, `geom_chr(component = ...)` |
| Declare or adjust tracks and add native layers | `geom_track()`, `track_layout()` |
| Gene and transcript structures | `geom_genemodel(mode = ...)` |
| Locus points, text, intervals and leaders | `geom_locus(geom = ...)` |
| Node tables and paired relationships | `chr_nodes()`, `geom_chrlink(type = ...)` |
| Anchor complete insets or extract a standard grob | `geom_locus_inset()`, `as_ideogram_grob()` |
| Chromosome axes for ordinary plots | `scale_x_chromosome()`, `scale_y_chromosome()` |
| Read data, select local views and summarize windows | `read_karyotype()`, `read_chr_features()`, `chr_view()`, `view_chr_data()`, `bin_genome()` |
| Coordinate projections for extensions | `project_chr_point()`, `project_chr_interval()`, `project_chr_track()` |
| Interactive viewing, selection and source-record export | `as_ideogram_widget()` |

See individual R help pages, such as `?geom_track`, for parameters and examples. Use standard `theme()`, `scale_*()` and `guides()` to style plots.

## Example gallery

| Content | Runnable scripts / figures |
| --- | --- |
| Local linear and circular views, shared tracks and gene models | [Script](inst/examples/optimized-gallery.R) · [Linear PNG](inst/examples/optimized-local-linear.png) · [Linear PDF](inst/examples/optimized-local-linear.pdf) · [Circular PNG](inst/examples/optimized-local-circle.png) |
| Ideograms combined with ordinary bar charts and boxplots | [Script](inst/examples/original-data-composition.R) · [PNG](inst/examples/original-data-composition.png) |
| Complete plots inside ideograms, and ideograms inside other plots | [Script](inst/examples/original-data-gallery.R) · [PNG](inst/examples/gallery-bidirectional-insets.png) |
| Native chromosome x/y axes | [Script](inst/examples/original-data-axis-integration.R) · [x-axis figure](inst/examples/original-data-chromosome-x-axis.png) · [y-axis figure](inst/examples/original-data-chromosome-y-axis.png) |
| Local windows, genes and transcripts | [Local views](inst/examples/local-view-gallery.R) · [Gene models](inst/examples/gene-structure-gallery.R) |
| Internal fills, dense labels, synteny and multiple genomes | [Internal annotations](inst/examples/internal-annotation-gallery.R) · [Labels](inst/examples/label-layout-gallery.R) · [Links](inst/examples/connection-gallery.R) · [Multiple genomes](inst/examples/multi-genome-gallery.R) |
| Whole-genome circular and interactive plots | [Circular](inst/examples/circular-gallery.R) · [Interactive](inst/examples/interactive-gallery.R) |

Each script starts with run instructions. Associated data are bundled in [inst/extdata](inst/extdata); source records and preparation details appear in script comments and provenance files in that directory.

## Acknowledgements and license

This project's chromosome geometry and some example datasets derive from [RIdeogram](https://github.com/TickingClock1992/RIdeogram). We thank its original authors: Zhaodong Hao, Dekang Lv, Ying Ge, Jisen Shi, Dolf Weijers, Guangchuang Yu and Jinhui Chen. Original publication: Hao et al. (2020), *RIdeogram: drawing SVG graphics to visualize and map genome-wide data on the idiograms*, [PeerJ Computer Science 6:e251](https://doi.org/10.7717/peerj-cs.251).

`ggideogram` is an independent ggplot2 extension with a reworked data model, layout, native layers and composition interfaces, distributed under the [Artistic License 2.0](LICENSE). Other bundled datasets include source records; the license for MCScanX materials is provided in [MCScanX-LICENSE.txt](inst/extdata/MCScanX-LICENSE.txt).
