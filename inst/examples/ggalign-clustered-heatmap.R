# Additional packages: ggalign >= 1.2.0, svglite and rsvg.
# From a source checkout:
# Rscript -e 'pkgload::load_all("."); source("inst/examples/ggalign-clustered-heatmap.R")'
suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(ggalign)
})

data("human_karyotype", package = "ggideogram")
data("gene_density", package = "ggideogram")
data("LTR_density", package = "ggideogram")
karyotype <- human_karyotype

# Assign each source window to a relative-position bin by its midpoint.
# Pool counts and observed window widths before converting to records per Mb.
relative_profile <- function(data, karyotype, bins = 10) {
  length_bp <- karyotype$End[match(as.character(data$Chr), as.character(karyotype$Chr))]
  midpoint <- (data$Start - 1 + data$End) / 2
  bin <- pmin(bins, floor(midpoint / length_bp * bins) + 1L)
  pooled <- aggregate(
    cbind(Count = data$Value, Width = data$End - data$Start + 1),
    list(Chr = as.character(data$Chr), Bin = bin), sum)
  edges <- seq(0, 100, length.out = bins + 1L)
  result <- matrix(NA_real_, nrow = nrow(karyotype), ncol = bins,
    dimnames = list(as.character(karyotype$Chr),
      sprintf("%02d-%02d", head(edges, -1), tail(edges, -1))))
  result[cbind(match(pooled$Chr, rownames(result)), pooled$Bin)] <-
    pooled$Count / pooled$Width * 1e6
  result
}

genes <- relative_profile(gene_density, karyotype)
ltr <- relative_profile(LTR_density, karyotype)
density <- cbind(genes, ltr)
colnames(density) <- paste(rep(c("Genes", "LTR"), each = ncol(genes)),
  rep(colnames(genes), 2), sep = ":")

# Standardize each feature/bin across chromosomes and compute the tree once.
z <- scale(density)
tree <- hclust(dist(z, method = "euclidean"), method = "complete")
z_limit <- ceiling(max(abs(z)))

totals <- aggregate(
  cbind(Count = gene_density$Value, Width = gene_density$End - gene_density$Start + 1),
  list(Chr = as.character(gene_density$Chr)), sum)
gene_rate <- matrix(totals$Count / totals$Width * 1e6, ncol = 1,
  dimnames = list(totals$Chr, "Genes"))
gene_rate <- gene_rate[rownames(z), , drop = FALSE]

base_theme <- theme_classic(base_size = 7, base_family = "Arial") + theme(
  text = element_text(colour = "black"),
  line = element_line(colour = "black", linewidth = .225),
  axis.text = element_text(colour = "black"),
  axis.ticks.length = grid::unit(1.2, "mm"),
  strip.background = element_rect(fill = "white", colour = "black", linewidth = .225),
  strip.text = element_text(colour = "black", size = 7),
  legend.title = element_text(size = 6), legend.text = element_text(size = 6),
  legend.background = element_blank(), legend.key = element_blank(),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = grid::unit(1.5, "mm"),
  plot.margin = margin(0, 0, 0, 0))
matrix_theme <- base_theme + theme(
  axis.line = element_blank(), axis.ticks = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.spacing.x = grid::unit(2, "mm"))
tree_theme <- base_theme + theme(
  axis.text = element_blank(), axis.title = element_blank(),
  axis.line = element_blank(), axis.ticks = element_blank())

# One native heatmap owns the shared row domain. Its chromosome guide, tree
# and summary track all follow the same chromosome keys after clustering.
figure <- ggheatmap(z, filling = NULL,
    width = grid::unit(100, "mm"), height = grid::unit(124, "mm")) +
  geom_tile(aes(x = .x, y = .discrete_y, fill = value),
    colour = "white", linewidth = .1) +
  scale_fill_gradient2(low = "#4477AA", mid = "#F7F7F7", high = "#EE6677",
    midpoint = 0, limits = c(-z_limit, z_limit), breaks = seq(-z_limit, z_limit, 2),
    name = "Z-score", guide = guide_colourbar(theme = theme(
      legend.key.width = grid::unit(1.8, "mm"),
      legend.key.height = grid::unit(20, "mm"),
      legend.title = element_text(margin = margin(b = 1.5, unit = "mm"))))) +
  scale_y_chromosome(karyotype, limits = NULL,
    labels = function(chr) paste0("chr", chr),
    axis_length = grid::unit(24, "mm"), body_width = grid::unit(1.5, "mm"),
    linewidth = .225) +
  scale_x_continuous(labels = function(key) sub(".*:", "", key)) +
  labs(x = "Chromosome position (%)", y = NULL) + matrix_theme +
  anno_top() + align_group(rep(c("Genes", "LTR"), each = ncol(genes))) +
  anno_left(size = grid::unit(22, "mm")) +
  align_dendro(distance = NULL, method = tree, linewidth = .225) + tree_theme +
  anno_right(size = grid::unit(23, "mm")) + ggalign(data = gene_rate) +
  geom_col(aes(x = value), orientation = "y", width = .65, fill = "#4477AA") +
  scale_x_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "Genes / Mb", y = NULL) + base_theme + theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.y = element_blank(), plot.margin = margin(l = 2, unit = "mm")) +
  quad_active() + layout_theme(
    legend.position = "right", panel.spacing.x = grid::unit(2, "mm"),
    plot.margin = margin(3, 3, 3, 3, unit = "mm"))

output_dir <- file.path("work", "ggalign-examples")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
svg <- file.path(output_dir, "ggalign-clustered-heatmap.svg")
png <- file.path(output_dir, "ggalign-clustered-heatmap.png")
save_ideogram(svg, figure, width = 210, height = 160, units = "mm",
  device = svglite::svglite, bg = "white")
svg_text <- readLines(svg, warn = FALSE)
svg_text <- append(svg_text,
  '<metadata xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:creator>zbh</dc:creator></metadata>',
  after = 2)
writeLines(svg_text, svg)
rsvg::rsvg_png(svg, png, width = 210 / 25.4 * 300, height = 160 / 25.4 * 300)
