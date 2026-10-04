# Composition example built from the four datasets bundled with ggideogram.
#
# From a source checkout:
#   Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/original-data-composition.R")'

suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(patchwork)
})
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))

# Device and typography are output choices, not geometry constants. Change
# these named values without changing any chromosome or track calculation.
figure <- list(
  width_mm = 240,
  height_mm = 185,
  dpi = 200,
  panel_widths = c(2.2, 1),
  base_size = 9,
  row_gap = 5,
  marker_size = 0.65,
  marker_stroke = 0.18,
  marker_separation = 0.18,
  ideogram_margin_pt = c(top = 5.5, right = 5.5, bottom = 5.5, left = 18),
  output_dir = file.path("work", "api-optimization", "examples", "original-data-composition"),
  output_stem = "original-data-composition"
)

data("human_karyotype", package = "ggideogram")
data("gene_density", package = "ggideogram")
data("LTR_density", package = "ggideogram")
data("Random_RNAs_500", package = "ggideogram")

gene_density$Position <- (gene_density$Start - 1 + gene_density$End) / 2
gene_density$Width <- gene_density$End - gene_density$Start + 1
gene_density$Rate <- gene_density$Value / (gene_density$Width / 1e6)
LTR_density$Position <- (LTR_density$Start - 1 + LTR_density$End) / 2
Random_RNAs_500$Position <-
  (Random_RNAs_500$Start - 1 + Random_RNAs_500$End) / 2

# Original marker names are translated once to standard filled ggplot2 shapes.
# The renderer itself has no shape registry or dataset-specific branch.
shape_lookup <- c(circle = 21, box = 22, triangle = 24)
rna_key <- unique(Random_RNAs_500[c("Type", "Shape", "color")])
unknown_shapes <- setdiff(rna_key$Shape, names(shape_lookup))
if (length(unknown_shapes)) {
  stop("No standard ggplot2 shape mapping for: ",
       paste(unknown_shapes, collapse = ", "))
}
rna_shapes <- stats::setNames(shape_lookup[rna_key$Shape], rna_key$Type)
rna_colours <- stats::setNames(paste0("#", rna_key$color), rna_key$Type)

tracks <- track_layout(markers = geom_track(side = "right", width = 0.65, gap = 0.2), genes = geom_track(side = "right",
    width = 3, gap = 0.5, limits = range(gene_density$Value, na.rm = TRUE)), ltr = geom_track(side = "right",
    width = 3, gap = 0.5, limits = range(LTR_density$Value, na.rm = TRUE)))

# Choose the chromosome grid from the data, declared tracks and actual slot
# aspect. Restricting candidates to exact divisors avoids anonymous empty cells
# in the final row for this complete-genome example.
semantic <- as_ideogram_data(human_karyotype)
chromosome_count <- nrow(semantic$karyotype)
ideogram_columns <- gallery_grid(semantic, tracks,
  figure$width_mm * figure$panel_widths[1] / sum(figure$panel_widths) - 20,
  figure$height_mm - 42, row_gap = figure$row_gap)$ncol

# One visible bp ruler per chromosome row documents the longitudinal alignment
# shared by both side tracks. Their common value ranges are stated below panel A.
row_axis_chromosomes <- human_karyotype$Chr[
  seq.int(1, chromosome_count, by = ideogram_columns)
]

repel <- position_chr_repel(
  figure$marker_separation, units = "body_width"
)
panel_a <- ggideogram(
  human_karyotype,
  ncol = ideogram_columns,
  row_gap = figure$row_gap,
  tracks = tracks,
  axis = row_axis_chromosomes,
  axis_breaks = seq(0, 250e6, 100e6),
  axis_side = "left",
  base_family = "sans"
) +
  geom_track(data = gene_density, mapping = aes(chr = Chr, x = Position, y = Value, group = Chr),
      track = "genes", geom = ggplot2::geom_line(colour = "#277DA1", linewidth = 0.25)) +
  geom_track(data = LTR_density, mapping = aes(chr = Chr, x = Position, y = Value, width = End -
      Start + 1), track = "ltr", geom = ggplot2::geom_col(position = "identity", fill = "#43AA8B",
      alpha = 0.7)) +
  geom_locus(geom = "link", data = Random_RNAs_500, mapping = aes(chr = Chr, position = Position),
      position = repel, track = "markers", colour = "#777777", linewidth = 0.12, alpha = 0.45) +
  geom_locus(geom = "point", data = Random_RNAs_500, mapping = aes(chr = Chr, position = Position,
      shape = Type, fill = Type), position = repel, track = "markers", size = figure$marker_size,
      stroke = figure$marker_stroke, colour = "#303030") +
  scale_shape_manual(values = rna_shapes) +
  scale_fill_manual(values = rna_colours) +
  labs(
    title = "Human chromosomes · GRCh38",
    subtitle = "Genes, LTRs and sampled RNAs",
    caption = sprintf(paste0("Counts/window; nominal 1 Mb bins.\n",
      "Blue: genes (%g-%g); green: LTRs (%g-%g).\n",
      "Terminal windows are shorter."),
      min(gene_density$Value), max(gene_density$Value),
      min(LTR_density$Value), max(LTR_density$Value)),
    shape = "RNA type",
    fill = "RNA type"
  ) +
  guides(
    shape = guide_legend(title.position = "top", nrow = 1, override.aes = list(size = 1.8)),
    fill = guide_legend(title.position = "top", nrow = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = gallery_legend_spacing,
    legend.box = "horizontal",
    legend.title = element_text(size = figure$base_size),
    legend.text = element_text(size = figure$base_size),
    plot.title = element_text(face = "bold", size = figure$base_size + 1),
    plot.subtitle = element_text(size = figure$base_size),
    plot.caption = element_text(size = figure$base_size - 1, hjust = 0),
    plot.margin = do.call(
      margin,
      as.list(unname(figure$ideogram_margin_pt))
    )
  )

rna_counts <- as.data.frame(table(factor(
  Random_RNAs_500$Type, levels = rna_key$Type
)))
names(rna_counts) <- c("Type", "Count")

common_theme <- theme_minimal(base_size = figure$base_size) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot"
  )

panel_b <- ggplot(rna_counts, aes(Type, Count, fill = Type)) +
  geom_col(width = 0.65, alpha = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = rna_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Sampled RNA records", x = NULL, y = "Records") +
  common_theme

gene_density$Chr <- factor(
  gene_density$Chr, levels = rev(human_karyotype$Chr)
)
panel_c <- ggplot(gene_density, aes(Rate, Chr)) +
  geom_boxplot(
    width = 0.62,
    outlier.size = 0.35,
    linewidth = 0.28,
    fill = "#DCEAF3",
    colour = "#315A7D"
  ) +
  labs(
    title = "Gene density",
    x = "Genes / Mb",
    y = NULL
  ) +
  common_theme

right_column <- panel_b / panel_c + plot_layout(heights = c(0.4, 1.6))
composition <- panel_a | right_column
composition <- composition +
  plot_layout(widths = figure$panel_widths) +
  plot_annotation(tag_levels = "A")

dir.create(figure$output_dir, recursive = TRUE, showWarnings = FALSE)
png_file <- file.path(
  figure$output_dir, paste0(figure$output_stem, ".png")
)
pdf_file <- file.path(
  figure$output_dir, paste0(figure$output_stem, ".pdf")
)

ggsave(
  png_file, composition,
  width = figure$width_mm, height = figure$height_mm,
  units = "mm", dpi = figure$dpi, bg = "white"
)
ggsave(
  pdf_file, composition,
  width = figure$width_mm, height = figure$height_mm,
  units = "mm", bg = "white"
)

message("Wrote ", normalizePath(png_file))
message("Wrote ", normalizePath(pdf_file))
