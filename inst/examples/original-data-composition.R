# Composition example built from the four datasets bundled with ggideogram.
#
# From a source checkout:
#   Rscript -e 'pkgload::load_all("."); source("inst/examples/original-data-composition.R")'

suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(patchwork)
})

# Device and typography are output choices, not geometry constants. Change
# these named values without changing any chromosome or track calculation.
figure <- list(
  width_mm = 280,
  height_mm = 190,
  dpi = 200,
  panel_widths = c(1.2, 1),
  base_size = 9,
  marker_size = 0.65,
  marker_stroke = 0.18,
  marker_separation = 0.18,
  ideogram_margin_pt = c(top = 5.5, right = 5.5, bottom = 5.5, left = 18),
  output_dir = file.path("inst", "examples"),
  output_stem = "original-data-composition"
)

data("human_karyotype", package = "ggideogram")
data("gene_density", package = "ggideogram")
data("LTR_density", package = "ggideogram")
data("Random_RNAs_500", package = "ggideogram")

gene_density$Position <- (gene_density$Start + gene_density$End) / 2
LTR_density$Position <- (LTR_density$Start + LTR_density$End) / 2
Random_RNAs_500$Position <-
  (Random_RNAs_500$Start + Random_RNAs_500$End) / 2

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

tracks <- track_layout(
  markers = track(side = "right", width = 0.65, gap = 0.20),
  genes = track(
    side = "right", width = 1.15, gap = 0.25,
    limits = range(gene_density$Value, na.rm = TRUE)
  ),
  ltr = track(
    side = "right", width = 1.15, gap = 0.25,
    limits = range(LTR_density$Value, na.rm = TRUE)
  )
)

# Choose the chromosome grid from the data, declared tracks and actual slot
# aspect. Restricting candidates to exact divisors avoids anonymous empty cells
# in the final row for this complete-genome example.
semantic <- as_ideogram_data(human_karyotype)
chromosome_count <- nrow(semantic$karyotype)
column_candidates <- seq_len(chromosome_count)
column_candidates <- column_candidates[
  chromosome_count %% column_candidates == 0
]
target_aspect <-
  figure$width_mm * figure$panel_widths[1] / sum(figure$panel_widths) /
  figure$height_mm
candidate_aspect <- vapply(column_candidates, function(columns) {
  candidate_layout <- ideogram_layout(
    semantic, ncol = columns, tracks = tracks
  )
  diff(candidate_layout$bounds$x) / diff(candidate_layout$bounds$y)
}, numeric(1))
ideogram_columns <- column_candidates[
  which.min(abs(log(candidate_aspect / target_aspect)))
]

# One visible bp ruler per chromosome row documents the longitudinal alignment
# shared by both side tracks without repeating colliding labels 24 times. The
# transverse value scales stay local to their declared tracks; drawing 24 tiny
# copies of those axes would be less informative than the ordinary summary
# panels B and C.
row_axis_chromosomes <- human_karyotype$Chr[
  seq.int(1, chromosome_count, by = ideogram_columns)
]

repel <- position_chr_repel(
  figure$marker_separation, units = "body_width"
)
panel_a <- ggideogram(
  human_karyotype,
  ncol = ideogram_columns,
  tracks = tracks,
  axis = row_axis_chromosomes,
  axis_side = "left",
  base_family = "sans"
) +
  geom_track_line(
    data = gene_density,
    mapping = aes(
      chr = Chr, position = Position, value = Value, group = Chr
    ),
    track = "genes",
    colour = "#277DA1",
    linewidth = 0.25
  ) +
  geom_track_col(
    data = LTR_density,
    mapping = aes(
      chr = Chr, position = Position, value = Value,
      width = End - Start
    ),
    track = "ltr",
    position = "identity",
    fill = "#43AA8B",
    alpha = 0.70
  ) +
  geom_chr_link(
    data = Random_RNAs_500,
    mapping = aes(chr = Chr, position = Position),
    position = repel,
    track = "markers",
    colour = "#777777",
    linewidth = 0.12,
    alpha = 0.45
  ) +
  geom_chr_marker(
    data = Random_RNAs_500,
    mapping = aes(
      chr = Chr, position = Position, shape = Type, fill = Type
    ),
    position = repel,
    track = "markers",
    size = figure$marker_size,
    stroke = figure$marker_stroke,
    colour = "#303030"
  ) +
  scale_shape_manual(values = rna_shapes) +
  scale_fill_manual(values = rna_colours) +
  labs(
    title = "Human chromosome ideogram",
    subtitle = "bp-aligned gene and LTR tracks; native ggplot2 RNA markers",
    shape = "RNA type",
    fill = "RNA type"
  ) +
  guides(
    shape = guide_legend(title.position = "top", nrow = 1),
    fill = guide_legend(title.position = "top", nrow = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = figure$base_size),
    legend.text = element_text(size = figure$base_size),
    plot.title = element_text(face = "bold", size = figure$base_size + 1),
    plot.subtitle = element_text(size = figure$base_size),
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
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = rna_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "RNA marker counts", x = NULL, y = "Records") +
  common_theme

gene_density$Chr <- factor(
  gene_density$Chr, levels = rev(human_karyotype$Chr)
)
panel_c <- ggplot(gene_density, aes(Value, Chr)) +
  geom_boxplot(
    width = 0.62,
    outlier.size = 0.35,
    linewidth = 0.28,
    fill = "#DCEAF3",
    colour = "#315A7D"
  ) +
  labs(
    title = "Gene density by chromosome",
    x = "Genes per 1 Mb window",
    y = NULL
  ) +
  common_theme

right_column <- panel_b / panel_c + plot_layout(heights = c(0.65, 1.35))
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
