# Native chromosome-axis examples built only from data bundled with
# ggideogram. Run from a source checkout after loading the package:
#
#   Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/original-data-axis-integration.R")'

suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
})

# Every device and appearance choice is named here and can be changed without
# touching the axis engine or the biological data.
axis_gallery <- list(
  output_dir = file.path("work", "api-optimization", "examples", "original-data-axis-integration"),
  dpi = 300,
  base_size = 9,
  body_fill = "#FAFAFA",
  body_colour = "black",
  body_linewidth = 0.25,
  plot_margin_pt = c(top = 7, right = 8, bottom = 7, left = 8),
  axis_text_size = 9,
  x_axis = list(
    width_mm = 185,
    height_mm = 100,
    axis_length_mm = 16,
    body_width_mm = 1.5,
    label_gap_mm = 0.8,
    bar_width = 0.72,
    bar_fill = "#2C7FB8",
    value_expand = c(lower = 0, upper = 0.06)
  ),
  y_axis = list(
    width_mm = 190,
    height_mm = 150,
    axis_length_mm = 22,
    body_width_mm = 1.45,
    label_gap_mm = 0.8,
    box_width = 0.58,
    box_fill = "#DCEAF4",
    box_colour = "#356A92",
    box_linewidth = 0.36,
    outlier_shape = 16,
    outlier_size = 0.55,
    outlier_alpha = 0.72,
    value_expand = c(lower = 0.01, upper = 0.04)
  )
)

data("human_karyotype", package = "ggideogram")
data("gene_density", package = "ggideogram")
data("Random_RNAs_500", package = "ggideogram")
gene_density$Width <- gene_density$End - gene_density$Start + 1
gene_density$Rate <- gene_density$Value / (gene_density$Width / 1e6)

chromosome_order <- as.character(human_karyotype$Chr)
if (any(!Random_RNAs_500$Chr %in% chromosome_order) ||
    any(!gene_density$Chr %in% chromosome_order)) {
  stop("Bundled observations contain a chromosome absent from human_karyotype.")
}

# table() retains zero-count chromosomes because the factor levels are the
# complete bundled karyotype order.
rna_counts <- as.data.frame(table(factor(
  Random_RNAs_500$Chr,
  levels = chromosome_order
)), stringsAsFactors = FALSE)
names(rna_counts) <- c("Chr", "Records")
rna_counts$Chr <- as.character(rna_counts$Chr)
stopifnot(sum(rna_counts$Records) == nrow(Random_RNAs_500))

save_axis_example <- function(plot, stem, dimensions) {
  dir.create(axis_gallery$output_dir, recursive = TRUE, showWarnings = FALSE)
  png_file <- file.path(axis_gallery$output_dir, paste0(stem, ".png"))
  pdf_file <- file.path(axis_gallery$output_dir, paste0(stem, ".pdf"))
  ggsave(
    png_file, plot,
    width = dimensions$width_mm,
    height = dimensions$height_mm,
    units = "mm", dpi = axis_gallery$dpi, bg = "white"
  )
  ggsave(
    pdf_file, plot,
    width = dimensions$width_mm,
    height = dimensions$height_mm,
    units = "mm", bg = "white"
  )
  message("Wrote ", normalizePath(png_file))
  message("Wrote ", normalizePath(pdf_file))
  invisible(c(png = png_file, pdf = pdf_file))
}

publication_theme <- theme_classic(base_size = axis_gallery$base_size) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.ticks = element_line(colour = "black", linewidth = 0.25),
    plot.title.position = "plot",
    plot.margin = do.call(margin, as.list(unname(axis_gallery$plot_margin_pt)))
  )

# The chromosome guide is the x axis of an otherwise ordinary geom_col plot.
x_axis_plot <- ggplot(rna_counts, aes(Chr, Records)) +
  geom_col(
    width = axis_gallery$x_axis$bar_width,
    fill = axis_gallery$x_axis$bar_fill
  ) +
  scale_x_chromosome(
    human_karyotype,
    limits = chromosome_order,
    axis_length = grid::unit(axis_gallery$x_axis$axis_length_mm, "mm"),
    body_width = grid::unit(axis_gallery$x_axis$body_width_mm, "mm"),
    label_gap = grid::unit(axis_gallery$x_axis$label_gap_mm, "mm"),
    fill = axis_gallery$body_fill,
    colour = axis_gallery$body_colour,
    linewidth = axis_gallery$body_linewidth
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = axis_gallery$x_axis$value_expand)
  ) +
  labs(
    x = "Chromosome",
    y = "Sampled records"
  ) +
  publication_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = axis_gallery$axis_text_size)
  )

# The chromosome guide is the y axis of an otherwise ordinary geom_boxplot
# plot. All original 1 Mb windows remain in the input; boxplot outliers are
# rendered as small native ggplot2 points rather than suppressed.
y_axis_plot <- ggplot(gene_density, aes(Rate, Chr)) +
  geom_boxplot(
    orientation = "y",
    width = axis_gallery$y_axis$box_width,
    fill = axis_gallery$y_axis$box_fill,
    colour = axis_gallery$y_axis$box_colour,
    linewidth = axis_gallery$y_axis$box_linewidth,
    outlier.shape = axis_gallery$y_axis$outlier_shape,
    outlier.size = axis_gallery$y_axis$outlier_size,
    outlier.alpha = axis_gallery$y_axis$outlier_alpha
  ) +
  scale_y_chromosome(
    human_karyotype,
    limits = rev(chromosome_order),
    axis_length = grid::unit(axis_gallery$y_axis$axis_length_mm, "mm"),
    body_width = grid::unit(axis_gallery$y_axis$body_width_mm, "mm"),
    label_gap = grid::unit(axis_gallery$y_axis$label_gap_mm, "mm"),
    fill = axis_gallery$body_fill,
    colour = axis_gallery$body_colour,
    linewidth = axis_gallery$body_linewidth
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = axis_gallery$y_axis$value_expand)
  ) +
  labs(
    x = "Genes / Mb",
    y = "Chromosome"
  ) +
  publication_theme +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = axis_gallery$axis_text_size)
  )

stopifnot(inherits(x_axis_plot$layers[[1]]$geom, "GeomBar"))
stopifnot(inherits(y_axis_plot$layers[[1]]$geom, "GeomBoxplot"))
stopifnot(inherits(ggplotGrob(x_axis_plot), "gtable"))
stopifnot(inherits(ggplotGrob(y_axis_plot), "gtable"))

save_axis_example(
  x_axis_plot,
  "original-data-chromosome-x-axis",
  axis_gallery$x_axis
)
save_axis_example(
  y_axis_plot,
  "original-data-chromosome-y-axis",
  axis_gallery$y_axis
)

writeLines(c(
  "Human chromosomes, GRCh38. RNA counts are 500 randomly selected GENCODE annotations, grouped by chromosome.",
  "Gene density uses all bundled windows, normalized to genes per Mb using each actual window width. Boxes show median and quartiles; whiskers extend to 1.5 IQR and points show outliers."
), file.path(axis_gallery$output_dir, "captions.txt"))
