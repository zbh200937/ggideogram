# Six reproducible gallery figures built only from data bundled with
# ggideogram.  Run from a source checkout after loading the package:
#
#   Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/original-data-gallery.R")'

suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(patchwork)
})
source(system.file('examples', 'gallery-style.R', package = 'ggideogram'))

# These are named output and presentation choices, not renderer constants.
# Chromosome columns, axis chromosomes, track limits, inset locus and reverse
# inset corner are derived below from the bundled data and the declared slot.
gallery <- list(
  output_dir = file.path("work", "api-optimization", "examples", "original-data-gallery"),
  dpi = 200,
  base_size = 9,
  focus_chromosomes = 4L,
  original = list(width_mm = 185, height_mm = 135),
  tracks = list(width_mm = 185, height_mm = 135),
  insets = list(
    width_mm = 185,
    height_mm = 135,
    child_width_mm = 42,
    child_height_mm = 35
  ),
  marker = list(size = 0.70, stroke = 0.18, separation = 0.18),
  reverse_inset_candidates = data.frame(
    id = c("top_left", "top_right", "bottom_left", "bottom_right"),
    left = c(0.02, 0.60, 0.02, 0.60),
    bottom = c(0.56, 0.56, 0.08, 0.08),
    right = c(0.42, 1.00, 0.42, 1.00),
    top = c(1.00, 1.00, 0.50, 0.50)
  )
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

# A value of one half centres each heatmap tile in a unit-width overlay track;
# the mapped `fill` remains the observed density value.
gene_density$Overlay <- 0.5

shape_lookup <- c(circle = 16, box = 15, triangle = 17)
rna_key <- unique(Random_RNAs_500[c("Type", "Shape", "color")])
unknown_shapes <- setdiff(rna_key$Shape, names(shape_lookup))
if (length(unknown_shapes)) {
  stop("No standard ggplot2 shape mapping for: ",
       paste(unknown_shapes, collapse = ", "))
}
rna_shapes <- stats::setNames(shape_lookup[rna_key$Shape], rna_key$Type)
rna_colours <- stats::setNames(paste0("#", rna_key$color), rna_key$Type)

gene_colour <- "#277DA1"
ltr_colour <- "#43AA8B"
highlight_colour <- "#D55E00"
context_colour <- "#B8B8B8"
heatmap_colours <- c("#F2F5F8", "#9BBACD", "#315A7D")

common_theme <- theme_minimal(base_size = gallery$base_size) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot"
  )

row_axis_chromosomes <- function(karyotype, columns) {
  karyotype$Chr[seq.int(1L, nrow(karyotype), by = columns)]
}

save_gallery_plot <- function(plot, stem, dimensions) {
  dir.create(gallery$output_dir, recursive = TRUE, showWarnings = FALSE)
  png_file <- file.path(gallery$output_dir, paste0(stem, ".png"))
  pdf_file <- file.path(gallery$output_dir, paste0(stem, ".pdf"))
  ggsave(
    png_file, plot,
    width = dimensions$width_mm, height = dimensions$height_mm,
    units = "mm", dpi = gallery$dpi, bg = "white"
  )
  ggsave(
    pdf_file, plot,
    width = dimensions$width_mm, height = dimensions$height_mm,
    units = "mm", device = grDevices::cairo_pdf, bg = "white"
  )
  message("Wrote ", normalizePath(png_file))
  message("Wrote ", normalizePath(pdf_file))
  invisible(c(png = png_file, pdf = pdf_file))
}

add_rna_scales <- function(plot, legend = TRUE) {
  plot <- plot +
    scale_shape_manual(values = rna_shapes, name = "RNA type") +
    scale_colour_manual(values = rna_colours, name = "RNA type")
  if (!legend) plot <- plot + guides(shape = "none", colour = "none")
  plot
}

# -------------------------------------------------------------------------
# 1. Classic RIdeogram semantics rebuilt with standard ggplot2 components.

classic_tracks <- track_layout(density = geom_track(side = "overlay", width = 0.9, limits = c(0, 1)), markers = geom_track(side = "right",
    width = 0.55, gap = 0.18))
classic_grid <- gallery_grid(
  human_karyotype, classic_tracks,
  gallery$original$width_mm - 18, gallery$original$height_mm - 32,
  orientations = "vertical", row_gap = 5
)
classic_columns <- classic_grid$ncol
classic_axes <- row_axis_chromosomes(human_karyotype, classic_columns)
classic_repel <- position_chr_repel(
  gallery$marker$separation, units = "body_width"
)

classic_plot <- ggideogram(
  human_karyotype,
  ncol = classic_columns,
  row_gap = 5,
  orientation = classic_grid$orientation,
  tracks = classic_tracks,
  axis = classic_axes,
  axis_breaks = seq(0, max(human_karyotype$End), by = 50e6),
  axis_side = "left",
  base_family = "sans"
) +
  geom_track(data = gene_density, mapping = aes(chr = Chr, x = Position, y = Overlay,
      width = Width, fill = Rate), track = "density", geom = ggplot2::geom_tile(height = 1, colour = NA)) +
  # The outline is redrawn as a public component after the overlay tiles.
  geom_chr(component = "body", fill = NA, cytoband = FALSE, linewidth = 0.36) +
  geom_locus(geom = "link", data = Random_RNAs_500, mapping = aes(chr = Chr, position = Position),
      position = classic_repel, track = "markers", colour = "#777777", linewidth = 0.12,
      alpha = 0.45) +
  geom_locus(geom = "point", data = Random_RNAs_500, mapping = aes(chr = Chr, position = Position,
      shape = Type, colour = Type), position = classic_repel, track = "markers", size = gallery$marker$size,
      stroke = gallery$marker$stroke) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    limits = range(gene_density$Rate, na.rm = TRUE),
    name = "Genes / Mb"
  ) +
  labs(
    title = "Human chromosomes · original RIdeogram data",
    subtitle = "Gene density and 500 randomly sampled RNA annotations",
    caption = "Density uses actual window width; RNA markers are a sample of GENCODE annotations."
  ) +
  guides(
    fill = guide_colourbar(title.position = "top"),
    shape = guide_legend(title.position = "top", nrow = 1, override.aes = list(size = 1.8)),
    colour = guide_legend(title.position = "top", nrow = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = gallery_legend_spacing,
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(size = gallery$base_size),
    legend.text = element_text(size = gallery$base_size),
    plot.title = element_text(face = "bold", size = gallery$base_size + 1),
    plot.subtitle = element_text(size = gallery$base_size),
    plot.caption = element_text(size = 8, hjust = 0)
  )
classic_plot <- add_rna_scales(classic_plot)

save_gallery_plot(
  classic_plot, "gallery-classic-rebuilt", gallery$original
)

# -------------------------------------------------------------------------
# 2. Ordinary geoms connected inward through bp-aligned tracks.

focus_count <- min(gallery$focus_chromosomes, nrow(human_karyotype))
focus_chr <- human_karyotype$Chr[seq_len(focus_count)]
focus_karyotype <- human_karyotype[human_karyotype$Chr %in% focus_chr, ]
focus_gene <- gene_density[gene_density$Chr %in% focus_chr, ]
focus_ltr <- LTR_density[LTR_density$Chr %in% focus_chr, ]
focus_rna <- Random_RNAs_500[Random_RNAs_500$Chr %in% focus_chr, ]

inward_tracks <- track_layout(density = geom_track(side = "overlay", width = 0.78, limits = c(0, 1)), markers = geom_track(side = "right",
    width = 0.55, gap = 0.18), genes = geom_track(side = "right", width = 1.8, gap = 0.24,
    limits = c(0, max(focus_gene$Value, na.rm = TRUE))), ltr = geom_track(side = "right",
    width = 1.8, gap = 0.24, limits = c(0, max(focus_ltr$Value, na.rm = TRUE))))
inward_grid <- gallery_grid(
  focus_karyotype, inward_tracks,
  gallery$tracks$width_mm - 18, gallery$tracks$height_mm - 38,
  orientations = c("vertical", "horizontal"), padding = 1
)
inward_columns <- inward_grid$ncol
inward_axes <- row_axis_chromosomes(focus_karyotype, inward_columns)
inward_repel <- position_chr_repel(
  gallery$marker$separation, units = "body_width"
)

inward_plot <- ggideogram(
  focus_karyotype,
  ncol = inward_columns,
  orientation = inward_grid$orientation,
  tracks = inward_tracks,
  axis = inward_axes,
  axis_side = "left",
  padding = 1,
  base_family = "sans"
) +
  geom_track(data = focus_gene, mapping = aes(chr = Chr, x = Position, y = Overlay,
      width = Width, fill = Rate), track = "density", geom = ggplot2::geom_tile(height = 1, colour = NA)) +
  geom_chr(component = "body", fill = NA, cytoband = FALSE, linewidth = 0.38) +
  geom_track(data = focus_gene, mapping = aes(chr = Chr, x = Position, y = Value, group = Chr),
      track = "genes", geom = ggplot2::geom_line(colour = gene_colour, linewidth = 0.32)) +
  geom_track(data = focus_ltr, mapping = aes(chr = Chr, x = Position, y = Value, width = End -
      Start + 1), track = "ltr", geom = ggplot2::geom_col(position = "identity", fill = ltr_colour,
      alpha = 0.72)) +
  geom_locus(geom = "link", data = focus_rna, mapping = aes(chr = Chr, position = Position),
      position = inward_repel, track = "markers", colour = "#777777", linewidth = 0.14,
      alpha = 0.45) +
  geom_locus(geom = "point", data = focus_rna, mapping = aes(chr = Chr, position = Position,
      shape = Type, colour = Type), position = inward_repel, track = "markers", size = gallery$marker$size,
      stroke = gallery$marker$stroke) +
  geom_track(track = "genes", axis = gallery_axis(chr = focus_chr[1], position = "start", breaks = c(0,
      max(focus_gene$Value, na.rm = TRUE)), labels = scales::label_number(accuracy = 1))) +
  geom_track(track = "ltr", axis = gallery_axis(chr = focus_chr[min(2L, length(focus_chr))], position = "start",
      breaks = c(0, max(focus_ltr$Value, na.rm = TRUE)), labels = scales::label_number(accuracy = 1))) +
  scale_fill_gradientn(
    colours = heatmap_colours,
    limits = range(focus_gene$Rate, na.rm = TRUE),
    name = "Genes / Mb"
  ) +
  labs(
    title = "Ordinary ggplot2 geoms as chromosome components",
    subtitle = paste0(
      "Counts/window: blue = genes, green = LTRs (chr ",
      paste(focus_chr, collapse = ", "), ")"
    ),
    caption = "Nominal 1 Mb bins; terminal bins are shorter. RNA markers are randomly sampled annotations."
  ) +
  guides(
    fill = guide_colourbar(title.position = "top"),
    shape = guide_legend(title.position = "top", nrow = 1, override.aes = list(size = 1.8)),
    colour = guide_legend(title.position = "top", nrow = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = gallery_legend_spacing,
    legend.box = "horizontal",
    legend.box.just = "left",
    plot.title = element_text(face = "bold", size = gallery$base_size + 1),
    plot.caption = element_text(size = 8, hjust = 0),
    plot.subtitle = element_text(
      size = gallery$base_size,
      margin = margin(b = 8)
    )
  )
inward_plot <- add_rna_scales(inward_plot)

save_gallery_plot(
  inward_plot, "gallery-inward-tracks", gallery$tracks
)

# -------------------------------------------------------------------------
# 3. Complete plots connected in both directions.

windows <- merge(
  gene_density[c("Chr", "Start", "End", "Position", "Value")],
  LTR_density[c("Chr", "Start", "End", "Value")],
  by = c("Chr", "Start", "End"),
  suffixes = c("_gene", "_ltr"),
  sort = FALSE
)
names(windows)[names(windows) == "Position"] <- "Position_gene"
selected_window <- windows[which.max(windows$Value_gene), , drop = FALSE]
selected_chr <- as.character(selected_window$Chr)
selected_karyotype <- human_karyotype[
  human_karyotype$Chr == selected_chr, , drop = FALSE
]
selected_gene <- gene_density[gene_density$Chr == selected_chr, ]
selected_rna <- Random_RNAs_500[Random_RNAs_500$Chr == selected_chr, ]

locus_values <- data.frame(
  Metric = factor(c("Genes", "LTRs"), levels = c("Genes", "LTRs")),
  Count = c(selected_window$Value_gene, selected_window$Value_ltr)
)
locus_child <- ggplot(locus_values, aes(Metric, Count, fill = Metric)) +
  geom_col(width = 0.68, show.legend = FALSE) +
  scale_fill_manual(values = c(Genes = gene_colour, LTRs = ltr_colour)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Selected window",
    subtitle = paste0("chr", selected_chr, "\n",
                      format(selected_window$Start, big.mark = ","), "–",
                      format(selected_window$End, big.mark = ",")),
    x = NULL,
    y = "Count"
  ) +
  theme_minimal(base_size = gallery$base_size - 1) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.margin = margin(4, 5, 4, 5)
  )

locus_inset <- selected_window[c("Chr", "Start", "End")]
locus_inset$Plot <- I(list(locus_child))
locus_tracks <- track_layout(genes = geom_track(side = "right", width = 1.6, gap = 0.24, limits = c(0,
    max(selected_gene$Value, na.rm = TRUE))), locus = geom_track(side = "right", width = 4.5,
    gap = 0.4))
locus_midpoint <- (selected_window$Start - 1 + selected_window$End) / 2
chromosome_midpoint <-
  (selected_karyotype$Start + selected_karyotype$End) / 2
locus_vjust <- if (locus_midpoint > chromosome_midpoint) 0 else 1

host_plot <- ggideogram(
  selected_karyotype,
  ncol = 1,
  tracks = locus_tracks,
  axis = selected_chr,
  axis_side = "left",
  base_family = "sans"
) +
  geom_track(data = selected_gene, mapping = aes(chr = Chr, x = Position, y = Value, group = Chr),
      track = "genes", geom = ggplot2::geom_line(colour = gene_colour, linewidth = 0.38)) +
  geom_locus(geom = "interval", data = selected_window, mapping = aes(chr = Chr, start = Start,
      end = End), colour = highlight_colour, linewidth = 1.2) +
  geom_track(track = "genes", axis = gallery_axis(chr = selected_chr, position = "start", breaks = c(0,
      max(selected_gene$Value, na.rm = TRUE)), labels = scales::label_number(accuracy = 1))) +
  geom_locus_inset(data = locus_inset, mapping = aes(chr = Chr, start = Start, end = End,
      plot = Plot), track = "locus", width = grid::unit(gallery$insets$child_width_mm, "mm"),
      height = grid::unit(gallery$insets$child_height_mm, "mm"), vjust = locus_vjust) +
  labs(
    tag = "A",
    title = "Complete ggplot at a locus",
    subtitle = "The child keeps its own axes and theme"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = gallery$base_size + 1),
    plot.subtitle = element_text(size = gallery$base_size),
    plot.margin = margin(6, 6, 6, 6)
  )

mini_tracks <- track_layout(markers = geom_track(side = "right", width = 0.5, gap = 0.16))
mini_plot <- ggideogram(
  selected_karyotype,
  ncol = 1,
  tracks = mini_tracks,
  axis = selected_chr,
  axis_side = "left",
  axis_size = 1.8,
  name_size = 2.4,
  padding = 0.25,
  base_family = "sans"
) +
  geom_locus(geom = "point", data = selected_rna, mapping = aes(chr = Chr, position = Position,
      shape = Type, colour = Type), track = "markers", size = 0.55, stroke = 0.15) +
  theme(
    legend.position = "none",
    plot.background = element_rect(
      fill = scales::alpha("white", 0.94),
      colour = "#C8C8C8",
      linewidth = 0.3
    ),
    plot.margin = margin(3, 3, 3, 3)
  )
mini_plot <- add_rna_scales(mini_plot, legend = FALSE)

windows$Selected <- windows$Chr == selected_chr
scatter_plot <- ggplot(windows, aes(Value_gene, Value_ltr)) +
  geom_point(
    data = windows[!windows$Selected, ],
    colour = context_colour,
    size = 0.55,
    alpha = 0.38
  ) +
  geom_point(
    data = windows[windows$Selected, ],
    colour = gene_colour,
    size = 0.72,
    alpha = 0.72
  ) +
  labs(
    tag = "B",
    title = "Ideogram inside a ggplot",
    subtitle = paste0("Blue points = chromosome ", selected_chr),
    x = "Genes / window",
    y = "LTRs / window",
    caption = "Nominal 1 Mb bins; terminal bins are shorter."
  ) +
  common_theme +
  theme(plot.margin = margin(6, 6, 6, 6))

normalise <- function(x) {
  span <- diff(range(x, na.rm = TRUE))
  if (span == 0) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / span
}
nx <- normalise(windows$Value_gene)
ny <- normalise(windows$Value_ltr)
candidate_counts <- vapply(
  seq_len(nrow(gallery$reverse_inset_candidates)),
  function(index) {
    box <- gallery$reverse_inset_candidates[index, ]
    sum(nx >= box$left & nx <= box$right &
        ny >= box$bottom & ny <= box$top, na.rm = TRUE)
  },
  integer(1)
)
inset_box <- gallery$reverse_inset_candidates[
  which.min(candidate_counts), , drop = FALSE
]

reverse_plot <- scatter_plot + patchwork::inset_element(
  mini_plot,
  left = inset_box$left,
  bottom = inset_box$bottom,
  right = inset_box$right,
  top = inset_box$top,
  align_to = "panel",
  clip = TRUE
)

inset_composition <- host_plot | reverse_plot
inset_composition <- inset_composition +
  plot_layout(widths = c(1, 1.12))

save_gallery_plot(
  inset_composition, "gallery-bidirectional-insets", gallery$insets
)

# -------------------------------------------------------------------------
# 4. The existing outward patchwork example remains the full-genome gallery
#    panel. It is sourced in a clean child environment.

source(
  file.path("inst", "examples", "original-data-composition.R"),
  local = new.env(parent = globalenv())
)

# -------------------------------------------------------------------------
# 5-6. Native x/y chromosome axes. These are separate figures because their
#      defining property is that each is one ordinary ggplot with one native
#      scale/guide, not a multi-panel composition.

source(
  file.path("inst", "examples", "original-data-axis-integration.R"),
  local = new.env(parent = globalenv())
)

# Stable outward boundaries used by other graphics systems.
stopifnot(inherits(as_ideogram_grob(classic_plot), "gtable"))
stopifnot(inherits(ggplotGrob(inward_plot), "gtable"))
