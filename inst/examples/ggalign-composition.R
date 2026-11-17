# Run after loading ggideogram from an installation or source checkout.
suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(ggalign)
})

data("human_karyotype", package = "ggideogram")
data("Random_RNAs_500", package = "ggideogram")

karyotype <- human_karyotype[seq_len(6), ]
rna <- subset(Random_RNAs_500, Chr %in% karyotype$Chr)
counts <- table(factor(as.character(rna$Chr), levels = as.character(karyotype$Chr)))
rna_matrix <- matrix(as.numeric(counts), ncol = 1,
  dimnames = list(names(counts), "RNAs"))
summary <- data.frame(Chr = names(counts), Count = as.numeric(counts))

p <- ggideogram(karyotype, ncol = nrow(karyotype)) +
  geom_locus(data = rna,
    aes(chr = Chr, position = (Start - 1 + End) / 2, colour = Type), size = 1.2)
bars <- ggplot(summary, aes(Chr, Count)) +
  geom_col(fill = "#4477AA") + labs(y = "Bundled RNA records") + theme_minimal()

# Complete panels retain their own coordinate systems.
panels <- ggalign::align_plots(p, bars, ncol = 2, widths = c(2, 1))

# Row names carry chromosome identity through ggalign ordering.
aligned <- ggalign::stack_discretev(rna_matrix) +
  ggalign::align_order() + ggalign::ggalign() +
  geom_col(aes(x = .discrete_x, y = value), fill = "#4477AA") +
  scale_x_chromosome(karyotype, limits = NULL) +
  labs(y = "Bundled RNA records") + theme(axis.text.x = element_text())

# The composition package builds a standard grob for the anchored viewport.
child <- ggalign::align_plots(bars, bars + coord_flip(), ncol = 2)
embedded <- ggideogram(karyotype[1, ],
  tracks = list(summary = geom_track(width = 35))) +
  geom_locus_inset(ggalign::ggalignGrob(child),
    data = data.frame(Chr = karyotype$Chr[1], Pos = karyotype$End[1] / 2),
    mapping = aes(chr = Chr, position = Pos), track = "summary",
    width = grid::unit(.8, "npc"), height = grid::unit(.7, "npc"))

output_dir <- file.path("work", "ggalign-examples")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plots <- list(panels = panels, axis = aligned, inset = embedded)
for (name in names(plots)) {
  ggplot2::ggsave(file.path(output_dir, paste0("ggalign-", name, ".png")),
    plots[[name]], width = 210, height = 115, units = "mm", dpi = 150, bg = "white")
}
