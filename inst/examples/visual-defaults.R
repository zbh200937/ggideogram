# Run from the package checkout:
# Rscript -e 'pkgload::load_all("."); source("inst/examples/visual-defaults.R")'
# Real Arabidopsis MCScanX block 15, bundled with ggideogram.
library(ggideogram)
library(ggplot2)

pairs <- read.delim(system.file("extdata", "arabidopsis-synteny-pairs.tsv",
  package = "ggideogram"))
pairs <- pairs[pairs$Block == 15, ]
genes <- unique(rbind(
  data.frame(Chr = pairs$Chr1, Pos = (pairs$Start1 - 1 + pairs$End1) / 2,
    Gene = pairs$Gene1),
  data.frame(Chr = pairs$Chr2, Pos = (pairs$Start2 - 1 + pairs$End2) / 2,
    Gene = pairs$Gene2)))
karyotype <- data.frame(Chr = "1", Start = 0, End = 30427671)
view <- chr_view(karyotype, "1", min(genes$Pos) - 20000, max(genes$Pos) + 20000)

p <- ggideogram(view, orientation = "horizontal", axis = TRUE,
  axis_units = "Mb", axis_colour = "black", name_colour = "black") +
  geom_locus(data = genes, aes(chr = Chr, position = Pos),
    colour = "#1B9E77", gap = 0, size = 2) +
  geom_locus(data = genes, aes(chr = Chr, position = Pos, label = Gene),
    geom = "text", position = "spread", fontface = "italic") +
  theme_classic(base_size = 11) +
  theme(plot.margin = margin(5, 5, 5, 15, unit = "mm"))

out <- "work/visual-defaults"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (width in c(140, 220)) {
  ggsave(file.path(out, paste0("arabidopsis-", width, ".png")), p,
    width = width, height = 70, units = "mm", dpi = 180, bg = "white")
  ggsave(file.path(out, paste0("arabidopsis-", width, ".pdf")), p,
    width = width, height = 70, units = "mm", device = grDevices::cairo_pdf)
}
