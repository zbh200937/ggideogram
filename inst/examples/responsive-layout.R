# Run from the package checkout:
# Rscript -e 'pkgload::load_all("."); source("inst/examples/responsive-layout.R")'
library(ggideogram)
library(ggplot2)

# Real MCScanX Arabidopsis gene pairs; one middle-ranked pair per block.
pairs <- read.delim(system.file("extdata", "arabidopsis-synteny-pairs.tsv",
  package = "ggideogram"))
pairs <- pairs[pairs$Chr1 != pairs$Chr2, ]
pairs <- do.call(rbind, lapply(split(pairs, pairs$Block), function(d) d[ceiling(nrow(d) / 2), ]))
pairs$Pos1 <- (pairs$Start1 - 1 + pairs$End1) / 2
pairs$Pos2 <- (pairs$Start2 - 1 + pairs$End2) / 2
genes <- rbind(data.frame(Chr = pairs$Chr1, Pos = pairs$Pos1),
  data.frame(Chr = pairs$Chr2, Pos = pairs$Pos2))
karyotype <- data.frame(Chr = c("1", "5"), Start = 0, End = c(30427671, 26975502))
p <- ggideogram(karyotype, orientation = "horizontal", axis = TRUE,
  axis_units = "Mb", axis_n = 3) +
  geom_chrlink(data = pairs,
    aes(chr1 = Chr1, position1 = Pos1, chr2 = Chr2, position2 = Pos2),
    colour = "#9DA5AD", linewidth = 0.3) +
  geom_locus(data = genes, aes(chr = Chr, position = Pos),
    colour = "#4477AA", size = 2, gap = 0) +
  theme(plot.margin = margin(5, 5, 5, 10, unit = "mm"))

out <- "work/responsive-layout/examples"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_pair <- function(plot, name, width, height) {
  ggsave(file.path(out, paste0(name, ".png")), plot,
    width = width, height = height, units = "mm", dpi = 160, bg = "white")
  ggsave(file.path(out, paste0(name, ".pdf")), plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf)
}

# Reuse one built grob across two different device proportions.
g <- ggplotGrob(p)
save_pair(g, "01-wide", 200, 65)
save_pair(g, "02-tall", 110, 130)

# Reuse the same plot in unequal panel slots and as an inset.
panels <- patchwork::wrap_plots(p, p, widths = c(1, 2))
save_pair(panels, "03-unequal-panels", 230, 100)
inset <- p + patchwork::inset_element(p, left = 0.58, bottom = 0.30,
  right = 0.98, top = 0.80)
save_pair(inset, "04-inset", 200, 140)
