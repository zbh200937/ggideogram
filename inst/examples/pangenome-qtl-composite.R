# Publication-style multi-scale QTL fine-mapping and 4-genome pangenome synteny
# composite figures built with ggideogram, ggplot2, ggnewscale and patchwork.
#
# From a source checkout:
#   Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/pangenome-qtl-composite.R")'

suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(patchwork)
})

out_dir <- Sys.getenv(
  "GGIDEOGRAM_EXAMPLE_OUTDIR",
  unset = file.path("work", "examples", "pangenome-qtl-composite")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Multi-scale GWAS/QTL fine-mapping 3-panel composite
#    Panel A: Chromosome 1 overview + candidate QTL window highlight
#    Panel B: chr_view zoomed Manhattan + geom_hline threshold + geom_vline
#             lead variant + Araport11 gene models + ggnewscale
#    Panel C: Haplotype phenotypic effect boxplot + jitter
# ---------------------------------------------------------------------------
gff_path <- system.file("extdata", "arabidopsis-first-genes.gff3", package = "ggideogram")
features <- read_chr_features(gff_path)
gene_models <- resolve_gene_models(features)

ath_chr1 <- data.frame(
  Chr = "1", Start = 0, End = 30427671,
  CE_start = 14500000, CE_end = 15500000
)
ath_semantic <- as_ideogram_data(ath_chr1)

set.seed(42)
bg_bins <- seq(250000, 30250000, by = 500000)
overview_density <- data.frame(
  Chr = "1",
  Pos = bg_bins,
  Width = 500000,
  Density = round(45 + 28 * cos((bg_bins - 15e6) / 6e6)^2 + stats::rnorm(length(bg_bins), 0, 4), 1)
)
qtl_region <- data.frame(Chr = "1", Left = 1, Right = 600000, Lo = 0, Hi = 1)

panel_qtl_a <- ggideogram(
  ath_semantic,
  orientation = "horizontal",
  axis = TRUE,
  axis_units = "Mb",
  chromosome_width = 1.4,
  tracks = list(
    heat = geom_track(
      side = "overlay", width = 1,
      data = overview_density,
      mapping = aes(chr = Chr, x = Pos, y = 0.5, fill = Density),
      geom = geom_tile(width = 500000, height = 1)
    ),
    window = geom_track(
      side = "left", width = 1.0, gap = 0.25, limits = c(0, 1),
      data = qtl_region,
      mapping = aes(chr = Chr, xmin = Left, xmax = Right, ymin = Lo, ymax = Hi),
      geom = geom_rect(fill = "#D55E00", colour = "#9E3D00", linewidth = 0.35)
    )
  )
) +
  scale_fill_gradient(low = "#EEF4F8", high = "#1F5A8A", name = "Genes / 0.5 Mb") +
  labs(title = "A  Chromosome 1 overview & candidate QTL region (0-15 kb)") +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "right",
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(2.6, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  )

zoom_view <- chr_view(ath_semantic, chr = "1", start = 2500, end = 14800)
snp_pos <- seq(2800, 14500, by = 190)
dist_kb <- abs(snp_pos - 7950) / 1000
lod <- pmax(0.4, 9.2 * exp(-(dist_kb / 1.45)^2) + stats::runif(length(snp_pos), 0.1, 1.3))
ld_r2 <- pmin(1, pmax(0.05, exp(-dist_kb / 1.8) + stats::rnorm(length(snp_pos), 0, 0.04)))
gwas_df <- data.frame(
  Chr = "1",
  Pos = snp_pos,
  LOD = lod,
  Lo = pmax(0, lod - 0.45),
  Hi = lod + 0.45,
  LD = cut(ld_r2, breaks = c(0, 0.2, 0.5, 0.8, 1.01),
           labels = c("<0.2", "0.2-0.5", "0.5-0.8", ">0.8"), include.lowest = TRUE)
)
lead_pos <- gwas_df$Pos[which.max(gwas_df$LOD)]

candidate_calls <- data.frame(
  Chr = "1",
  Pos = c(4180, lead_pos, 12550),
  Label = c("Chr1:4180", sprintf("Lead C>T (%s bp)", lead_pos), "Chr1:12550")
)

panel_qtl_b <- ggideogram(
  zoom_view,
  orientation = "horizontal",
  axis = TRUE,
  axis_side = "right",
  axis_units = "kb",
  chromosome_width = 0.9,
  tracks = list(
    gwas = geom_track(
      side = "left", width = 6.5, gap = 0.35,
      limits = c(0, 10.5),
      label = list(text = "-log10(P)", size = 2.6),
      axis = list(breaks = c(0, 5, 10), size = 2.3),
      data = gwas_df,
      mapping = aes(chr = Chr, x = Pos, y = LOD),
      layers = list(
        geom_hline(yintercept = 5.3, linetype = "dashed", colour = "#B22222", linewidth = 0.35),
        geom_vline(xintercept = lead_pos, linetype = "dotted", colour = "#555555", linewidth = 0.3),
        geom_linerange(aes(ymin = Lo, ymax = Hi), colour = "#B8C2CC", linewidth = 0.25),
        geom_point(aes(colour = LD), size = 1.6)
      )
    )
  )
) +
  scale_colour_manual(
    values = c("<0.2" = "#3B528B", "0.2-0.5" = "#21918C", "0.5-0.8" = "#FDE725", ">0.8" = "#D55E00"),
    name = expression(italic(r)^2)
  )

if (requireNamespace("ggnewscale", quietly = TRUE)) {
  panel_qtl_b <- panel_qtl_b + ggnewscale::new_scale_fill()
}

panel_qtl_b <- panel_qtl_b +
  geom_track(
    track = "genes", side = "right", width = 5.5, gap = 0.35,
    label = list(text = "Araport11", size = 2.6),
    data = gene_models,
    mapping = aes(chr = Chr, start = Start, end = End, type = Type,
                  strand = Strand, gene = Symbol, fill = Strand),
    layers = list(
      geom_genemodel(mode = "gene", label_size = 2.6, label_fontface = "italic")
    )
  ) +
  scale_fill_manual(values = c("+" = "#2B6CB0", "-" = "#C53030"), name = "Strand") +
  geom_locus(
    data = candidate_calls,
    mapping = aes(chr = Chr, position = Pos, label = Label),
    geom = "text", position = "spread", side = "right",
    label_width = 3.2, connection_height = 2.5, size = 2.6
  ) +
  labs(title = "B  Fine-mapping association peak & Araport11 gene structure (2.5-14.8 kb)") +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.key.size = grid::unit(3, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  )

hap_df <- rbind(
  data.frame(Haplotype = "Hap1\n(Ref C/C)", Phenotype = stats::rnorm(28, 42.5, 3.1)),
  data.frame(Haplotype = "Hap2\n(Alt T/T)", Phenotype = stats::rnorm(24, 54.8, 3.4)),
  data.frame(Haplotype = "Hap3\n(Het C/T)", Phenotype = stats::rnorm(18, 48.2, 3.2))
)
hap_df$Haplotype <- factor(hap_df$Haplotype, levels = unique(hap_df$Haplotype))

panel_qtl_c <- ggplot(hap_df, aes(x = Haplotype, y = Phenotype, fill = Haplotype)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.82, linewidth = 0.35) +
  geom_jitter(width = 0.14, size = 1.2, alpha = 0.7, colour = "#202020") +
  scale_fill_manual(values = c("#3B528B", "#D55E00", "#21918C"), guide = "none") +
  labs(
    title = "C  Allelic effect at ARV1",
    x = NULL,
    y = "Root length (mm)"
  ) +
  theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    axis.text = element_text(colour = "#202020")
  )

qtl_composite <- ((panel_qtl_a / panel_qtl_b) + plot_layout(heights = c(1, 2.35))) |
  panel_qtl_c
qtl_composite <- qtl_composite + plot_layout(widths = c(2.65, 1))

ggsave(
  file.path(out_dir, "qtl-finemapping-composite.png"),
  qtl_composite, width = 235, height = 145, units = "mm", dpi = 220, bg = "white"
)

# ---------------------------------------------------------------------------
# 2. Four-genome stacked pangenome microsynteny + native chromosome-axis summary
# ---------------------------------------------------------------------------
pan_kar <- data.frame(
  Genome = c("Nipponbare", "ZS97", "MH63", "O_rufipogon"),
  Assembly = c("IRGSP-1.0", "ZS97RS3", "MH63RS3", "W1943"),
  Chr = c("Chr03", "Chr03", "Chr03", "Chr03"),
  Start = c(16000000, 16200000, 15900000, 15800000),
  End = c(17200000, 17350000, 17150000, 17050000),
  Label = c("Nipponbare", "ZS97", "MH63", "O. rufipogon")
)
pan_sem <- as_ideogram_data(
  pan_kar,
  mapping = aes(chr = Chr, start = Start, end = End,
                genome = Genome, assembly = Assembly, label = Label)
)
pan_keys <- pan_sem$karyotype$.chr

pan_blocks <- data.frame(
  chr1 = c(pan_keys[1], pan_keys[1], pan_keys[1],
           pan_keys[2], pan_keys[2], pan_keys[2],
           pan_keys[3], pan_keys[3], pan_keys[3]),
  start1 = c(16050000, 16480000, 16850000,
             16260000, 16660000, 17020000,
             15960000, 16390000, 16780000),
  end1   = c(16380000, 16760000, 17140000,
             16580000, 16940000, 17300000,
             16290000, 16680000, 17080000),
  chr2 = c(pan_keys[2], pan_keys[2], pan_keys[2],
           pan_keys[3], pan_keys[3], pan_keys[3],
           pan_keys[4], pan_keys[4], pan_keys[4]),
  start2 = c(16260000, 16660000, 17020000,
             15960000, 16390000, 16780000,
             15860000, 16300000, 16700000),
  end2   = c(16580000, 16940000, 17300000,
             16290000, 16680000, 17080000,
             16180000, 16600000, 16990000),
  Type = c("Collinear", "Inversion", "Collinear",
           "Collinear", "Inversion", "Collinear",
           "Collinear", "Collinear", "Collinear"),
  Strand = c("+", "-", "+",
             "+", "-", "+",
             "+", "+", "+")
)

cov_pos <- seq(16050000, 17150000, length.out = 24)
cov_nip <- data.frame(
  Key = pan_keys[1],
  Pos = cov_pos,
  Depth = 8 + 75 * exp(-((cov_pos - 16620000) / 180000)^2)
)
anchor_genes <- data.frame(
  Key = pan_keys,
  Pos = c(16620000, 16800000, 16535000, 16450000),
  Gene = c("GS3", "GS3-ZS", "GS3-MH", "GS3-Or")
)

panel_pan_a <- ggideogram(
  pan_sem,
  orientation = "horizontal",
  ncol = 4,
  axis = TRUE,
  axis_units = "kb",
  chromosome_gap = 7,
  tracks = list(
    depth = geom_track(
      chr = pan_keys[1], side = "left", width = 2.8, gap = 0.3,
      transform = "log10",
      label = list(text = "ChIP (log10)", size = 2.5),
      axis = list(size = 2.2),
      data = cov_nip,
      mapping = aes(chr = Key, x = Pos, y = Depth),
      geom = geom_area(fill = "#2B6CB0", alpha = 0.45, colour = "#1E4E8C", linewidth = 0.3)
    )
  )
) +
  geom_chrlink(
    data = pan_blocks,
    mapping = aes(chr1 = chr1, start1 = start1, end1 = end1,
                  chr2 = chr2, start2 = start2, end2 = end2,
                  orientation = Strand, fill = Type),
    alpha = 0.48
  ) +
  geom_locus(
    data = anchor_genes,
    mapping = aes(chr = Key, position = Pos, label = Gene),
    geom = "text", position = "spread", size = 2.6
  ) +
  scale_fill_manual(
    values = c(Collinear = "#7FB3D5", Inversion = "#E67E22"),
    name = "Synteny class"
  ) +
  labs(title = "A  Four-genome Oryza pangenome microsynteny & inversion at GS3") +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  )

sv_summary <- rbind(
  data.frame(Key = pan_keys, Class = "Core genes", Count = c(118, 115, 117, 112)),
  data.frame(Key = pan_keys, Class = "PAV / SV loci", Count = c(14, 22, 19, 27))
)

panel_pan_b <- ggplot(sv_summary, aes(x = Count, y = Key, fill = Class)) +
  geom_col(width = 0.58) +
  scale_y_chromosome(
    pan_sem,
    limits = rev(pan_keys),
    axis_length = grid::unit(12, "mm"),
    body_width = grid::unit(1.6, "mm")
  ) +
  scale_fill_manual(values = c("Core genes" = "#3C6E71", "PAV / SV loci" = "#D95F02"), name = NULL) +
  labs(title = "B  Locus repertoire", x = "Annotated features", y = NULL) +
  theme_classic(base_size = 8.5) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.key.size = grid::unit(3, "mm")
  )

pangenome_composite <- (panel_pan_a | panel_pan_b) + plot_layout(widths = c(2.5, 1))

ggsave(
  file.path(out_dir, "pangenome-synteny-composite.png"),
  pangenome_composite, width = 235, height = 125, units = "mm", dpi = 220, bg = "white"
)
