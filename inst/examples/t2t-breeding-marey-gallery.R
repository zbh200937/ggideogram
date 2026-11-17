# T2T phased haplotype rearrangements (SyRI + bedGraph) and plant breeding
# CSSL/RIL introgression BinMap + genetic-physical (cM <-> Mb) Marey composite.
#
# From a source checkout:
#   Rscript -e 'pkgload::load_all(".", export_all = FALSE); source("inst/examples/t2t-breeding-marey-gallery.R")'

suppressPackageStartupMessages({
  library(ggideogram)
  library(ggplot2)
  library(patchwork)
})

out_dir <- Sys.getenv(
  "GGIDEOGRAM_EXAMPLE_OUTDIR",
  unset = file.path("work", "examples", "t2t-breeding-marey-gallery")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. T2T Phased Haplotype Landscape (HapA vs HapB) with SyRI & bedGraph
# ---------------------------------------------------------------------------
t2t_kar <- data.frame(
  Genome = c("HapA", "HapB"),
  Chr = c("Chr03", "Chr03"),
  Start = c(0, 0),
  End = c(36.4e6, 35.8e6),
  CE_start = c(16.2e6, 15.8e6),
  CE_end = c(19.0e6, 18.5e6),
  Label = c("HapA · Chr03 (T2T)", "HapB · Chr03 (T2T)")
)
t2t_data <- as_ideogram_data(
  t2t_kar,
  mapping = aes(genome = Genome, chr = Chr, start = Start, end = End, label = Label),
  centromere = aes(start = CE_start, end = CE_end)
)

syri_tmp <- tempfile(fileext = ".syri.out")
bg_tmp <- tempfile(fileext = ".bedgraph")
on.exit(unlink(c(syri_tmp, bg_tmp)), add = TRUE)

writeLines(c(
  "Chr03\t1\t6800000\t-\t-\tChr03\t1\t6750000\tSYN1\t-\tSYN\t-",
  "Chr03\t7100000\t11400000\t-\t-\tChr03\t11200000\t6950000\tINV1\t-\tINV\t-",
  "Chr03\t11700000\t15900000\t-\t-\tChr03\t11500000\t15600000\tSYN2\t-\tSYN\t-",
  "Chr03\t19300000\t24800000\t-\t-\tChr03\t18800000\t24200000\tSYN3\t-\tSYN\t-",
  "Chr03\t25200000\t27900000\t-\t-\tChr03\t29800000\t32400000\tTRANS1\t-\tTRANS\t-",
  "Chr03\t28200000\t31000000\t-\t-\tChr03\t24600000\t27300000\tDUP1\t-\tDUP\tcopygain",
  "Chr03\t31300000\t36400000\t-\t-\tChr03\t32600000\t35800000\tSYN4\t-\tSYN\t-"
), syri_tmp)

syri_links <- read_chr_links(
  syri_tmp,
  format = "syri",
  target_genome = "HapA",
  query_genome = "HapB",
  karyotype = t2t_data
)
syri_links$SV_Class <- factor(
  syri_links$Type,
  levels = c("SYN", "INV", "TRANS", "DUP")
)

bins_a <- seq(0, 36.0e6, by = 6e5)
mid_a <- bins_a + 3e5
cg_a <- pmin(98, pmax(18, 28 + 64 * exp(-((mid_a - 17.6e6) / 3.4e6)^2) +
  18 * exp(-((mid_a - 9.2e6) / 1.8e6)^2)))
bins_b <- seq(0, 35.2e6, by = 6e5)
mid_b <- bins_b + 3e5
cg_b <- pmin(98, pmax(18, 27 + 65 * exp(-((mid_b - 17.1e6) / 3.3e6)^2) +
  16 * exp(-((mid_b - 9.0e6) / 1.8e6)^2)))

bg_lines <- c(
  "track type=bedGraph name=\"CpG_methylation\"",
  sprintf("HapA_Chr03\t%.0f\t%.0f\t%.2f", bins_a, pmin(bins_a + 6e5, 36.4e6), cg_a),
  sprintf("HapB_Chr03\t%.0f\t%.0f\t%.2f", bins_b, pmin(bins_b + 6e5, 35.8e6), cg_b)
)
writeLines(bg_lines, bg_tmp)
cpg_bg <- read_chr_features(
  bg_tmp,
  format = "bedgraph",
  chr_map = stats::setNames(
    c(chr_key("HapA", "Chr03"), chr_key("HapB", "Chr03")),
    c("HapA_Chr03", "HapB_Chr03")
  ),
  karyotype = t2t_data
)

landmarks <- data.frame(
  Key = c(
    chr_key("HapA", "Chr03"), chr_key("HapA", "Chr03"), chr_key("HapA", "Chr03"),
    chr_key("HapB", "Chr03"), chr_key("HapB", "Chr03")
  ),
  Start = c(15.9e6, 3.1e6, 28.2e6, 15.5e6, 3.0e6),
  End = c(19.3e6, 3.8e6, 28.7e6, 18.8e6, 3.6e6),
  Feature = c("CENH3 / satellite", "45S rDNA", "5S rDNA", "CENH3 / satellite", "45S rDNA")
)

telomeres <- data.frame(
  Key = rep(c(chr_key("HapA", "Chr03"), chr_key("HapB", "Chr03")), each = 2),
  Pos = c(0, 36.4e6, 0, 35.8e6),
  Side = c("left", "left", "right", "right")
)

p_t2t <- ggideogram(
  t2t_data,
  orientation = "horizontal",
  chromosome_width = 1.25,
  chromosome_gap = 10,
  fill = "#F3F5F8",
  colour = "#495057",
  name_size = 3.1,
  tracks = list(
    body_feat = geom_track(
      side = "overlay", width = 0.92,
      data = landmarks,
      mapping = aes(chr = Key, start = Start, end = End, fill = Feature),
      layers = list(geom_chr(component = "fill"))
    ),
    cpg_a = geom_track(
      side = "left", width = 2.6, gap = 0.45,
      chr = chr_key("HapA", "Chr03"),
      limits = c(0, 100),
      data = cpg_bg,
      mapping = aes(chr = Chr, x = Mid, y = Value),
      label = list(text = "5mCG (%)", size = 2.5),
      axis = list(breaks = c(0, 50, 100), size = 2.1),
      layers = list(
        geom_area(fill = "#9ECAE1", alpha = 0.55, colour = NA),
        geom_smooth(method = "loess", formula = y ~ x, span = 0.25, se = FALSE,
          colour = "#08519C", linewidth = 0.55)
      )
    ),
    cpg_b = geom_track(
      side = "right", width = 2.6, gap = 0.45,
      chr = chr_key("HapB", "Chr03"),
      limits = c(0, 100),
      data = cpg_bg,
      mapping = aes(chr = Chr, x = Mid, y = Value),
      label = list(text = "5mCG (%)", size = 2.5),
      axis = list(breaks = c(0, 50, 100), size = 2.1),
      layers = list(
        geom_area(fill = "#C7E9C0", alpha = 0.55, colour = NA),
        geom_smooth(method = "loess", formula = y ~ x, span = 0.25, se = FALSE,
          colour = "#006D2C", linewidth = 0.55)
      )
    )
  )
) +
  scale_fill_manual(
    values = c(
      "CENH3 / satellite" = "#B2182B",
      "45S rDNA" = "#2166AC",
      "5S rDNA" = "#E08214"
    ),
    name = "T2T landmark"
  ) +
  ggnewscale::new_scale_fill() +
  geom_chrlink(
    data = syri_links,
    mapping = aes(
      chr1 = Chr1, start1 = Start1, end1 = End1,
      chr2 = Chr2, start2 = Start2, end2 = End2,
      orientation = Orientation, fill = SV_Class
    ),
    curvature = 0.22,
    alpha = 0.58,
    colour = NA
  ) +
  scale_fill_manual(
    values = c(
      "SYN" = "#B8C4D0",
      "INV" = "#D6604D",
      "TRANS" = "#4393C3",
      "DUP" = "#762A83"
    ),
    labels = c(
      "SYN" = "Syntenic (SYN)",
      "INV" = "Inversion (INV)",
      "TRANS" = "Translocation (TRANS)",
      "DUP" = "Duplication (DUP)"
    ),
    name = "SyRI block"
  ) +
  geom_locus(
    data = telomeres[telomeres$Side == "left", ],
    mapping = aes(chr = Key, position = Pos),
    side = "left", gap = 0.15, shape = 21, size = 2.0,
    fill = "#FFD92F", colour = "#333333", stroke = 0.3
  ) +
  geom_locus(
    data = telomeres[telomeres$Side == "right", ],
    mapping = aes(chr = Key, position = Pos),
    side = "right", gap = 0.15, shape = 21, size = 2.0,
    fill = "#FFD92F", colour = "#333333", stroke = 0.3
  ) +
  geom_chr(
    chr = chr_key("HapA", "Chr03"), component = "axis",
    side = "left", units = "Mb", gap = 3.4, size = 2.4
  ) +
  geom_chr(
    chr = chr_key("HapB", "Chr03"), component = "axis",
    side = "right", units = "Mb", gap = 3.4, size = 2.4
  ) +
  labs(
    title = "Gap-free T2T phased haplotype architecture and structural rearrangements",
    subtitle = "SyRI rearrangements (SYN / INV / TRANS / DUP), CENH3 & rDNA arrays, telomeres, and bedGraph 5mCG profiles"
  ) +
  theme(
    plot.title = element_text(size = 10.5, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.2, colour = "#495057", hjust = 0),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(1, "mm"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.key.size = grid::unit(3.6, "mm")
  )

ggsave(
  file.path(out_dir, "t2t-haplotype-syri.png"),
  p_t2t, width = 225, height = 136, units = "mm", dpi = 300, bg = "white"
)

# ---------------------------------------------------------------------------
# 2. Plant Breeding CSSL Introgression BinMap (discrete y + geom_rug) &
#    Genetic-Physical (cM <-> Mb) Marey Composite
# ---------------------------------------------------------------------------
chr6_kar <- data.frame(
  Chr = "Chr06", Start = 0, End = 31.2e6,
  CE_start = 13.5e6, CE_end = 16.4e6
)

line_levels <- c("Recurrent (Nip)", "Donor (Kasalath)", "CSSL-01", "CSSL-02", "CSSL-03", "CSSL-04")
binmap_df <- rbind(
  data.frame(Chr = "Chr06", Start = 1, End = 31.2e6, Line = "Recurrent (Nip)", Genotype = "Recurrent"),
  data.frame(Chr = "Chr06", Start = 1, End = 31.2e6, Line = "Donor (Kasalath)", Genotype = "Donor"),
  data.frame(Chr = "Chr06", Start = c(1, 4.2e6 + 1, 11.8e6 + 1), End = c(4.2e6, 11.8e6, 31.2e6),
             Line = "CSSL-01", Genotype = c("Recurrent", "Donor", "Recurrent")),
  data.frame(Chr = "Chr06", Start = c(1, 7.6e6 + 1, 14.5e6 + 1), End = c(7.6e6, 14.5e6, 31.2e6),
             Line = "CSSL-02", Genotype = c("Recurrent", "Donor", "Recurrent")),
  data.frame(Chr = "Chr06", Start = c(1, 7.6e6 + 1, 10.6e6 + 1, 21.0e6 + 1),
             End = c(7.6e6, 10.6e6, 21.0e6, 31.2e6),
             Line = "CSSL-03", Genotype = c("Recurrent", "Het", "Recurrent", "Donor")),
  data.frame(Chr = "Chr06", Start = c(1, 18.5e6 + 1, 26.8e6 + 1), End = c(18.5e6, 26.8e6, 31.2e6),
             Line = "CSSL-04", Genotype = c("Recurrent", "Donor", "Recurrent"))
)
binmap_df$Line <- factor(binmap_df$Line, levels = line_levels)

kasp_markers <- data.frame(
  Chr = "Chr06",
  Pos = c(1.5, 3.8, 4.5, 7.4, 7.9, 9.1, 10.4, 11.6, 14.3, 18.3, 21.2, 24.5, 27.0, 29.8) * 1e6
)

set.seed(12)
scan_pos <- seq(0.5e6, 30.8e6, by = 0.5e6)
fst_df <- data.frame(
  Chr = "Chr06",
  Pos = scan_pos,
  Fst = pmin(0.85, pmax(0.03,
    0.10 + 0.58 * exp(-((scan_pos - 9.1e6) / 1.4e6)^2) +
    0.22 * exp(-((scan_pos - 23.0e6) / 2.2e6)^2) + rnorm(length(scan_pos), 0, 0.035)))
)

p_binmap <- ggideogram(
  chr6_kar,
  orientation = "horizontal",
  chromosome_width = 0.9,
  axis = TRUE, axis_side = "left", axis_units = "Mb", axis_size = 2.4,
  name_size = 2.9,
  tracks = list(
    fst = geom_track(
      side = "left", width = 2.6, gap = 0.4, limits = c(0, 0.9),
      data = fst_df,
      mapping = aes(chr = Chr, x = Pos, y = Fst),
      label = list(text = "Fst scan", size = 2.5),
      axis = list(breaks = c(0, 0.4, 0.8), size = 2.0),
      layers = list(
        geom_point(size = 0.7, colour = "#7F8C8D", alpha = 0.65),
        geom_smooth(method = "loess", formula = y ~ x, span = 0.22, se = TRUE,
          colour = "#B2182B", fill = "#F4A582", linewidth = 0.55, alpha = 0.35),
        geom_hline(yintercept = 0.45, linetype = "dashed", colour = "#B2182B", linewidth = 0.3),
        geom_vline(xintercept = 9.1e6, linetype = "dotted", colour = "#202020", linewidth = 0.35)
      )
    ),
    hap = geom_track(
      side = "right", width = 6.2, gap = 0.45,
      mapping = aes(chr = Chr),
      label = list(text = "Breeding lines", size = 2.5),
      axis = list(position = "start", size = 2.15, tick_length = 1.0),
      layers = list(
        geom_rect(data = binmap_df,
          mapping = aes(xmin = Start - 1, xmax = End, y = Line, fill = Genotype),
          height = 0.84, colour = "white", linewidth = 0.18),
        geom_hline(yintercept = 2.5, colour = "#202020", linewidth = 0.32, linetype = "dashed"),
        geom_rug(data = kasp_markers, mapping = aes(x = Pos),
          sides = "b", length = 0.10, colour = "#202020", linewidth = 0.28)
      )
    )
  )
) +
  scale_fill_manual(
    values = c("Recurrent" = "#4393C3", "Donor" = "#D6604D", "Het" = "#FDB863"),
    name = "Genotype"
  ) +
  labs(
    title = "A  CSSL introgression BinMap (discrete y), KASP rug & LOESS Fst sweep"
  ) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(size = 7.8, face = "bold"),
    legend.text = element_text(size = 7.3),
    legend.key.size = grid::unit(3.4, "mm")
  )

# Panel B: Genetic linkage group (cM) vs physical chromosome (Mb) + Marey curve
marey_kar <- data.frame(
  Genome = c("Genetic", "Physical"),
  Chr = c("LG06", "Chr06"),
  Start = c(0, 0),
  End = c(124, 31.2e6),
  CE_start = c(NA, 13.5e6),
  CE_end = c(NA, 16.4e6),
  Label = c("LG06 (Genetic, cM)", "Chr06 (Physical, Mb)")
)
marey_sem <- as_ideogram_data(
  marey_kar,
  mapping = aes(genome = Genome, chr = Chr, start = Start, end = End, label = Label),
  centromere = aes(start = CE_start, end = CE_end)
)

phys_anchors <- seq(0.8e6, 30.4e6, length.out = 28)
recomb_density <- 1.15 - 0.92 * exp(-((phys_anchors - 15.0e6) / 4.2e6)^2)
gen_cM <- cumsum(recomb_density)
gen_cM <- round((gen_cM - min(gen_cM)) / (max(gen_cM) - min(gen_cM)) * 118 + 3, 2)

anchor_links <- data.frame(
  G1 = chr_key("Genetic", "LG06"),
  Pos1 = gen_cM,
  G2 = chr_key("Physical", "Chr06"),
  Pos2 = phys_anchors
)
marey_pts <- data.frame(
  Key = chr_key("Physical", "Chr06"),
  Pos = phys_anchors,
  cM = gen_cM
)

p_marey <- ggideogram(
  marey_sem,
  orientation = "horizontal",
  scale_length = "per_chr",
  chromosome_width = 0.95,
  chromosome_gap = 8.5,
  name_size = 2.8,
  tracks = list(
    marey = geom_track(
      side = "right", width = 3.2, gap = 0.45,
      chr = chr_key("Physical", "Chr06"),
      limits = c(0, 130),
      data = marey_pts,
      mapping = aes(chr = Key, x = Pos, y = cM),
      label = list(text = "Marey (cM)", size = 2.5),
      axis = list(breaks = c(0, 60, 120), size = 2.0),
      layers = list(
        geom_point(size = 0.85, colour = "#1B7837", alpha = 0.8),
        geom_smooth(method = "loess", formula = y ~ x, span = 0.35, se = TRUE,
          colour = "#00441B", fill = "#A6DBA0", linewidth = 0.55)
      )
    )
  )
) +
  geom_chrlink(
    data = anchor_links,
    mapping = aes(chr1 = G1, position1 = Pos1, chr2 = G2, position2 = Pos2),
    type = "point", colour = "#5C6B73", alpha = 0.55, linewidth = 0.28, curvature = 0.15
  ) +
  geom_chr(
    chr = chr_key("Genetic", "LG06"), component = "axis",
    side = "left", units = "cM", breaks = seq(0, 120, by = 20), size = 2.35
  ) +
  geom_chr(
    chr = chr_key("Physical", "Chr06"), component = "axis",
    side = "right", units = "Mb", gap = 4.0, size = 2.35
  ) +
  labs(
    title = "B  Genetic-to-physical synteny (cM vs Mb) & Marey recombination curve"
  ) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0)
  )

p_breeding <- p_binmap / p_marey + patchwork::plot_layout(heights = c(1.15, 1))
ggsave(
  file.path(out_dir, "breeding-binmap-marey-composite.png"),
  p_breeding, width = 225, height = 175, units = "mm", dpi = 300, bg = "white"
)
