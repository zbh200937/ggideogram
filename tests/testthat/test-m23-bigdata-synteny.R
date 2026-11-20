test_that("merge_chr_links chains co-directional and inverted fragments accurately", {
  links <- data.frame(
    Chr1 = c("A1", "A1", "A1", "A1",  "A2", "A2", "A2",  "A2", "A2"),
    Start1 = c(100, 500, 1200, 500000,  100, 1100, 2100,  100, 1100),
    End1   = c(600, 1100, 1800, 501000, 1000, 2000, 3000, 1000, 2000),
    Chr2 = c("B1", "B1", "B1", "B1",  "B2", "B2", "B2",  "B2", "B2"),
    Start2 = c(200, 600, 1300, 500000,  9000, 7950, 6900,  1000, 5000),
    End2   = c(700, 1200, 1900, 501000, 9900, 8950, 7900,  1900, 5900),
    Orientation = c("+", "+", "+", "+",  "-", "-", "-",  "-", "-"),
    Type = c("SYN", "SYN", "SYN", "SYN",  "INV", "INV", "INV",  "INVTR", "INVTR"),
    stringsAsFactors = FALSE
  )

  # Separate by Type: first 3 SYN chain (overlap 500..600 -> union 100..1100 + 1200..1800 = 1001 + 601 = 1602);
  # 4th SYN is beyond max_gap; 3 INV chain downward on B2; 2 INVTR jump UP on B2 (non-monotonic for "-") so stay separate.
  m <- merge_chr_links(links, max_gap = 200, group = "Type")
  expect_identical(attr(m, "coordinates"), "1-based-closed")
  expect_equal(nrow(m), 5L)

  syn1 <- m[m$Chr1 == "A1" & m$Start1 == 100, ]
  expect_equal(syn1$End1, 1800)
  expect_equal(syn1$Start2, 200)
  expect_equal(syn1$End2, 1900)
  expect_equal(syn1$N, 3L)
  expect_equal(syn1$Span1, 1701)
  expect_equal(syn1$Span2, 1701)
  expect_equal(syn1$Aligned_bp1, (1100 - 100 + 1) + (1800 - 1200 + 1))
  expect_equal(syn1$Aligned_bp2, (1200 - 200 + 1) + (1900 - 1300 + 1))

  inv1 <- m[m$Type == "INV", ]
  expect_equal(nrow(inv1), 1L)
  expect_equal(inv1$Start1, 100)
  expect_equal(inv1$End1, 3000)
  expect_equal(inv1$Start2, 6900)
  expect_equal(inv1$End2, 9900)
  expect_equal(inv1$Orientation, "-")
  expect_equal(inv1$N, 3L)
  expect_equal(inv1$Aligned_bp1, 901 + 901 + 901)
  expect_equal(inv1$Aligned_bp2, 901 + 1001 + 1001)

  expect_equal(sum(m$Type == "INVTR"), 2L)

  # min_blocks and min_length filters
  m_filt <- merge_chr_links(links, max_gap = 200, min_blocks = 2L, min_length = 2000, group = "Type")
  expect_equal(nrow(m_filt), 1L)
  expect_equal(m_filt$Type, "INV")

  # Empty and error cases
  m_empty <- merge_chr_links(links[0, ], group = "Type")
  expect_equal(nrow(m_empty), 0L)
  expect_true(all(c("Type", "N", "Span1", "Aligned_bp1") %in% names(m_empty)))

  bad <- links[1, ]
  bad$Start1 <- 0
  expect_error(merge_chr_links(bad), "one-based closed")
})

test_that("read_chr_links vectorized parsers and unknown_chr='drop' handle PAF, SyRI and MCScanX", {
  kar <- data.frame(
    Chr = c(chr_key("Ref", c("Chr1", "Chr2")), chr_key("Qry", c("Chr1", "Chr2"))),
    Start = 0, End = 10000, stringsAsFactors = FALSE
  )

  # PAF with comments, whitespace lines, long cg:Z: tags and an unplaced scaffold
  paf_tmp <- tempfile(fileext = ".paf")
  on.exit(unlink(paf_tmp), add = TRUE)
  long_cg <- paste0("cg:Z:", paste(rep("10M1I10M1D", 400), collapse = ""))
  writeLines(c(
    "# header comment",
    "   # indented comment",
    "   ",
    paste("Chr1", 10000, 100, 1100, "+", "Chr1", 10000, 200, 1200, 950, 1000, 60, "tp:A:P", long_cg, sep = "\t"),
    paste("scaf_99", 5000, 10, 910, "-", "Chr2", 10000, 500, 1400, 880, 900, 50, long_cg, sep = "\t")
  ), paf_tmp)

  expect_error(
    read_chr_links(paf_tmp, target_genome = "Ref", query_genome = "Qry", karyotype = kar),
    "Unknown link chromosomes"
  )
  paf_drop <- read_chr_links(paf_tmp, target_genome = "Ref", query_genome = "Qry",
                             karyotype = kar, unknown_chr = "drop")
  expect_equal(nrow(paf_drop), 1L)
  expect_equal(paf_drop$Start1, 201)
  expect_equal(paf_drop$End1, 1200)
  expect_equal(paf_drop$Start2, 101)
  expect_equal(paf_drop$End2, 1100)
  expect_equal(paf_drop$Identity, 0.95)

  # Malformed short PAF still errors
  bad_paf <- tempfile(fileext = ".paf")
  on.exit(unlink(bad_paf), add = TRUE)
  writeLines("Chr1\t1000\t0\t100\t+", bad_paf)
  expect_error(read_chr_links(bad_paf, format = "paf"), "at least 12")

  # SyRI with SNP noise, comments and unknown scaffold
  syri_tmp <- tempfile(fileext = ".out")
  on.exit(unlink(syri_tmp), add = TRUE)
  writeLines(c(
    "# syri output",
    "Chr1\t100\t100\tA\tG\tChr1\t100\t100\tSNP1\tSYN1\tSNP\t-",
    "Chr1\t100\t1500\t-\t-\tChr1\t200\t1600\tSYN1\t-\tSYN\t-",
    "Chr2\t3000\t2000\t-\t-\tChr2\t4000\t5000\tINV1\t-\tINV\t-",
    "scaf_unplaced\t10\t500\t-\t-\tChr1\t10\t500\tSYN2\t-\tSYN\t-"
  ), syri_tmp)

  expect_error(
    read_chr_links(syri_tmp, format = "syri", target_genome = "Ref",
                   query_genome = "Qry", karyotype = kar),
    "Unknown link chromosomes"
  )
  syri_drop <- read_chr_links(syri_tmp, format = "syri", target_genome = "Ref",
                              query_genome = "Qry", karyotype = kar, unknown_chr = "drop")
  expect_equal(nrow(syri_drop), 2L)
  expect_equal(syri_drop$Orientation, c("+", "-"))
  expect_equal(syri_drop$Start1[2], 2000)
  expect_equal(syri_drop$End1[2], 3000)

  # MCScanX with inferred negative slope when header lacks plus/minus, plus unknown gene drop
  genes <- data.frame(
    Chr = c("Chr1", "Chr1", "Chr2", "Chr2", "scaf1"),
    Gene = c("g1", "g2", "q1", "q2", "u1"),
    Start = c(100, 1000, 8000, 2000, 10),
    End = c(300, 1200, 8200, 2200, 200),
    stringsAsFactors = FALSE
  )
  mcs_tmp <- tempfile(fileext = ".collinearity")
  on.exit(unlink(mcs_tmp), add = TRUE)
  writeLines(c(
    "## Alignment 0: score=200.0 e_value=1e-10 N=2 Chr1&Chr2",
    "  0-  0:\tg1\tq1\t1e-20",
    "  0-  1:\tg2\tq2\t1e-25",
    "## Alignment 1: score=100.0 e_value=1e-5 N=2 Chr1&scaf1 plus",
    "  1-  0:\tg1\tMISSING_GENE\t1e-5",
    "  1-  1:\tg2\tu1\t1e-5"
  ), mcs_tmp)

  expect_error(read_chr_links(mcs_tmp, genes = genes), "unknown gene")
  mcs_col <- read_chr_links(mcs_tmp, genes = genes, collapse_blocks = TRUE,
                            target_genome = "Ref", query_genome = "Qry",
                            karyotype = kar, unknown_chr = "drop")
  expect_equal(nrow(mcs_col), 1L)
  expect_equal(mcs_col$Orientation, "-")
  expect_equal(mcs_col$Start1, 100)
  expect_equal(mcs_col$End1, 1200)
  expect_equal(mcs_col$Start2, 2000)
  expect_equal(mcs_col$End2, 8200)
  expect_equal(mcs_col$N, 2L)
})

test_that("read_chr_features vectorized GTF/GFF3 attributes and unknown_chr='drop' match exact semantics", {
  kar <- data.frame(Chr = "Chr1", Start = 0, End = 5000, stringsAsFactors = FALSE)

  gtf_tmp <- tempfile(fileext = ".gtf")
  on.exit(unlink(gtf_tmp), add = TRUE)
  writeLines(c(
    "Chr1\tEns\tgene\t100\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"FT\";",
    "Chr1\tEns\tgene\t600\t900\t.\t-\t.\tgene_id \"G2\"; gene_name \"gene; one\";",
    "Chr1\tEns\tgene\t1000\t1400\t.\t+\t.\tgene_id \"G3\"; gene_name \"A\\\"B\";",
    "Chr1\tEns\tgene\t1500\t1900\t.\t+\t.\tgene_id G4; gene_name \"G4_A\"; gene_name \"G4_B\";",
    "scaf_unplaced\tEns\tgene\t10\t100\t.\t+\t.\tgene_id \"G5\"; gene_name \"DROPME\";"
  ), gtf_tmp)

  expect_error(read_chr_features(gtf_tmp, karyotype = kar), "Unknown annotation chromosomes")
  gtf_df <- read_chr_features(gtf_tmp, karyotype = kar, unknown_chr = "drop")
  expect_equal(nrow(gtf_df), 4L)
  expect_equal(gtf_df$gene_id, c("G1", "G2", "G3", "G4"))
  expect_equal(gtf_df$gene_name, c("FT", "gene; one", "A\"B", "G4_A,G4_B"))

  # Bounds validation still applies to retained chromosomes
  gtf_oob <- tempfile(fileext = ".gtf")
  on.exit(unlink(gtf_oob), add = TRUE)
  writeLines("Chr1\tEns\tgene\t100\t9999\t.\t+\t.\tgene_id \"G1\";", gtf_oob)
  expect_error(read_chr_features(gtf_oob, karyotype = kar, unknown_chr = "drop"), "outside their chromosome bounds")

  # GFF3 percent-decoding and duplicate keys
  gff_tmp <- tempfile(fileext = ".gff3")
  on.exit(unlink(gff_tmp), add = TRUE)
  writeLines(c(
    "##gff-version 3",
    "Chr1\tSRC\tgene\t100\t500\t.\t+\t.\tID=G1;Name=NormalName",
    "Chr1\tSRC\tgene\t600\t900\t.\t+\t.\tID=G2;Name=My%20Gene%3B1;Parent=P1;Parent=P2",
    "Chr1\tSRC\tgene\t950\t1200\t.\t+\t.\tID=G4;Parent=P3;Note=mid;Parent=P4",
    "scaf_x\tSRC\tgene\t1\t50\t.\t+\t.\tID=G3;Name=Drop"
  ), gff_tmp)
  gff_df <- read_chr_features(gff_tmp, karyotype = kar, unknown_chr = "drop")
  expect_equal(nrow(gff_df), 3L)
  expect_equal(gff_df$Name[1:2], c("NormalName", "My Gene;1"))
  expect_equal(gff_df$Parent[2:3], c("P1,P2", "P3,P4"))

  # as_chr_features unknown_chr = "drop"
  raw_feat <- data.frame(Chr = c("Chr1", "scaf_x"), Start = c(10, 10), End = c(50, 50))
  expect_error(as_chr_features(raw_feat, karyotype = kar), "Unknown annotation chromosomes")
  kept <- as_chr_features(raw_feat, karyotype = kar, unknown_chr = "drop")
  expect_equal(kept$Chr, "Chr1")
})

test_that("chr_layout='serial', scale_length='per_genome', row_align and Bezier S-ribbons work cleanly", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Extreme length contrast: Wheat (800 Mb), Rice (40 Mb), tiny scaffold (0.2 Mb)
  k3 <- data.frame(
    Genome = c("Wheat", "Wheat", "Rice", "Rice", "Rice", "Arabidopsis", "Arabidopsis"),
    Chr = c("1A", "2A", "Os1", "Os2", "OsTiny", "At1", "At2"),
    Start = 0,
    End = c(8e8, 7e8, 4.3e7, 3.6e7, 2e5, 3e7, 2e7),
    CE_start = c(2.5e8, 3e8, 1.6e7, 1.3e7, NA, 1.5e7, 4e6),
    CE_end   = c(2.7e8, 3.2e8, 1.8e7, 1.5e7, NA, 1.7e7, 5e6),
    stringsAsFactors = FALSE
  )
  links <- data.frame(
    Chr1 = c(chr_key("Wheat", "1A"), chr_key("Wheat", "2A"), chr_key("Rice", "Os1"), chr_key("Rice", "Os1")),
    Start1 = c(1e8, 2e8, 5e6, 2.5e7),
    End1   = c(4e8, 5e8, 1.8e7, 3.5e7),
    Chr2 = c(chr_key("Rice", "Os1"), chr_key("Rice", "Os2"), chr_key("Arabidopsis", "At1"), chr_key("Rice", "Os2")),
    Start2 = c(4e6, 8e6, 2e6, 1e7),
    End2   = c(2e7, 2.6e7, 1.4e7, 2e7),
    Orientation = c("+", "-", "+", "+"),
    stringsAsFactors = FALSE
  )

  # Must NOT error on tiny OsTiny chromosome in serial layout; per_genome equalizes species spans
  p <- ggideogram(
    k3,
    mapping = ggplot2::aes(chr = Chr, start = Start, end = End,
                           centromere_start = CE_start, centromere_end = CE_end,
                           genome = Genome),
    orientation = "horizontal",
    chr_layout = "serial",
    order_by = "genome",
    scale_length = "per_genome",
    row_align = "center",
    row_gap = 6
  ) +
    geom_chrlink(
      ggplot2::aes(chr1 = Chr1, start1 = Start1, end1 = End1,
                   chr2 = Chr2, start2 = Start2, end2 = End2,
                   orientation = Orientation),
      data = links,
      curve_style = "auto"
    )

  lay <- p$coordinates$layout
  expect_identical(lay$chr_layout, "serial")
  expect_identical(lay$scale_length, "per_genome")
  expect_false(lay$length_comparable)
  expect_identical(lay$serial$row_align, "center")
  expect_equal(lay$nrow, 3L)

  # All three rows have total display_length == max_chr_length
  row_len_sum <- tapply(lay$chrom$.display_length, lay$chrom$.row, sum)
  expect_equal(as.numeric(row_len_sum), rep(lay$max_chr_length, 3L), tolerance = 1e-8)

  # Round-trip projection / unprojection works accurately in serial layout
  pts <- data.frame(Chr = c(chr_key("Wheat", "1A"), chr_key("Rice", "Os2"), chr_key("Arabidopsis", "At2")),
                    Pos = c(3e8, 1.5e7, 1e7), stringsAsFactors = FALSE)
  pr <- project_chr_point(p, pts, "Chr", "Pos")
  unpr <- unproject_chr_point(p, pr, "Chr", ".x", ".y")
  expect_equal(unpr$.position, pts$Pos, tolerance = 1e-5)

  # Cross-row auto links use cubic S-curve bezier bounded between row y-edges
  link_layer <- utils::tail(p$layers, 1)[[1]]
  expect_identical(class(link_layer$geom)[1], "GeomPolygon")
  ld <- link_layer$data
  # Pair 1 is cross-row (Wheat 1A -> Rice Os1); y stays strictly within [row1_top, row0_bottom]
  p1 <- ld[ld$.pair_id == 1L, ]
  r0_min <- lay$chrom$.y_min[lay$chrom$.chr == chr_key("Wheat", "1A")]
  r1_max <- lay$chrom$.y_max[lay$chrom$.chr == chr_key("Rice", "Os1")]
  expect_true(nrow(p1) > 4L)
  expect_true(all(p1$.pair_y <= r0_min + 1e-8 & p1$.pair_y >= r1_max - 1e-8))
  # Pair 4 is same-row (Rice Os1 -> Rice Os2); arches outward outside the chromosome y-range
  p4 <- ld[ld$.pair_id == 4L, ]
  expect_true(nrow(p4) > 4L)
  expect_true(max(p4$.pair_y) > r1_max + 0.1 || min(p4$.pair_y) < lay$chrom$.y_min[lay$chrom$.chr == chr_key("Rice", "Os1")] - 0.1)

  # Renders cleanly without warnings (including middle-row inside-arm labels & bold genome titles)
  expect_no_warning(ggplot2::ggplotGrob(p))

  # Vertical serial layout + tracks + overlay + single-row label placement
  k1 <- data.frame(Chr = paste0("Chr", 1:4), Start = 0, End = c(4e7, 3e7, 2.5e7, 2e7))
  sig <- data.frame(Chr = rep(k1$Chr, each = 5), Pos = rep(seq(2e6, 1.8e7, length.out = 5), 4),
                    Val = seq(0.1, 2, length.out = 20))
  p_single <- ggideogram(
    k1, orientation = "horizontal", chr_layout = "serial",
    tracks = list(gwas = geom_track(side = "left", width = 4, data = sig,
      mapping = ggplot2::aes(chr = Chr, x = Pos, y = Val),
      layers = list(ggplot2::geom_point())))
  )
  expect_no_warning(ggplot2::ggplotGrob(p_single))

  p_vert <- ggideogram(k3, mapping = ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome),
                       orientation = "vertical", chr_layout = "serial", order_by = "genome",
                       scale_length = "per_genome") +
    geom_chrlink(ggplot2::aes(chr1 = Chr1, start1 = Start1, end1 = End1,
                              chr2 = Chr2, start2 = Start2, end2 = End2,
                              orientation = Orientation), data = links[1:3, ])
  expect_no_warning(ggplot2::ggplotGrob(p_vert))

  # Circular + serial must error
  expect_error(ggideogram(k1, orientation = "circular", chr_layout = "serial"),
               "linear layouts")
})

test_that("circular layout supports per_genome scaling, N=200 chromosomes and high-throughput links", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # 1. N = 200 chromosomes with default gap_angle = NULL must not fail
  k200 <- data.frame(Chr = sprintf("S%03d", 1:200), Start = 0, End = rep(1e6, 200))
  lay200 <- ideogram_layout(k200, orientation = "circular")
  expect_equal(nrow(lay200$chrom), 200L)
  expect_true(sum(lay200$circular$gap_angle[-200L]) <= 340 * 0.25 + 1e-6)
  p200 <- ggideogram(k200, orientation = "circular", show_names = FALSE)
  expect_no_warning(ggplot2::ggplotGrob(p200))

  # 2. Circular per_genome equalizes genome sweeps
  k_pg <- data.frame(
    Genome = c("Big", "Big", "Small", "Small"),
    Chr = c("B1", "B2", "S1", "S2"),
    Start = 0, End = c(5e8, 5e8, 1e7, 1e7)
  )
  lay_pg <- ideogram_layout(
    as_ideogram_data(k_pg, mapping = ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome)),
    orientation = "circular", order_by = "genome", scale_length = "per_genome"
  )
  expect_false(lay_pg$length_comparable)
  sweeps <- tapply(lay_pg$chrom$.arc_end - lay_pg$chrom$.arc_start, lay_pg$chrom$.genome, sum)
  expect_equal(unname(sweeps[["Big"]]), unname(sweeps[["Small"]]), tolerance = 1e-6)

  # 3. Large circular links (>200) prune columns, preserve source_pairs, keep native GeomPolygon/GeomSegment,
  #    and support unknown_chr = "drop" for both direct links and chr_nodes()
  k5 <- data.frame(Chr = paste0("Chr", 1:5), Start = 0, End = 1e7)
  n_links <- 250L
  big_links <- data.frame(
    Chr1 = rep(paste0("Chr", 1:5), length.out = n_links),
    Start1 = round(seq(1e5, 8e6, length.out = n_links)),
    End1 = round(seq(2e5, 8.1e6, length.out = n_links)),
    Chr2 = rep(paste0("Chr", c(2, 3, 4, 5, 1)), length.out = n_links),
    Start2 = round(seq(1e5, 8e6, length.out = n_links)),
    End2 = round(seq(2e5, 8.1e6, length.out = n_links)),
    Orientation = rep(c("+", "-"), length.out = n_links),
    ExtraMeta = paste0("meta_", seq_len(n_links)),
    stringsAsFactors = FALSE
  )
  with_unknown <- rbind(
    big_links,
    data.frame(Chr1 = "ScafUnplaced", Start1 = 10, End1 = 200,
               Chr2 = "Chr1", Start2 = 10, End2 = 200,
               Orientation = "+", ExtraMeta = "drop", stringsAsFactors = FALSE)
  )

  p_circ <- ggideogram(k5, orientation = "circular", show_names = FALSE) +
    geom_chrlink(
      ggplot2::aes(chr1 = Chr1, start1 = Start1, end1 = End1,
                   chr2 = Chr2, start2 = Start2, end2 = End2,
                   orientation = Orientation, fill = Orientation),
      data = with_unknown,
      unknown_chr = "drop"
    )
  poly_layer <- utils::tail(p_circ$layers, 1)[[1]]
  expect_identical(class(poly_layer$geom)[1], "GeomPolygon")
  expect_false("ExtraMeta" %in% names(poly_layer$data))
  expect_equal(nrow(attr(poly_layer$data, "source_pairs")), n_links)
  expect_no_warning(ggplot2::ggplotGrob(p_circ))

  # Point arcs with chr_nodes() and unknown_chr = "drop"
  nodes_df <- data.frame(
    ID = c("N1", "N2", "N_UNK"),
    Chr = c("Chr1", "Chr3", "ScafX"),
    Pos = c(2e6, 5e6, 100),
    stringsAsFactors = FALSE
  )
  edges_df <- data.frame(From = c("N1", "N1"), To = c("N2", "N_UNK"), stringsAsFactors = FALSE)
  p_pt <- ggideogram(k5, orientation = "circular", show_names = FALSE) +
    geom_chrlink(
      ggplot2::aes(from = From, to = To),
      data = edges_df,
      nodes = chr_nodes(nodes_df, ggplot2::aes(id = ID, chr = Chr, position = Pos)),
      unknown_chr = "drop"
    )
  seg_layer <- utils::tail(p_pt$layers, 1)[[1]]
  expect_identical(class(seg_layer$geom)[1], "GeomSegment")
  expect_equal(length(unique(seg_layer$data$.pair_id)), 1L)
  expect_equal(unique(seg_layer$data$To), "N2")
  expect_no_warning(ggplot2::ggplotGrob(p_pt))
})

test_that("bin_genome fast count and grouped coverage paths match exact interval semantics", {
  k <- data.frame(Chr = c("Chr1", "Chr2"), Start = 0, End = c(1000, 1000))
  feats <- data.frame(
    Chr = c("Chr1", "Chr1", "Chr1", "Chr2", "Chr2"),
    Start = c(10, 50, 450, 100, 600),
    End   = c(100, 150, 550, 200, 700),
    Grp   = c("A", "A", "B", "A", NA),
    stringsAsFactors = FALSE
  )

  b_default <- bin_genome(feats, k, window = 500)
  b_count <- bin_genome(feats, k, window = 500, method = "count")
  expect_equal(b_default$Value, b_count$Value)
  expect_equal(b_default$Value, c(3, 0, 1, 1))

  # Grouped coverage with overlap (10..100 U 50..150 = 141 bp in bin [0, 500]) and NA group preserved
  b_cov <- bin_genome(feats, k, window = 500, method = "coverage", group = "Grp")
  a_bin1 <- b_cov[b_cov$Chr == "Chr1" & b_cov$Start == 1 & !is.na(b_cov$Grp) & b_cov$Grp == "A", ]
  expect_equal(a_bin1$Value, 141 / 500, tolerance = 1e-10)
  b_bin2 <- b_cov[b_cov$Chr == "Chr1" & b_cov$Start == 501 & !is.na(b_cov$Grp) & b_cov$Grp == "B", ]
  expect_equal(b_bin2$Value, 50 / 500, tolerance = 1e-10)
  na_bin2 <- b_cov[b_cov$Chr == "Chr2" & b_cov$Start == 501 & is.na(b_cov$Grp), ]
  expect_equal(na_bin2$Value, 101 / 500, tolerance = 1e-10)
})
