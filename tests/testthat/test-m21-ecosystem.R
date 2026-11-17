test_that("retired legacy shims are removed from namespace exports", {
  exports <- getNamespaceExports("ggideogram")
  expect_false("ideogram" %in% exports)
  expect_false("GFFex" %in% exports)
  expect_false("save_ideogram" %in% exports)
  expect_false("as_ideogram_grob" %in% exports)
  expect_true(all(c("ggideogram", "resolve_gene_models", "read_chr_links",
                    "read_chr_features", "track_layout", "view_chr_data") %in% exports))
})

test_that("resolve_gene_models handles GFF3 multi-parent, non-coding, direct-gene and GTF inputs", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  gff_path <- system.file("extdata", "arabidopsis-first-genes.gff3", package = "ggideogram")
  feat <- read_chr_features(gff_path)
  models <- resolve_gene_models(feat)
  expect_true(all(c("Transcript", "Gene", "Symbol", "Part") %in% names(models)))
  expect_false(any(grepl("^gene:", models$Gene)))
  expect_false(any(grepl("^transcript:", models$Transcript)))
  expect_true("NAC001" %in% models$Symbol || "AT1G01010" %in% models$Gene)
  expect_true(all(models$Part[models$Type == "CDS"] == "CDS"))
  expect_true(all(models$Part[models$Type == "five_prime_UTR"] == "UTR"))
  expect_true(all(is.na(models$Part[models$Type == "exon"])))

  # Direct two-level gene -> exon/CDS + lnc_RNA + shared exon
  d_gff <- data.frame(
    Chr = "Chr1",
    Start = c(100, 100, 120, 500, 500, 520, 900, 910),
    End = c(300, 300, 200, 800, 800, 700, 1200, 1150),
    Type = c("gene", "lnc_RNA", "exon", "gene", "mRNA", "CDS", "gene", "exon"),
    Strand = c("+", "+", "+", "-", "-", "-", "+", "+"),
    ID = c("gene:G1", "transcript:T1", "E1", "gene:G2", "transcript:T2", "C1", "G3", "E3"),
    Parent = c(NA, "gene:G1", "transcript:T1,transcript:T2", NA, "gene:G2", "transcript:T2", NA, "G3"),
    Name = c("LNC1", NA, NA, "PROT2", NA, NA, "PROK3", NA),
    stringsAsFactors = FALSE
  )
  res_gff <- resolve_gene_models(d_gff)
  e1_rows <- res_gff[!is.na(res_gff$ID) & res_gff$ID == "E1", ]
  expect_equal(nrow(e1_rows), 2L)
  expect_equal(e1_rows$Transcript, c("T1", "T2"))
  expect_equal(e1_rows$Gene, c("G1", "G2"))
  expect_equal(e1_rows$Symbol, c("LNC1", "PROT2"))
  e3_rows <- res_gff[!is.na(res_gff$ID) & res_gff$ID == "E3", ]
  expect_equal(nrow(e3_rows), 1L)
  expect_equal(e3_rows$Gene, "G3")
  expect_equal(e3_rows$Symbol, "PROK3")

  # GTF input with 5UTR / 3UTR
  gtf_tmp <- tempfile(fileext = ".gtf")
  on.exit(unlink(gtf_tmp), add = TRUE)
  writeLines(c(
    "Chr1\tEns\tgene\t100\t900\t.\t+\t.\tgene_id \"G10\"; gene_name \"FT\";",
    "Chr1\tEns\ttranscript\t100\t900\t.\t+\t.\tgene_id \"G10\"; transcript_id \"G10.1\"; gene_name \"FT\";",
    "Chr1\tEns\texon\t100\t350\t.\t+\t.\tgene_id \"G10\"; transcript_id \"G10.1\"; gene_name \"FT\";",
    "Chr1\tEns\t5UTR\t100\t150\t.\t+\t.\tgene_id \"G10\"; transcript_id \"G10.1\"; gene_name \"FT\";",
    "Chr1\tEns\tCDS\t151\t350\t.\t+\t.\tgene_id \"G10\"; transcript_id \"G10.1\"; gene_name \"FT\";"
  ), gtf_tmp)
  gtf_feat <- read_chr_features(gtf_tmp)
  gtf_mod <- resolve_gene_models(gtf_feat)
  expect_false("gene" %in% gtf_mod$Type)
  expect_equal(unique(gtf_mod$Symbol), "FT")
  expect_equal(gtf_mod$Part[gtf_mod$Type == "5UTR"], "UTR")

  kar <- data.frame(Chr = "Chr1", Start = 0, End = 2000)
  p <- ggideogram(kar, orientation = "horizontal",
    tracks = list(g = geom_track(side = "right", data = gtf_mod,
      mapping = ggplot2::aes(chr = Chr, start = Start, end = End, type = Type,
        strand = Strand, transcript = Transcript, label = Symbol, fill = Part),
      layers = list(geom_genemodel(mode = "transcript")))))
  expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("read_chr_features supports bedGraph including negative and decimal values", {
  bg_tmp <- tempfile(fileext = ".bedgraph.gz")
  on.exit(unlink(bg_tmp), add = TRUE)
  gz <- gzfile(bg_tmp, "wt")
  writeLines(c(
    "# comment line",
    "browser position Chr1:1-1000",
    "track type=bedGraph name=\"log2ChIP\"",
    "Chr1\t0\t500\t-1.25",
    "Chr1\t500\t1000\t2.75"
  ), gz)
  close(gz)

  kar <- data.frame(Chr = "Chr1", Start = 0, End = 1000)
  bg <- read_chr_features(bg_tmp, karyotype = kar)
  expect_equal(bg$Start, c(1, 501))
  expect_equal(bg$End, c(500, 1000))
  expect_equal(bg$Value, c(-1.25, 2.75))
  expect_equal(bg$Mid, c(250, 750))

  bad_bg <- tempfile(fileext = ".bg")
  on.exit(unlink(bad_bg), add = TRUE)
  writeLines("Chr1\t0\t100\tnot_a_number", bad_bg)
  expect_error(read_chr_features(bad_bg), "finite numeric")
})

test_that("read_chr_links parses PAF, SyRI and MCScanX with exact 1-based closed coordinates", {
  # 1. PAF (.paf.gz)
  paf_tmp <- tempfile(fileext = ".paf.gz")
  on.exit(unlink(paf_tmp), add = TRUE)
  gz <- gzfile(paf_tmp, "wt")
  writeLines(c(
    "QChr1\t10000\t0\t2000\t+\tTChr1\t12000\t500\t2500\t1900\t2000\t60\ttp:A:P",
    "QChr1\t10000\t3000\t4500\t-\tTChr1\t12000\t4000\t5500\t1400\t1500\t45\ttp:A:P",
    "QChr1\t10000\t5000\t5100\t+\tTChr1\t12000\t6000\t6100\t90\t100\t10\ttp:A:S"
  ), gz)
  close(gz)

  kar_pair <- as_ideogram_data(
    data.frame(Genome = c("Ref", "Qry"), Chr = c("Chr1", "Chr1"), Start = 0, End = c(12000, 10000)),
    mapping = ggplot2::aes(genome = Genome, chr = Chr, start = Start, end = End)
  )
  paf <- read_chr_links(
    paf_tmp, min_length = 500, min_mapq = 30,
    chr_map = c(TChr1 = "Chr1", QChr1 = "Chr1"),
    target_genome = "Ref", query_genome = "Qry", karyotype = kar_pair
  )
  expect_equal(nrow(paf), 2L)
  expect_equal(paf$Start1, c(501, 4001))
  expect_equal(paf$End1, c(2500, 5500))
  expect_equal(paf$Start2, c(1, 3001))
  expect_equal(paf$End2, c(2000, 4500))
  expect_equal(paf$Orientation, c("+", "-"))
  expect_equal(paf$Identity, c(1900 / 2000, 1400 / 1500))

  # 2. SyRI (.out) including reversed query coordinates on INV and unaligned '-' rows
  syri_tmp <- tempfile(fileext = ".out")
  on.exit(unlink(syri_tmp), add = TRUE)
  writeLines(c(
    "chrA\tstartA\tendA\tseqA\tseqB\tchrB\tstartB\tendB\tid\tparent\ttype\tcopy",
    "Chr1\t1\t3000\t-\t-\tChr1\t1\t2900\tSYN1\t-\tSYN\t-",
    "Chr1\t3500\t6000\t-\t-\tChr1\t5800\t3300\tINV1\t-\tINV\t-",
    "Chr1\t6001\t7000\t-\t-\t-\t-\t-\tNOTAL1\t-\tNOTAL\t-"
  ), syri_tmp)
  syri <- read_chr_links(syri_tmp, format = "auto",
    target_genome = "Ref", query_genome = "Qry", karyotype = kar_pair)
  expect_equal(nrow(syri), 2L)
  expect_equal(syri$Start2, c(1, 3300))
  expect_equal(syri$End2, c(2900, 5800))
  expect_equal(syri$Orientation, c("+", "-"))

  # 3. MCScanX (.collinearity + 4-column .gff)
  col_tmp <- tempfile(fileext = ".collinearity")
  gff4_tmp <- tempfile(fileext = ".gff")
  on.exit(unlink(c(col_tmp, gff4_tmp)), add = TRUE)
  writeLines(c(
    "Chr1\tG_A1\t100\t400",
    "Chr1\tG_A2\t600\t900",
    "Chr1\tG_B1\t1200\t1500",
    "Chr1\tG_B2\t200\t500"
  ), gff4_tmp)
  writeLines(c(
    "############### Parameters ###############",
    "## Alignment 0: score=180.0 e_value=1e-20 N=2 Chr1&Chr1 minus",
    "  0-  0:\tG_A1\tG_B1\t1e-12",
    "  0-  1:\tG_A2\tG_B2\t2e-10"
  ), col_tmp)

  pairs_mc <- read_chr_links(col_tmp, genes = gff4_tmp, collapse_blocks = FALSE)
  expect_equal(nrow(pairs_mc), 2L)
  expect_equal(pairs_mc$Orientation, c("-", "-"))

  blocks_mc <- read_chr_links(col_tmp, genes = gff4_tmp, collapse_blocks = TRUE)
  expect_equal(nrow(blocks_mc), 1L)
  expect_equal(blocks_mc$Start1, 100)
  expect_equal(blocks_mc$End1, 900)
  expect_equal(blocks_mc$Start2, 200)
  expect_equal(blocks_mc$End2, 1500)
  expect_equal(blocks_mc$Orientation, "-")
  expect_equal(blocks_mc$N, 2L)
})

test_that("geom_track natively supports discrete y channels, category axes and sparse chromosomes", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  kar <- data.frame(Chr = c("Chr1", "Chr2"), Start = 0, End = c(100, 80))
  df <- data.frame(
    Chr = c("Chr1", "Chr1", "Chr1", "Chr2", "Chr2"),
    Mid = c(25, 25, 25, 40, 40),
    Start = c(1, 1, 1, 21, 21),
    End = c(50, 50, 50, 60, 60),
    Line = factor(c("P1", "NIL_A", "NIL_B", "P1", "NIL_B"),
                  levels = c("P1", "NIL_A", "NIL_B")),
    State = c("AA", "BB", "AB", "AA", "BB")
  )
  for (ori in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(kar, orientation = ori,
      opening_angle = if (ori == "circular") 50 else NULL,
      tracks = list(hap = geom_track(
        side = if (ori == "circular") "outer" else "right",
        width = 4, data = df,
        mapping = ggplot2::aes(chr = Chr, x = Mid, width = 50, y = Line, fill = State),
        axis = TRUE, label = "Lines",
        layers = list(
          ggplot2::geom_tile(),
          ggplot2::geom_hline(yintercept = 1.5, linetype = "dashed")
        )
      )))
    expect_no_warning(g <- ggplot2::ggplotGrob(p))
    expect_s3_class(g, "gtable")
    b <- ggplot2::ggplot_build(p)
    lbls <- unlist(lapply(b$data, function(d) if ("label" %in% names(d)) as.character(d$label) else NULL))
    expect_true(all(c("P1", "NIL_A", "NIL_B") %in% lbls))
  }

  # geom_rect on discrete y auto-computes lane ymin/ymax
  p_rect <- ggideogram(kar, orientation = "horizontal",
    tracks = list(hap = geom_track(side = "right", width = 4, data = df,
      mapping = ggplot2::aes(chr = Chr, xmin = Start - 1, xmax = End, y = Line, fill = State),
      geom = ggplot2::geom_rect())))
  expect_no_warning(ggplot2::ggplotGrob(p_rect))

  # Mixing continuous and discrete y in the same track is rejected clearly
  expect_error(
    ggideogram(kar, tracks = list(bad = geom_track(side = "right", data = df,
      mapping = ggplot2::aes(chr = Chr),
      layers = list(
        ggplot2::geom_point(ggplot2::aes(x = Mid, y = Line)),
        ggplot2::geom_point(ggplot2::aes(x = Mid, y = Start))
      )))),
    "cannot mix discrete and continuous"
  )
})

test_that("geom_smooth, geom_rug and cM / none axis units work across orientations", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  kar <- data.frame(Chr = c("LG1", "LG2"), Start = 0, End = c(120, 90))
  set.seed(7)
  pts <- data.frame(
    Chr = rep(c("LG1", "LG2"), each = 20),
    Pos = c(seq(3, 117, length.out = 20), seq(3, 87, length.out = 20)),
    Val = c(seq(1, 15, length.out = 20) + rnorm(20, 0, 0.5),
            seq(14, 2, length.out = 20) + rnorm(20, 0, 0.5))
  )
  for (ori in c("vertical", "horizontal", "circular")) {
    p <- ggideogram(kar, orientation = ori, reverse_chr = "LG2",
      axis = ori != "circular", axis_units = "cM",
      opening_angle = if (ori == "circular") 40 else NULL,
      tracks = list(sm = geom_track(
        side = if (ori == "circular") "outer" else "right",
        width = 3.5, limits = c(0, 18), data = pts,
        mapping = ggplot2::aes(chr = Chr, x = Pos, y = Val),
        layers = list(
          ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE),
          ggplot2::geom_rug(sides = "bt")
        )
      )))
    expect_no_warning(ggplot2::ggplotGrob(p))
    b <- ggplot2::ggplot_build(p)
    sm_data <- b$data[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomSmooth"), logical(1)))]]
    expect_equal(sort(unique(sm_data$chr)), c("LG1", "LG2"))
  }

  # Invalid rug sides (only 'lr' without 'b' or 't') raises clear error
  expect_error(
    ggideogram(kar, tracks = list(r = geom_track(side = "right", data = pts,
      mapping = ggplot2::aes(chr = Chr, x = Pos),
      geom = ggplot2::geom_rug(sides = "lr")))),
    "sides"
  )

  # Unitless axis ('none') formats plain numbers without trailing space
  p_u <- ggideogram(kar, orientation = "horizontal", axis = TRUE, axis_units = "none")
  bu <- ggplot2::ggplot_build(p_u)
  u_lbls <- unlist(lapply(bu$data, function(d) if ("label" %in% names(d)) d$label else NULL))
  expect_true("20" %in% u_lbls)
  expect_false(any(grepl(" $", u_lbls)))
})
