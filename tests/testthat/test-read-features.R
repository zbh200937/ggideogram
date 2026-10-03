write_feature_fixture <- function(lines, extension) {
  path <- tempfile(fileext = extension)
  writeLines(lines, path)
  path
}

test_that("BED and GFF features share one-based coordinates without changing source rows", {
  bed <- write_feature_fixture(c("track name=example", "01\t0\t1\tfirst\t500\t+",
                                 "01\t9\t20\tlast\t0\t-"), ".bed")
  gff <- write_feature_fixture(c("##gff-version 3",
    "01\ttest\tgene\t1\t1\t.\t+\t.\tID=first;Name=gene%20one",
    "01\ttest\texon\t10\t20\t.\t-\t.\tID=last;Parent=tx1,tx2",
    "##FASTA", ">01", "ACGT"), ".gff3")
  on.exit(unlink(c(bed, gff)))
  k <- data.frame(Chr = "01", Start = 0, End = 20)
  b <- read_chr_features(bed, karyotype = k)
  g <- read_chr_features(gff, karyotype = k)
  expect_identical(b$Chr, c("01", "01"))
  expect_equal(b$Start, c(1, 10))
  expect_equal(b$End, g$End)
  expect_equal(b$Start, g$Start)
  expect_equal(g$Name[1], "gene one")
  expect_equal(g$Parent[2], "tx1,tx2")
  expect_match(g$Attributes[1], "gene%20one", fixed = TRUE)
  expect_identical(b$Strand, c("+", "-"))
  expect_identical(attr(b, "source_coordinates"), "0-based-half-open")
  p <- ggideogram(k) + geom_chr_interval(data = b,
    ggplot2::aes(chr = Chr, start = Start, end = End))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("GTF identifiers and quoted text survive extraction", {
  gtf <- write_feature_fixture(c(
    'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "gene; one";',
    'chr1\ttest\texon\t20\t30\t.\t+\t.\tgene_id "g1"; transcript_id "t2";'), ".gtf")
  on.exit(unlink(gtf))
  x <- read_chr_features(gtf)
  expect_equal(x$gene_id, c("g1", "g1"))
  expect_equal(x$transcript_id, c("t1", "t2"))
  expect_equal(x$gene_name, c("gene; one", NA))
})

test_that("table conversion preserves metadata and validates coordinates and assembly", {
  x <- data.frame(seq = "chr1", lo = 0, hi = 10, score = 2.5)
  k <- data.frame(Chr = "1", Start = 0, End = 20, Assembly = "v1")
  mapping <- c(chr = "seq", start = "lo", end = "hi")
  converted <- as_chr_features(x, mapping,
    coordinate_system = "0-based-half-open", chr_map = c(chr1 = "1"),
    karyotype = k, assembly = "v1")
  expect_equal(converted$Chr, "1")
  expect_equal(converted$Start, 1)
  expect_equal(converted$lo, 0)
  expect_equal(converted$score, 2.5)
  expect_equal(attr(converted, "assembly"), "v1")
  expect_error(as_chr_features(x, mapping), "positive integer")
  expect_error(as_chr_features(transform(x, hi = 0), mapping,
    coordinate_system = "0-based-half-open"), "zero-width")
  expect_error(as_chr_features(transform(x, lo = 1), mapping, karyotype = k), "Unknown")
  expect_error(as_chr_features(transform(x, lo = 1, hi = 21), mapping,
    chr_map = c(chr1 = "1"), karyotype = k), "outside")
  expect_error(as_chr_features(transform(x, lo = 1), mapping,
    chr_map = c(chr1 = "1"), karyotype = k, assembly = "v2"), "assembly")
})

test_that("annotation gzip input is supported", {
  path <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(path))
  con <- gzfile(path, "wt")
  writeLines("chr1\t0\t10", con)
  close(con)
  expect_equal(read_chr_features(path)$Start, 1)
})

test_that("assemblies are checked on matched chromosomes with explicit genome keys", {
  k <- data.frame(Genome = c("A", "B"), Chr = "1", Start = 0, End = 1000,
                  Build = c("build_A", "build_B"))
  s <- as_ideogram_data(k,
    ggplot2::aes(chr = Chr, start = Start, end = End, genome = Genome, assembly = Build))
  d <- data.frame(Chr = chr_key("A", "1", "build_A"), Start = 1, End = 10)
  expect_equal(as_chr_features(d, karyotype = s, assembly = "build_A")$Chr, d$Chr)
  expect_error(as_chr_features(d, karyotype = s, assembly = "build_B"), "assembly")
  plain <- data.frame(Chr = c("A1", "B1"), Start = 0, End = 1000,
                      Assembly = c("build_A", "build_B"))
  d$Chr <- "A1"
  expect_equal(as_chr_features(d, karyotype = plain, assembly = "build_A")$Start, 1)
  expect_error(as_chr_features(d, karyotype = plain, assembly = "build_B"), "assembly")
})

test_that("closed annotations respect nonzero source plotting boundaries", {
  k <- data.frame(Chr = "1", Start = 10, End = 1000)
  d <- data.frame(Chr = "1", Start = c(11, 1000), End = c(11, 1000))
  features <- as_chr_features(d, karyotype = k)
  expect_equal(features$Start, d$Start)
  expect_error(as_chr_features(transform(d, Start = 10), karyotype = k), "outside")
  p <- ggideogram(k, tracks = track_layout(body = track("overlay", width = 1))) +
    geom_chr_fill(data = features, track = "body", ggplot2::aes(chr = Chr, start = Start, end = End))
  expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("GRanges and Seqinfo preserve metadata and use known full sequence lengths", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("Seqinfo")
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1"), start = c(1, 20), end = c(10, 30),
    strand = c("+", "-"), gene_id = c("g1", "g2")), keep.extra.columns = TRUE)
  info <- Seqinfo::Seqinfo("chr1", seqlengths = 1000, genome = "test-v1")
  GenomicRanges::seqinfo(gr) <- info
  features <- as_chr_features(gr)
  expect_equal(features$Start, c(1, 20))
  expect_equal(features$gene_id, c("g1", "g2"))
  expect_equal(as.character(features$strand), c("+", "-"))
  expect_equal(attr(features, "assembly"), "test-v1")
  expect_equal(as_ideogram_data(gr)$karyotype$.end, 1000)
  expect_equal(as_ideogram_data(info)$karyotype$.end, 1000)
  expect_error(as_chr_features(gr, assembly = "wrong"), "conflicts")
  Seqinfo::seqlengths(info) <- NA_integer_
  expect_error(as_ideogram_data(info), "Lengths")
})
