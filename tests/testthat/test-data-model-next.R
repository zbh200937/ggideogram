test_that("ideogram_data maps standard columns without discarding source data", {
  x <- as_ideogram_data(human_karyotype)

  expect_s3_class(x, "ideogram_data")
  expect_identical(x$karyotype$Chr, human_karyotype$Chr)
  expect_identical(x$karyotype$.chr, human_karyotype$Chr)
  expect_equal(x$karyotype$.start, human_karyotype$Start)
  expect_equal(x$karyotype$.end, human_karyotype$End)
  expect_equal(x$karyotype$.centromere_start, human_karyotype$CE_start)
  expect_equal(x$karyotype$.centromere_end, human_karyotype$CE_end)
  expect_null(x$cytoband)
})
test_that("semantic mappings accept arbitrary source columns and expressions", {
  raw <- data.frame(
    seqname = c("A", "B"), origin = c(1, 11), span = c(100, 50),
    waist = c(40, 20), stringsAsFactors = FALSE
  )
  x <- as_ideogram_data(
    raw,
    mapping = ggplot2::aes(chr = .data$seqname,
                           start = .data$origin,
                           end = .data$origin + .data$span),
    centromere = ggplot2::aes(start = .data$origin + .data$waist,
                              end = .data$origin + .data$waist + 5)
  )

  expect_identical(x$karyotype$seqname, raw$seqname)
  expect_equal(x$karyotype$.start, c(1, 11))
  expect_equal(x$karyotype$.end, c(101, 61))
  expect_equal(x$karyotype$.centromere_start, c(41, 31))
  expect_equal(x$karyotype$.centromere_end, c(46, 36))
})

test_that("ideogram_data validates scientific coordinate invariants", {
  expect_error(as_ideogram_data(data.frame(Chr = character(), Start = numeric(),
                                            End = numeric())),
               "no rows")
  expect_error(as_ideogram_data(data.frame(Chr = c("A", "A"), Start = 0,
                                            End = c(10, 20))),
               "unique")
  expect_error(as_ideogram_data(data.frame(Chr = "A", Start = 10, End = 10)),
               "0 <= start < end", fixed = TRUE)
  expect_error(as_ideogram_data(data.frame(Chr = "A", Start = 0, End = Inf)),
               "finite")
  expect_error(as_ideogram_data(data.frame(Chr = "A", Start = 0, End = 10,
                                            CE_start = 2)),
               "not both")
  expect_error(as_ideogram_data(data.frame(Chr = "A", Start = 0, End = 10,
                                            CE_start = 8, CE_end = 11)),
               "out of range")
  expect_error(
    as_ideogram_data(data.frame(a = "A", b = 0, c = 10),
                      mapping = ggplot2::aes(chr = .data$a, start = .data$b)),
    "missing aesthetic.*`end`"
  )
})

test_that("cytobands retain source columns and are bounded by chromosomes", {
  kar <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
  bands <- data.frame(seqname = c("A", "A", "B"), lo = c(0, 50, 0),
                      hi = c(50, 100, 80), stain = c("gneg", "gpos", "gneg"))
  x <- as_ideogram_data(
    kar,
    cytoband = bands,
    cytoband_mapping = ggplot2::aes(chr = .data$seqname,
                                    start = .data$lo, end = .data$hi)
  )

  expect_identical(x$cytoband$stain, bands$stain)
  expect_identical(x$cytoband$.chr, bands$seqname)
  expect_error(
    as_ideogram_data(kar, cytoband = transform(bands, seqname = "Z"),
                      cytoband_mapping = ggplot2::aes(
                        chr = .data$seqname, start = .data$lo, end = .data$hi)),
    "unknown chromosome"
  )
  expect_error(
    as_ideogram_data(kar, cytoband = transform(bands, hi = 1000),
                      cytoband_mapping = ggplot2::aes(
                        chr = .data$seqname, start = .data$lo, end = .data$hi)),
    "invalid"
  )
})

test_that("coercing semantic data again is idempotent", {
  x <- as_ideogram_data(human_karyotype)
  expect_identical(as_ideogram_data(x), x)
  expect_error(as_ideogram_data(x, unknown = TRUE), "Unused argument")
})

test_that("ideogram_data has an informative compact print method", {
  x <- as_ideogram_data(human_karyotype)
  out <- capture.output(print(x))
  expect_match(out[1], "<ideogram_data>", fixed = TRUE)
  expect_true(any(grepl("chromosomes: 24", out, fixed = TRUE)))
  expect_true(any(grepl("centromeres: 24", out, fixed = TRUE)))
})
