test_that("size files preserve sequence identities, lengths and explicit order", {
  file <- tempfile(fileext = ".chrom.sizes")
  on.exit(unlink(file))
  writeLines(c("01\t1000", "10\t600", "02\t800", "scaffold_1\t40"), file)
  k <- read_karyotype(file)
  expect_identical(k$Chr, c("01", "10", "02", "scaffold_1"))
  expect_equal(k$Start, rep(0, 4))
  expect_equal(k$End, c(1000, 600, 800, 40))
  selected <- read_karyotype(file, chr = c("02", "01"))
  expect_identical(selected$Chr, c("02", "01"))
  expect_equal(selected$End, c(800, 1000))
  expect_false(any(grepl("CE_", names(selected))))
  expect_no_error(ggplot2::ggplotGrob(ggideogram(selected)))
})

test_that("FAI and compressed size tables have equivalent layouts", {
  sizes <- data.frame(chr = c("one", "two"), length = c(66, 28))
  fai <- data.frame(sizes, offset = c(5, 98), bases = c(30, 14), bytes = c(31, 15))
  file <- tempfile(fileext = ".fai")
  compressed <- tempfile(fileext = ".chrom.sizes.gz")
  on.exit(unlink(c(file, compressed)))
  utils::write.table(fai, file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  con <- gzfile(compressed, "wt")
  writeLines(c("one\t66", "two\t28"), con)
  close(con)
  expect_equal(read_karyotype(file), read_karyotype(sizes))
  expect_equal(read_karyotype(compressed), read_karyotype(sizes))
  expect_equal(read_karyotype(transform(fai, qualoffset = c(79, 188))),
               read_karyotype(sizes))
  expect_equal(ideogram_layout(read_karyotype(file))$chrom,
               ideogram_layout(read_karyotype(sizes))$chrom)
})

test_that("bundled chromosome lengths round trip without altering data", {
  for (k in list(human_karyotype, liriodendron_karyotype)) {
    restored <- read_karyotype(k[c("Chr", "End")])
    expect_identical(restored$Chr, as.character(k$Chr))
    expect_equal(restored$End, k$End)
    expect_equal(restored$Start, k$Start)
    expect_no_error(as_ideogram_data(restored))
  }
})

test_that("ambiguous formats and invalid chromosome lengths are rejected", {
  sizes <- data.frame(chr = "A", length = 100)
  expect_error(read_karyotype(transform(sizes, extra = 0)), "columns")
  expect_error(read_karyotype(sizes, format = "fai"), "columns")
  expect_error(read_karyotype(sizes[FALSE, ]), "no rows")
  expect_error(read_karyotype(rbind(sizes, sizes)), "unique")
  for (value in c("0", "-1", "1.5", "unknown", "Inf", NA, "9007199254740992")) {
    expect_error(read_karyotype(data.frame(chr = "A", length = value)), "Lengths")
  }
  expect_error(read_karyotype(data.frame(chr = " ", length = 100)), "blank")
  expect_error(read_karyotype(sizes, chr = "B"), "Not in")
  expect_error(read_karyotype(sizes, chr = c("A", "A")), "unique")
  expect_error(read_karyotype(sizes, chr = character()), "non-empty")
})
