#' Read chromosome lengths from a sequence index or size table
#'
#' Reads the sequence name and length columns of a headerless, tab-separated
#' `chrom.sizes` table or `.fai` index. Names and input order are preserved.
#' No centromere positions are inferred from sequence lengths.
#'
#' @param file A local file path (optionally gzip-compressed) or a data frame
#'   with positional columns in the same format.
#' @param format Input format. `"auto"` recognizes two-column `chrom.sizes`
#'   tables and five- or six-column FAI indexes. Only the first two FAI
#'   columns are used; the sequence file itself is not needed.
#' @param chr Optional chromosome selection and order. Names must be unique
#'   and present in the input.
#'
#' @return A data frame with `Chr`, `Start = 0`, and `End` equal to the input
#'   length, suitable for [ggideogram()] or [as_ideogram_data()]. `Start = 0`
#'   is the whole-chromosome plotting boundary, not a conversion of any
#'   annotation coordinates. Lengths must be positive whole numbers below
#'   `2^53`, the consecutive-integer precision limit of R doubles.
#' @examples
#' sizes <- data.frame(sequence = c("Chr1", "Chr2"), length = c(1000, 700))
#' karyotype <- read_karyotype(sizes)
#' ggideogram(karyotype)
#' read_karyotype(sizes, chr = c("Chr2", "Chr1"))
#' @seealso [read_cytoband()], [ideogram_layout()]
#' @export
read_karyotype <- function(file, format = c("auto", "chrom.sizes", "fai"),
                           chr = NULL) {
  format <- match.arg(format)
  input <- if (is.data.frame(file)) {
    file
  } else {
    if (!is.character(file) || length(file) != 1L || is.na(file)) {
      stopf("`file` must be a file path or data frame.")
    }
    utils::read.table(file, sep = "\t", header = FALSE, quote = "",
      comment.char = "", colClasses = "character", na.strings = NULL,
      stringsAsFactors = FALSE)
  }
  if (!nrow(input)) stopf("The chromosome length table has no rows.")
  columns <- ncol(input)
  allowed <- switch(format, auto = c(2L, 5L, 6L), chrom.sizes = 2L, fai = c(5L, 6L))
  if (!columns %in% allowed) {
    stopf("Format `%s` requires %s columns; got %d.", format,
          paste(allowed, collapse = " or "), columns)
  }
  names <- validate_chr(input[[1]], "sequence names")
  if (any(!nzchar(trimws(names)))) stopf("Sequence names must not be blank.")
  if (anyDuplicated(names)) {
    stopf("Sequence names must be unique; repeated: %s.",
          format_chr_rows(unique(names[duplicated(names)])))
  }
  lengths <- suppressWarnings(as.numeric(as.character(input[[2]])))
  bad <- !is.finite(lengths) | lengths <= 0 | lengths >= 2^53 |
    lengths != floor(lengths)
  if (any(bad)) {
    stopf("Lengths must be positive finite integers below 2^53; invalid for %s.",
          format_chr_rows(names[bad]))
  }
  result <- data.frame(Chr = names, Start = 0, End = lengths,
                       stringsAsFactors = FALSE)
  if (!is.null(chr)) {
    chr <- validate_chr(chr, "chr")
    if (!length(chr) || anyDuplicated(chr)) {
      stopf("`chr` must be a non-empty selection of unique names.")
    }
    unknown <- setdiff(chr, names)
    if (length(unknown)) {
      stopf("Not in the chromosome length table: %s.", format_chr_rows(unknown))
    }
    result <- result[match(chr, names), , drop = FALSE]
    rownames(result) <- NULL
  }
  result
}
