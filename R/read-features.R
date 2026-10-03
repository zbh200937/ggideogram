#' Read genomic annotations in BED, GFF3 or GTF format
#'
#' BED3--6 starts are converted from zero-based, half-open coordinates to
#' one-based, closed coordinates (`Start = chromStart + 1`, `End = chromEnd`).
#' GFF3 and GTF already use one-based, closed coordinates. Zero-width BED
#' insertion sites need an explicit point representation and are rejected.
#' Original GFF/GTF attributes are retained in `Attributes`; common identifiers
#' are also extracted without collapsing rows or parent relationships.
#'
#' @param file Local annotation file, optionally gzip-compressed.
#' @param format Format; `"auto"` uses the file extension (ignoring `.gz`).
#' @param karyotype Optional karyotype or `ideogram_data` used to validate
#'   sequence names, interval bounds and assembly metadata.
#' @param chr_map Optional named character vector, `c(old_name = "new_name")`.
#'   Unlisted names are retained.
#' @param assembly Optional assembly identifier. A conflicting known assembly
#'   for a matched chromosome in `karyotype` is an error. Mapped assembly
#'   metadata in an `ideogram_data` object is used; unknown versions are not
#'   inferred.
#'
#' @return A data frame with `Chr`, `Start`, `End` and annotation columns.
#'   Attributes `source_coordinates` and `assembly` record input conventions.
#' @seealso [as_chr_features()], [read_karyotype()]
#' @export
read_chr_features <- function(file, format = c("auto", "bed", "gff3", "gtf"),
                              karyotype = NULL, chr_map = NULL, assembly = NULL) {
  format <- match.arg(format)
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stopf("`file` must be one local annotation file path.")
  }
  if (format == "auto") {
    extension <- tolower(tools::file_ext(sub("\\.gz$", "", file, ignore.case = TRUE)))
    format <- switch(extension, bed = "bed", gff = "gff3", gff3 = "gff3", gtf = "gtf",
                     stopf("Cannot infer annotation format; supply `format`."))
  }
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  fasta <- match("##FASTA", trimws(lines))
  if (!is.na(fasta)) lines <- utils::head(lines, fasta - 1L)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  if (format == "bed") lines <- lines[!grepl("^(track|browser)\\s", lines)]
  if (!length(lines)) stopf("The annotation file contains no features.")
  input <- utils::read.table(text = lines, sep = "\t", header = FALSE,
    quote = "", comment.char = "", colClasses = "character", na.strings = NULL)
  if (format == "bed") {
    if (!ncol(input) %in% 3:6) stopf("BED input must have 3 to 6 columns.")
    names(input) <- c("Chr", "Start", "End", "Name", "Score", "Strand")[seq_len(ncol(input))]
    system <- "0-based-half-open"
  } else {
    if (ncol(input) != 9L) stopf("GFF3/GTF input must have 9 columns.")
    names(input) <- c("Chr", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
    keys <- if (format == "gff3") c("ID", "Parent", "Name") else c("gene_id", "transcript_id", "gene_name")
    for (key in keys) input[[key]] <- feature_attribute(input$Attributes, key, format)
    system <- "1-based-closed"
  }
  input$Start <- suppressWarnings(as.numeric(input$Start))
  input$End <- suppressWarnings(as.numeric(input$End))
  as_chr_features(input, coordinate_system = system, karyotype = karyotype,
                  chr_map = chr_map, assembly = assembly)
}

feature_attribute <- function(attributes, key, format) {
  vapply(attributes, function(value) {
    # scan respects quoted semicolons in GTF attribute values.
    fields <- if (format == "gtf") {
      scan(text = value, what = character(), sep = ";", quote = "\"", quiet = TRUE)
    } else strsplit(value, ";", fixed = TRUE)[[1]]
    fields <- trimws(fields)
    pattern <- if (format == "gff3") paste0("^", key, "=") else paste0("^", key, "\\s+")
    matches <- fields[grepl(pattern, fields)]
    if (!length(matches)) return(NA_character_)
    result <- sub(pattern, "", matches)
    if (format == "gff3") result <- utils::URLdecode(result)
    paste(result, collapse = ",")
  }, character(1), USE.NAMES = FALSE)
}

#' Convert an annotation table or GRanges to chromosome features
#'
#' Data-frame inputs retain their original columns, with the mapped coordinates
#' copied to `Chr`, `Start`, `End`. Coordinate conversion is explicit. GRanges
#' inputs retain their metadata and use their native one-based coordinates.
#'
#' @param x Data frame or GRanges object.
#' @param ... Arguments for the selected method.
#' @inheritParams read_chr_features
#' @return A data frame in one-based, closed coordinates.
#' @export
as_chr_features <- function(x, ...) UseMethod("as_chr_features")

#' @rdname as_chr_features
#' @param mapping Named vector connecting `chr`, `start`, `end` to columns.
#' @param coordinate_system Coordinate convention of the data-frame input.
#' @export
as_chr_features.data.frame <- function(x,
    mapping = c(chr = "Chr", start = "Start", end = "End"),
    coordinate_system = c("1-based-closed", "0-based-half-open"),
    karyotype = NULL, chr_map = NULL, assembly = NULL, ...) {
  check_unused_args(...)
  coordinate_system <- match.arg(coordinate_system)
  if (!is.character(mapping) || !setequal(names(mapping), c("chr", "start", "end")) ||
      length(mapping) != 3L || anyNA(mapping) || !all(mapping %in% names(x))) {
    stopf("`mapping` must map chr, start and end to existing columns.")
  }
  result <- x
  result$Chr <- validate_chr(x[[mapping[["chr"]]]], "chr")
  result$Start <- validate_coordinate(x[[mapping[["start"]]]], "start")
  result$End <- validate_coordinate(x[[mapping[["end"]]]], "end")
  if (coordinate_system == "0-based-half-open") result$Start <- result$Start + 1
  if (any(result$Start < 1 | result$End < result$Start |
          result$Start != floor(result$Start) | result$End != floor(result$End))) {
    stopf("Features require positive integer coordinates with Start <= End; zero-width BED is not an interval.")
  }
  if (!is.null(chr_map)) {
    if (!is.character(chr_map) || is.null(names(chr_map)) || anyNA(chr_map) ||
        anyNA(names(chr_map)) || anyDuplicated(names(chr_map)) ||
        any(!nzchar(names(chr_map))) || any(!nzchar(chr_map))) {
      stopf("`chr_map` must be a named character vector with unique source names.")
    }
    index <- match(result$Chr, names(chr_map))
    result$Chr[!is.na(index)] <- unname(chr_map[index[!is.na(index)]])
  }
  if (!is.null(assembly) && (!is.character(assembly) || length(assembly) != 1L ||
                             is.na(assembly) || !nzchar(assembly))) {
    stopf("`assembly` must be NULL or one non-empty identifier.")
  }
  if (!is.null(karyotype)) {
    k <- as_ideogram_data(karyotype)$karyotype
    index <- match(result$Chr, k$.chr)
    if (anyNA(index)) stopf("Unknown annotation chromosomes: %s.",
                            format_chr_rows(unique(result$Chr[is.na(index)])))
    if (any(result$Start - 1 < k$.start[index] | result$End > k$.end[index])) {
      stopf("Annotation intervals are outside their chromosome bounds.")
    }
    known <- if (".assembly" %in% names(k)) {
      unique(stats::na.omit(k$.assembly[index]))
    } else if ("Assembly" %in% names(k)) {
      unique(stats::na.omit(k$Assembly[index]))
    } else attr(karyotype, "assembly")
    if (!is.null(assembly) && length(known) && any(known != assembly)) {
      stopf("Annotation assembly `%s` does not match the karyotype.", assembly)
    }
  }
  attr(result, "source_coordinates") <- coordinate_system
  attr(result, "assembly") <- assembly
  result
}

#' @rdname as_chr_features
#' @export
as_chr_features.GRanges <- function(x, karyotype = NULL, chr_map = NULL,
                                    assembly = NULL, ...) {
  check_unused_args(...)
  info <- as.data.frame(GenomicRanges::seqinfo(x))
  known <- unique(stats::na.omit(info$genome))
  if (!is.null(assembly) && length(known) && any(known != assembly)) {
    stopf("`assembly` conflicts with the GRanges genome metadata.")
  }
  if (length(known) > 1L) stopf("GRanges contains multiple assemblies.")
  assembly <- assembly %||% if (length(known)) known else NULL
  result <- as_chr_features(as.data.frame(x),
    mapping = c(chr = "seqnames", start = "start", end = "end"),
    karyotype = karyotype, chr_map = chr_map, assembly = assembly)
  result
}

#' @rdname as_ideogram_data
#' @export
as_ideogram_data.Seqinfo <- function(x, ...) {
  info <- as.data.frame(x)
  karyotype <- read_karyotype(data.frame(Chr = rownames(info), End = info$seqlengths))
  karyotype$Assembly <- info$genome
  as_ideogram_data(karyotype, ...)
}

#' @rdname as_ideogram_data
#' @export
as_ideogram_data.GRanges <- function(x, ...) {
  as_ideogram_data(GenomicRanges::seqinfo(x), ...)
}
