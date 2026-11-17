#' Read genomic annotations in BED, bedGraph, GFF3 or GTF format
#'
#' BED3--6 and bedGraph starts are converted from zero-based, half-open
#' coordinates to one-based, closed coordinates (`Start = chromStart + 1`,
#' `End = chromEnd`). GFF3 and GTF already use one-based, closed coordinates.
#' Zero-width BED insertion sites need an explicit point representation and are
#' rejected. Original GFF/GTF attributes are retained in `Attributes`; common
#' identifiers are also extracted without collapsing rows or parent
#' relationships. Use [resolve_gene_models()] to expand GFF3/GTF parent
#' hierarchies for [geom_genemodel()].
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
#' @return A data frame with `Chr`, `Start`, `End` and annotation columns
#'   (`Value` and `Mid` for `bedgraph`). Attributes `source_coordinates` and
#'   `assembly` record input conventions.
#' @seealso [resolve_gene_models()], [read_chr_links()], [as_chr_features()], [read_karyotype()]
#' @export
read_chr_features <- function(file, format = c("auto", "bed", "bedgraph", "gff3", "gtf"),
                              karyotype = NULL, chr_map = NULL, assembly = NULL) {
  format <- match.arg(format)
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stopf("`file` must be one local annotation file path.")
  }
  if (format == "auto") {
    extension <- tolower(tools::file_ext(sub("\\.gz$", "", file, ignore.case = TRUE)))
    format <- switch(extension,
      bed = "bed", bedgraph = "bedgraph", bg = "bedgraph",
      gff = "gff3", gff3 = "gff3", gtf = "gtf",
      stopf("Cannot infer annotation format; supply `format`."))
  }
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  fasta <- match("##FASTA", trimws(lines))
  if (!is.na(fasta)) lines <- utils::head(lines, fasta - 1L)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  if (format %in% c("bed", "bedgraph")) {
    is_bed_coord <- grepl("^[^\t]+\t[0-9]+\t[0-9]+(\t|$)", lines)
    lines <- lines[is_bed_coord | !grepl("^(track|browser)($|[ \t])", lines)]
  }
  if (!length(lines)) stopf("The annotation file contains no features.")
  input <- utils::read.table(text = lines, sep = "\t", header = FALSE,
    quote = "", comment.char = "", colClasses = "character", na.strings = NULL)
  if (format == "bed") {
    if (!ncol(input) %in% 3:6) stopf("BED input must have 3 to 6 columns.")
    names(input) <- c("Chr", "Start", "End", "Name", "Score", "Strand")[seq_len(ncol(input))]
    system <- "0-based-half-open"
  } else if (format == "bedgraph") {
    if (ncol(input) != 4L) stopf("bedGraph input must have 4 columns (chrom, start, end, value).")
    names(input) <- c("Chr", "Start", "End", "Value")
    raw_val <- suppressWarnings(as.numeric(input$Value))
    if (anyNA(raw_val) || any(!is.finite(raw_val))) {
      stopf("bedGraph `Value` column must contain finite numeric values.")
    }
    input$Value <- raw_val
    system <- "0-based-half-open"
  } else {
    if (ncol(input) != 9L) stopf("GFF3/GTF input must have 9 columns.")
    names(input) <- c("Chr", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes")
    keys <- if (format == "gff3") c("ID", "Parent", "Name") else c("gene_id", "transcript_id", "gene_name")
    for (key in keys) input[[key]] <- feature_attribute(input$Attributes, key, format)
    system <- "1-based-closed"
  }
  if ("Strand" %in% names(input) && any(!input$Strand %in% c("+", "-", ".", "*", "?"))) {
    stopf("Feature `Strand` must be `+`, `-`, `.`, `*`, or `?`.")
  }
  input$Start <- suppressWarnings(as.numeric(input$Start))
  input$End <- suppressWarnings(as.numeric(input$End))
  result <- as_chr_features(input, coordinate_system = system, karyotype = karyotype,
                            chr_map = chr_map, assembly = assembly)
  if (format == "bedgraph") result$Mid <- (result$Start - 1 + result$End) / 2
  result
}

decode_gff3_percent <- function(x) {
  vapply(x, function(s) {
    m <- gregexpr("(%[0-9A-Fa-f]{2})+", s, perl = TRUE)
    regmatches(s, m) <- lapply(regmatches(s, m), function(chunks) {
      vapply(chunks, function(chunk) {
        bytes <- strtoi(substring(chunk, seq(2L, nchar(chunk), 3L),
                                  seq(3L, nchar(chunk), 3L)), 16L)
        bytes <- bytes[bytes != 0L]
        if (!length(bytes)) "" else enc2utf8(rawToChar(as.raw(bytes)))
      }, character(1))
    })
    s
  }, character(1), USE.NAMES = FALSE)
}

split_gtf_attributes <- function(value) {
  chars <- strsplit(value, "", fixed = TRUE)[[1]]
  if (!length(chars)) return(character())
  parts <- character()
  buf <- character()
  in_quote <- FALSE
  escaped <- FALSE
  for (ch in chars) {
    if (escaped) {
      buf <- c(buf, ch)
      escaped <- FALSE
    } else if (in_quote && ch == "\\") {
      buf <- c(buf, ch)
      escaped <- TRUE
    } else if (ch == "\"") {
      in_quote <- !in_quote
      buf <- c(buf, ch)
    } else if (!in_quote && ch == ";") {
      token <- trimws(paste(buf, collapse = ""))
      if (nzchar(token)) parts <- c(parts, token)
      buf <- character()
    } else {
      buf <- c(buf, ch)
    }
  }
  token <- trimws(paste(buf, collapse = ""))
  if (nzchar(token)) parts <- c(parts, token)
  parts
}

unquote_gtf_value <- function(x) {
  vapply(x, function(s) {
    s <- trimws(s)
    if (nchar(s) >= 2L && startsWith(s, "\"") && endsWith(s, "\"")) {
      s <- substr(s, 2L, nchar(s) - 1L)
      s <- gsub("\\\"", "\"", s, fixed = TRUE)
      s <- gsub("\\\\", "\\", s, fixed = TRUE)
    }
    s
  }, character(1), USE.NAMES = FALSE)
}

feature_attribute <- function(attributes, key, format) {
  vapply(attributes, function(value) {
    fields <- if (format == "gtf") {
      split_gtf_attributes(value)
    } else trimws(strsplit(value, ";", fixed = TRUE)[[1]])
    pattern <- if (format == "gff3") paste0("^", key, "=") else paste0("^", key, "\\s+")
    matches <- fields[grepl(pattern, fields)]
    if (!length(matches)) return(NA_character_)
    result <- sub(pattern, "", matches)
    result <- if (format == "gff3") decode_gff3_percent(result) else unquote_gtf_value(result)
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
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stopf("Install the optional Bioconductor package `GenomicRanges` to use GRanges inputs.")
  }
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
  if (!requireNamespace("Seqinfo", quietly = TRUE)) {
    stopf("Install the optional Bioconductor package `Seqinfo` to use Seqinfo inputs.")
  }
  info <- as.data.frame(x)
  karyotype <- read_karyotype(data.frame(Chr = rownames(info), End = info$seqlengths))
  karyotype$Assembly <- info$genome
  as_ideogram_data(karyotype, ...)
}

#' @rdname as_ideogram_data
#' @export
as_ideogram_data.GRanges <- function(x, ...) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stopf("Install the optional Bioconductor package `GenomicRanges` to use GRanges inputs.")
  }
  as_ideogram_data(GenomicRanges::seqinfo(x), ...)
}

#' Resolve GFF3 or GTF feature hierarchies for gene model tracks
#'
#' Expands multi-parent GFF3 `Parent` links (such as exons shared across
#' isoforms) and links subfeatures (`exon`, `CDS`, UTRs) to their transcript and
#' gene identifiers. Supports both three-level hierarchies (`gene -> mRNA ->
#' exon/CDS`), non-coding transcript classes (`lnc_RNA`, `rRNA`, `tRNA`, etc.),
#' two-level GFF3 hierarchies (`gene -> exon/CDS` without intermediate transcript
#' rows), and GTF tables (`gene_id` / `transcript_id`). All source columns are
#' preserved.
#'
#' @param features A data frame returned by [read_chr_features()] for a GFF3 or
#'   GTF file, or a compatible table containing `Chr`, `Start`, `End`, `Type`
#'   and either `ID`/`Parent` (GFF3) or `gene_id`/`transcript_id` (GTF).
#' @param transcript_types Feature `Type` values recognized as transcripts when
#'   resolving GFF3 hierarchies.
#' @param gene_types Feature `Type` values recognized as top-level gene
#'   containers.
#' @param strip_prefix Logical; if `TRUE` (default), strips common GFF3 type
#'   prefixes such as `"gene:"` and `"transcript:"` from resolved `Gene` and
#'   `Transcript` identifiers.
#'
#' @return A data frame ready for [geom_genemodel()], retaining original columns
#'   plus `Transcript`, `Gene`, `Symbol` and `Part` (`"CDS"`, `"UTR"` or `NA`).
#' @seealso [read_chr_features()], [geom_genemodel()]
#' @export
resolve_gene_models <- function(
    features,
    transcript_types = c("mRNA", "transcript", "lnc_RNA", "ncRNA", "rRNA",
      "tRNA", "snRNA", "snoRNA", "miRNA", "primary_transcript",
      "pseudogenic_transcript", "V_gene_segment", "C_gene_segment"),
    gene_types = c("gene", "ncRNA_gene", "pseudogene", "transposable_element_gene"),
    strip_prefix = TRUE) {
  if (!is.data.frame(features)) {
    stopf("`features` must be a data frame from `read_chr_features()`.")
  }
  required <- c("Chr", "Start", "End", "Type")
  missing_cols <- setdiff(required, names(features))
  if (length(missing_cols)) {
    stopf("`features` is missing required column(s): %s.",
          paste(missing_cols, collapse = ", "))
  }
  if (!is.logical(strip_prefix) || length(strip_prefix) != 1L || is.na(strip_prefix)) {
    stopf("`strip_prefix` must be TRUE or FALSE.")
  }
  if (!"Strand" %in% names(features)) {
    features$Strand <- rep(".", nrow(features))
  }

  has_gff3 <- all(c("ID", "Parent") %in% names(features)) &&
    (any(!is.na(features$ID)) || any(!is.na(features$Parent)) || !any(c("gene_id", "transcript_id") %in% names(features)))
  has_gtf <- any(c("gene_id", "transcript_id") %in% names(features))
  if (!has_gff3 && !has_gtf && "Attributes" %in% names(features)) {
    probe <- stats::na.omit(features$Attributes)
    if (length(probe) && any(grepl("(^|;)\\s*(ID|Parent)=", probe))) {
      for (k in c("ID", "Parent", "Name")) features[[k]] <- feature_attribute(features$Attributes, k, "gff3")
      has_gff3 <- TRUE
    } else if (length(probe) && any(grepl("(^|;)\\s*(gene_id|transcript_id)\\s+", probe))) {
      for (k in c("gene_id", "transcript_id", "gene_name")) features[[k]] <- feature_attribute(features$Attributes, k, "gtf")
      has_gtf <- TRUE
    }
  }
  if (!has_gff3 && !has_gtf) {
    stopf("`features` must contain GFF3 (`ID`, `Parent`) or GTF (`gene_id`, `transcript_id`) columns.")
  }

  clean_id <- function(x) {
    x <- as.character(x)
    if (isTRUE(strip_prefix)) {
      x <- sub("^(gene|transcript|mRNA|RNA|ncRNA|CDS|exon):", "", x)
    }
    x
  }
  classify_part <- function(type_vec) {
    low <- tolower(as.character(type_vec))
    out <- rep(NA_character_, length(low))
    out[low == "cds"] <- "CDS"
    out[low %in% c("utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr")] <- "UTR"
    out
  }

  if (!nrow(features)) {
    out <- features
    out$Transcript <- character(0)
    out$Gene <- character(0)
    out$Symbol <- character(0)
    out$Part <- character(0)
    return(out)
  }

  split_csv <- function(x) {
    lapply(ifelse(is.na(x), "", as.character(x)), function(s) {
      if (!nzchar(s)) return(character(0))
      toks <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
      toks[nzchar(toks)]
    })
  }

  if (has_gff3 && !(has_gtf && all(is.na(features$ID)) && all(is.na(features$Parent)))) {
    type_str <- as.character(features$Type)
    is_tx <- !is.na(type_str) & type_str %in% transcript_types &
      !is.na(features$ID) & nzchar(features$ID)
    is_gene <- !is.na(type_str) & type_str %in% gene_types &
      !is.na(features$ID) & nzchar(features$ID)

    tx_rows <- features[is_tx, , drop = FALSE]
    gene_rows <- features[is_gene, , drop = FALSE]
    tx_ids <- tx_rows$ID
    gene_ids <- gene_rows$ID

    name_col <- if ("Name" %in% names(features)) as.character(features$Name) else rep(NA_character_, nrow(features))
    id_col <- as.character(features$ID)
    valid_named <- !is.na(id_col) & nzchar(id_col) & !is.na(name_col) & nzchar(name_col)
    name_lookup <- stats::setNames(name_col[valid_named], id_col[valid_named])

    raw_links <- ifelse(is_tx, features$ID, ifelse(is_gene, NA_character_, features$Parent))
    tx_tokens <- split_csv(raw_links)
    models <- features[rep(seq_len(nrow(features)), lengths(tx_tokens)), , drop = FALSE]
    if (!nrow(models)) {
      models$Transcript <- character(0)
      models$Gene <- character(0)
      models$Symbol <- character(0)
      models$Part <- character(0)
      rownames(models) <- NULL
      return(models)
    }
    models$.raw_tx <- unlist(tx_tokens, use.names = FALSE)
    in_tx <- models$.raw_tx %in% tx_ids
    in_direct_gene <- !in_tx & (models$.raw_tx %in% gene_ids)
    models <- models[in_tx | in_direct_gene, , drop = FALSE]
    if (!nrow(models)) {
      models$.raw_tx <- NULL
      models$Transcript <- character(0)
      models$Gene <- character(0)
      models$Symbol <- character(0)
      models$Part <- character(0)
      rownames(models) <- NULL
      return(models)
    }

    tx_parent_raw <- ifelse(
      models$.raw_tx %in% tx_ids,
      tx_rows$Parent[match(models$.raw_tx, tx_ids)],
      models$.raw_tx
    )
    tx_parent_raw[is.na(tx_parent_raw) | !nzchar(tx_parent_raw)] <-
      models$.raw_tx[is.na(tx_parent_raw) | !nzchar(tx_parent_raw)]
    gene_tokens <- split_csv(tx_parent_raw)
    models <- models[rep(seq_len(nrow(models)), lengths(gene_tokens)), , drop = FALSE]
    raw_gene <- unlist(gene_tokens, use.names = FALSE)
    raw_tx <- models$.raw_tx
    models$.raw_tx <- NULL

    models$Transcript <- clean_id(raw_tx)
    models$Gene <- clean_id(raw_gene)
    sym <- unname(name_lookup[raw_gene])
    tx_sym <- unname(name_lookup[raw_tx])
    sym[is.na(sym) | !nzchar(sym)] <- tx_sym[is.na(sym) | !nzchar(sym)]
    sym[is.na(sym) | !nzchar(sym)] <- models$Gene[is.na(sym) | !nzchar(sym)]
    models$Symbol <- sym
    models$Part <- classify_part(models$Type)
    rownames(models) <- NULL
    return(models)
  }

  # GTF branch
  type_str <- as.character(features$Type)
  keep <- is.na(type_str) | !(type_str %in% gene_types)
  models <- features[keep, , drop = FALSE]
  gid <- if ("gene_id" %in% names(models)) as.character(models$gene_id) else rep(NA_character_, nrow(models))
  tid <- if ("transcript_id" %in% names(models)) as.character(models$transcript_id) else rep(NA_character_, nrow(models))
  tid[is.na(tid) | !nzchar(tid)] <- gid[is.na(tid) | !nzchar(tid)]
  gid[is.na(gid) | !nzchar(gid)] <- tid[is.na(gid) | !nzchar(gid)]
  valid <- !is.na(tid) & nzchar(tid) & !is.na(gid) & nzchar(gid)
  models <- models[valid, , drop = FALSE]
  gid <- gid[valid]
  tid <- tid[valid]
  models$Transcript <- clean_id(tid)
  models$Gene <- clean_id(gid)
  gname <- if ("gene_name" %in% names(models)) as.character(models$gene_name) else rep(NA_character_, nrow(models))
  gname[is.na(gname) | !nzchar(gname)] <- models$Gene[is.na(gname) | !nzchar(gname)]
  models$Symbol <- gname
  models$Part <- classify_part(models$Type)
  rownames(models) <- NULL
  models
}
