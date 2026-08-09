# Semantic data for the next-generation renderer.
#
# This file deliberately contains no drawing coordinates.  It turns user
# columns into stable, prefixed semantic fields while preserving every source
# column for later aesthetic mappings.

#' Create semantic ideogram data
#'
#' `as_ideogram_data()` separates biological coordinates from rendering.  The
#' returned object contains the original karyotype columns plus canonical
#' `.chr`, `.start`, `.end`, `.centromere_start` and `.centromere_end` fields;
#' it never stores canvas pixels, device dimensions or pre-built grobs.
#'
#' @param x A karyotype-like object.
#' @param ... Passed to the class method.
#'
#' @return An `ideogram_data` object.
#' @export
as_ideogram_data <- function(x, ...) {
  UseMethod("as_ideogram_data")
}
#' @param mapping An aesthetic mapping containing `chr`, `start` and `end`.
#'   Expressions are evaluated in `x`, so source columns do not have to be
#'   renamed.
#' @param centromere Optional mapping containing `start` and `end`.  When it is
#'   `NULL`, columns `CE_start` and `CE_end` are used together when present.
#' @param cytoband Optional cytoband/interval table.
#' @param cytoband_mapping Mapping containing `chr`, `start` and `end` for
#'   `cytoband`.  `NULL` uses columns `Chr`, `Start` and `End`.
#' @rdname as_ideogram_data
#' @export
as_ideogram_data.data.frame <- function(
    x,
    mapping = NULL,
    centromere = NULL,
    cytoband = NULL,
    cytoband_mapping = NULL,
    ...) {
  check_unused_args(...)
  if (!nrow(x)) {
    stopf("`x` has no rows.")
  }

  mapping <- mapping %||% ggplot2::aes(chr = .data$Chr,
                                       start = .data$Start,
                                       end = .data$End)
  canonical <- mapped_fields(x, mapping, c("chr", "start", "end"),
                             what = "mapping")
  karyotype <- x
  karyotype$.chr <- validate_chr(canonical$chr, "mapping$chr")
  karyotype$.start <- validate_coordinate(canonical$start, "mapping$start")
  karyotype$.end <- validate_coordinate(canonical$end, "mapping$end")

  bad_interval <- karyotype$.start < 0 | karyotype$.end <= karyotype$.start
  if (any(bad_interval)) {
    stopf(paste0(
      "Karyotype coordinates are invalid for %s.\n",
      "  Require 0 <= start < end."),
      format_chr_rows(karyotype$.chr[bad_interval]))
  }
  if (anyDuplicated(karyotype$.chr)) {
    repeated <- unique(karyotype$.chr[duplicated(karyotype$.chr)])
    stopf("Chromosome identifiers must be unique; repeated: %s.",
          format_chr_rows(repeated))
  }

  auto_centromere <- is.null(centromere) &&
    any(c("CE_start", "CE_end") %in% names(x))
  if (auto_centromere) {
    if (!all(c("CE_start", "CE_end") %in% names(x))) {
      stopf("`x` has one of `CE_start`/`CE_end` but not both.")
    }
    centromere <- ggplot2::aes(start = .data$CE_start, end = .data$CE_end)
  }

  if (is.null(centromere)) {
    karyotype$.centromere_start <- NA_real_
    karyotype$.centromere_end <- NA_real_
  } else {
    ce <- mapped_fields(x, centromere, c("start", "end"),
                        what = "centromere")
    karyotype$.centromere_start <- validate_coordinate(
      ce$start, "centromere$start", allow_na = TRUE)
    karyotype$.centromere_end <- validate_coordinate(
      ce$end, "centromere$end", allow_na = TRUE)
    one_missing <- xor(is.na(karyotype$.centromere_start),
                       is.na(karyotype$.centromere_end))
    if (any(one_missing)) {
      stopf("Centromere start and end must either both be present or both be missing for %s.",
            format_chr_rows(karyotype$.chr[one_missing]))
    }
    has_ce <- !is.na(karyotype$.centromere_start)
    bad_ce <- has_ce & (
      karyotype$.centromere_start < karyotype$.start |
      karyotype$.centromere_end < karyotype$.centromere_start |
      karyotype$.centromere_end > karyotype$.end
    )
    if (any(bad_ce)) {
      stopf(paste0(
        "Centromere coordinates are out of range for %s.\n",
        "  Require chromosome start <= centromere start <= ",
        "centromere end <= chromosome end."),
        format_chr_rows(karyotype$.chr[bad_ce]))
    }
  }

  bands <- canonical_cytoband(cytoband, cytoband_mapping, karyotype)
  structure(
    list(
      karyotype = karyotype,
      cytoband = bands,
      mapping = mapping,
      centromere_mapping = centromere,
      cytoband_mapping = cytoband_mapping
    ),
    class = "ideogram_data"
  )
}

#' @rdname as_ideogram_data
#' @export
as_ideogram_data.ideogram_data <- function(x, ...) {
  check_unused_args(...)
  x
}

#' @export
print.ideogram_data <- function(x, ...) {
  k <- x$karyotype
  cat("<ideogram_data>\n")
  cat("  chromosomes: ", nrow(k), "\n", sep = "")
  cat("  range      : ", format(min(k$.start)), " .. ",
      format(max(k$.end)), " bp\n", sep = "")
  cat("  centromeres: ", sum(!is.na(k$.centromere_start)), "\n", sep = "")
  cat("  cytobands  : ", if (is.null(x$cytoband)) 0L else nrow(x$cytoband),
      "\n", sep = "")
  invisible(x)
}

mapped_fields <- function(data, mapping, required, what) {
  if (!inherits(mapping, "uneval") && !inherits(mapping, "ggplot2::mapping")) {
    stopf("`%s` must be created with `ggplot2::aes()`.", what)
  }
  missing <- setdiff(required, names(mapping))
  if (length(missing)) {
    stopf("`%s` is missing aesthetic%s %s.", what,
          if (length(missing) > 1L) "s" else "",
          paste0("`", missing, "`", collapse = ", "))
  }

  values <- lapply(required, function(field) {
    value <- tryCatch(
      rlang::eval_tidy(mapping[[field]], data = data),
      error = function(error) {
        stopf("Could not evaluate `%s$%s`: %s", what, field,
              conditionMessage(error))
      }
    )
    recycle_semantic(value, nrow(data), sprintf("%s$%s", what, field))
  })
  names(values) <- required
  values
}

recycle_semantic <- function(value, n, what) {
  if (length(value) == n) return(value)
  if (length(value) == 1L && n > 0L) return(rep(value, n))
  if (length(value) == 0L && n == 0L) return(value)
  stopf("`%s` evaluated to length %d; expected length %d or 1.",
        what, length(value), n)
}

validate_chr <- function(x, what) {
  x <- as.character(x)
  if (anyNA(x) || any(!nzchar(x))) {
    stopf("`%s` must contain non-missing, non-empty chromosome identifiers.",
          what)
  }
  x
}

validate_coordinate <- function(x, what, allow_na = FALSE) {
  if (!is.numeric(x)) {
    stopf("`%s` must be numeric.", what)
  }
  bad <- is.infinite(x) | (!allow_na & is.na(x))
  if (any(bad)) {
    stopf("`%s` must contain %sfinite coordinates.", what,
          if (allow_na) "only missing or " else "")
  }
  as.numeric(x)
}

canonical_cytoband <- function(cytoband, mapping, karyotype) {
  if (is.null(cytoband)) return(NULL)
  if (!is.data.frame(cytoband)) {
    stopf("`cytoband` must be a data frame, not %s.", class(cytoband)[1])
  }
  mapping <- mapping %||% ggplot2::aes(chr = .data$Chr,
                                       start = .data$Start,
                                       end = .data$End)
  fields <- mapped_fields(cytoband, mapping, c("chr", "start", "end"),
                          what = "cytoband_mapping")
  bands <- cytoband
  bands$.chr <- validate_chr(fields$chr, "cytoband_mapping$chr")
  bands$.start <- validate_coordinate(fields$start, "cytoband_mapping$start")
  bands$.end <- validate_coordinate(fields$end, "cytoband_mapping$end")

  unknown <- setdiff(unique(bands$.chr), karyotype$.chr)
  if (length(unknown)) {
    stopf("`cytoband` contains unknown chromosome%s %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  index <- match(bands$.chr, karyotype$.chr)
  bad <- bands$.start < karyotype$.start[index] |
    bands$.end <= bands$.start |
    bands$.end > karyotype$.end[index]
  if (any(bad)) {
    stopf(paste0(
      "Cytoband coordinates are invalid for %s.\n",
      "  Require chromosome start <= band start < band end <= chromosome end."),
      format_chr_rows(unique(bands$.chr[bad])))
  }
  bands
}

format_chr_rows <- function(x) {
  paste0("`", unique(as.character(x)), "`", collapse = ", ")
}

check_unused_args <- function(...) {
  dots <- list(...)
  if (length(dots)) {
    stopf("Unused argument%s: %s.",
          if (length(dots) > 1L) "s" else "",
          paste0("`", names(dots) %||% rep("", length(dots)), "`",
                 collapse = ", "))
  }
  invisible(NULL)
}
