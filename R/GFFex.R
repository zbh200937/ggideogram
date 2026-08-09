# Feature density from a GFF, in the shape `overlaid` wants.
#
# Same output as RIdeogram::GFFex(), but in base R (the original pulled in tidyr
# only to split "(a,b]" interval labels back into numbers, which is a round trip
# through text we can skip) and without two of its failure modes: a chromosome
# whose length is an exact multiple of `window` produced duplicated breaks and
# made cut() error, and a chromosome named in the karyotype but absent from the
# GFF silently vanished from the result instead of appearing with zero counts.

#' Count features per window from a GFF file
#'
#' Produces a `Chr` / `Start` / `End` / `Value` data frame suitable for
#' [geom_track_line()], [geom_track_col()] and other chromosome tracks.
#'
#' @param input Path to a GFF file, or a data frame of one already read (columns
#'   are used positionally, as in GFF: 1 = sequence name, 3 = feature type,
#'   4 = start).
#' @param karyotype Path to a tab-separated karyotype file with a header, or a
#'   data frame with `Chr` and `End`.
#' @param feature Feature type to count, matched against column 3 of the GFF.
#' @param window Window size in base pairs.
#' @return A data frame with columns `Chr`, `Start`, `End`, `Value`, in
#'   karyotype order. Windows are half-open on the left: a feature at position
#'   `p` falls in the window with `Start <= p <= End`. The last window of each
#'   chromosome is truncated to the chromosome's `End`.
#' @export
GFFex <- function(input, karyotype, feature = "gene", window = 1000000) {
  gff <- if (is.data.frame(input)) {
    input
  } else {
    utils::read.table(input, stringsAsFactors = FALSE, header = FALSE,
                      comment.char = "#", sep = "\t", quote = "")
  }
  if (ncol(gff) < 4) stopf("`input` needs at least 4 GFF columns, got %d.", ncol(gff))

  kar <- if (is.data.frame(karyotype)) {
    karyotype
  } else {
    utils::read.table(karyotype, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  if (!all(c("Chr", "End") %in% names(kar))) {
    stopf("`karyotype` needs columns Chr and End.")
  }
  if (!is.numeric(window) || length(window) != 1L || window <= 0) {
    stopf("`window` must be a single positive number.")
  }

  chrs <- as.character(kar$Chr)
  seqname <- as.character(gff[[1]])
  keep <- seqname %in% chrs & as.character(gff[[3]]) == feature
  seqname <- seqname[keep]
  pos <- as.numeric(gff[[4]][keep])

  ends <- as.numeric(kar$End)
  if (anyNA(ends) || any(ends <= 0)) {
    stopf("`karyotype$End` must be positive and non-missing; bad for %s.",
          paste0("`", chrs[is.na(ends) | ends <= 0], "`", collapse = ", "))
  }

  parts <- lapply(seq_along(chrs), function(i) {
    chr <- chrs[i]
    end <- ends[i]
    e <- window_edges(end, window)
    upper <- e$upper
    lower <- e$lower
    p <- pos[seqname == chr]
    p <- p[!is.na(p) & p > 0 & p <= end]
    counts <- if (length(p)) {
      tabulate(findInterval(p, upper, left.open = TRUE) + 1L, nbins = length(upper))
    } else {
      integer(length(upper))
    }
    data.frame(Chr = chr, Start = lower + 1, End = upper, Value = counts,
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
