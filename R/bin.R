# Binning arbitrary genomic data into fixed windows.
#
# GFFex() does this for one specific case -- counting one feature type out of a
# GFF -- and every other kind of genomic data has to arrive pre-binned. This is
# the general form: any table with positions, any summary function.

#' Window edges for one chromosome
#'
#' Right edges are every whole multiple of `window` below `end`, then `end`
#' itself. Counting the multiples rather than calling `seq()` from `window` is
#' what lets a chromosome shorter than one window through -- `seq()` errors with
#' "wrong sign in 'by' argument" there. `unique()` keeps an exact multiple from
#' duplicating the last edge.
#'
#' @param end Chromosome length in base pairs.
#' @param window Window size in base pairs.
#' @return List of `lower` (0-based, exclusive) and `upper` (inclusive) edges.
#' @keywords internal
window_edges <- function(end, window) {
  upper <- unique(c(seq_len(floor(end / window)) * window, end))
  list(lower = c(0, upper[-length(upper)]), upper = upper)
}

#' Bin genomic data into fixed windows
#'
#' Summarises a table of positions or intervals into equal-width windows along
#' each chromosome, giving the `Chr` / `Start` / `End` / `Value` frame consumed
#' by chromosome track layers.
#'
#' A row is assigned to the window containing its midpoint when it has an `End`,
#' and its `Start` otherwise, so an interval is counted once however long it is.
#'
#' @param data Data frame with `Chr` and `Start`, and optionally `End`.
#' @param karyotype Data frame with `Chr` and `End`, or a path to a
#'   tab-separated karyotype file with a header. Chromosomes named here but
#'   absent from `data` come back filled with `empty` rather than disappearing.
#' @param window Window size in base pairs.
#' @param value Column of `data` to summarise. `NULL` counts rows.
#' @param FUN Summary function applied to each window's values. Defaults to
#'   [length()] when `value` is `NULL` and [mean()] otherwise.
#' @param empty Value for a window with no rows in it. Defaults to `0` for
#'   counts and `NA` for a summary of nothing.
#' @param ... Passed to `FUN`.
#' @return A data frame with `Chr`, `Start`, `End`, `Value`, in karyotype order.
#'   Windows are half-open on the left, matching [GFFex()].
#' @examples
#' kar <- data.frame(Chr = "A", Start = 0, End = 2500)
#' snps <- data.frame(Chr = "A", Start = c(10, 20, 1500), Qual = c(30, 50, 99))
#' bin_genome(snps, kar, window = 1000)
#' bin_genome(snps, kar, window = 1000, value = "Qual", FUN = max)
#' @export
bin_genome <- function(data, karyotype, window = 1e6, value = NULL, FUN = NULL,
                       empty = NULL, ...) {
  if (!is.data.frame(data)) {
    stopf("`data` must be a data frame, not %s.", class(data)[1])
  }
  miss <- setdiff(c("Chr", "Start"), names(data))
  if (length(miss)) {
    stopf("`data` is missing column%s %s.", if (length(miss) > 1) "s" else "",
          paste0("`", miss, "`", collapse = ", "))
  }
  if (!is.null(value) && !value %in% names(data)) {
    stopf("`value` names column `%s`, which is not in `data`.\n  Columns present: %s",
          value, paste0("`", names(data), "`", collapse = ", "))
  }
  if (!is.numeric(window) || length(window) != 1L || is.na(window) || window <= 0) {
    stopf("`window` must be a single positive number.")
  }
  kar <- if (is.data.frame(karyotype)) {
    karyotype
  } else {
    utils::read.table(karyotype, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  if (!all(c("Chr", "End") %in% names(kar))) {
    stopf("`karyotype` needs columns Chr and End.")
  }
  FUN <- FUN %||% if (is.null(value)) length else mean
  empty <- empty %||% if (is.null(value)) 0 else NA

  chrs <- as.character(kar$Chr)
  ends <- as.numeric(kar$End)
  if (anyNA(ends) || any(ends <= 0)) {
    stopf("`karyotype$End` must be positive and non-missing; bad for %s.",
          paste0("`", chrs[is.na(ends) | ends <= 0], "`", collapse = ", "))
  }

  seqname <- as.character(data$Chr)
  # Midpoint of the interval, so a feature is binned once wherever it is longest.
  at <- if ("End" %in% names(data)) {
    (as.numeric(data$Start) + as.numeric(data$End)) / 2
  } else {
    as.numeric(data$Start)
  }
  v <- if (is.null(value)) rep(1, nrow(data)) else data[[value]]

  parts <- lapply(seq_along(chrs), function(i) {
    e <- window_edges(ends[i], window)
    keep <- seqname == chrs[i] & !is.na(at) & at > 0 & at <= ends[i]
    out <- rep(empty, length(e$upper))
    if (any(keep)) {
      bin <- findInterval(at[keep], e$upper, left.open = TRUE) + 1L
      agg <- vapply(split(v[keep], factor(bin, levels = seq_along(e$upper))),
                    function(x) if (length(x)) as.numeric(FUN(x, ...)) else NA_real_,
                    numeric(1))
      # `empty` survives wherever the split produced nothing, so a count stays 0
      # and a mean-of-nothing stays NA rather than becoming a spurious number.
      out[!is.na(agg)] <- agg[!is.na(agg)]
    }
    data.frame(Chr = chrs[i], Start = e$lower + 1, End = e$upper, Value = out,
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
