# Fixed chromosome windows shared by genomic summaries and GFF feature counts.

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
#' Count and `FUN` summaries assign each row to one window. By default the
#' representative position is `(Start + End) / 2` when `End` is present, and
#' `Start` otherwise. A position exactly at a window's end belongs to that
#' window. Coverage and weighted means use complete interval overlaps.
#'
#' @param data Data frame with `Chr` and `Start`, and optionally `End`.
#'   Coordinates are positive integers using 1-based closed intervals. Unknown
#'   chromosomes, missing coordinates and intervals outside the karyotype are errors.
#' @param karyotype Data frame with `Chr` and `End`, or a path to a
#'   tab-separated karyotype file with a header. Chromosomes named here but
#'   absent from `data` come back filled with `empty` rather than disappearing.
#'   Optional `Start` gives the 0-based left boundary; it defaults to zero.
#' @param window Positive integer window size in base pairs. Windows begin at
#'   each chromosome's `Start + 1` and end at its `End`.
#' @param value Column of `data` to summarise. `NULL` counts rows.
#' @param FUN Summary function applied to each window's values. Defaults to
#'   [length()] when `value` is `NULL` and [mean()] otherwise.
#' @param empty Value for a window with no rows in it. Defaults to `0` for
#'   counts and `NA` for a summary of nothing. Missing results from occupied
#'   windows remain missing; they are not replaced by `empty`.
#' @param ... Passed to `FUN`.
#' @param method Explicit `"count"`, `"coverage"`, or `"weighted_mean"`.
#'   Count assigns each feature once by `count_position`. Coverage is the union of
#'   covered bases divided by actual window width. Weighted mean uses each
#'   interval's overlap length; overlapping observations contribute separately.
#'   `NULL` retains the original `FUN` behavior.
#' @param group Optional column names for independent grouped summaries with an
#'   explicit `method`. Every observed group receives every karyotype window.
#'   For weighted means, `na.rm = TRUE` in `...` excludes missing values and
#'   their weights; otherwise a missing overlapping value produces `NA`.
#' @param count_position Representative position for count and `FUN` summaries:
#'   `"midpoint"` (default), `"start"`, or `"end"`. With no `End` column all
#'   three use `Start`. Coverage and weighted means do not accept an explicit
#'   setting because they use the complete interval.
#' @return A data frame with `Chr`, `Start`, `End`, `Value`, in karyotype order.
#'   Window coordinates are 1-based closed intervals.
#'   Explicit methods additionally return `Width`, `N` and `N_valid`; counts
#'   include `Rate` per Mb. All results carry a `window_summary` attribute with
#'   `method`, `group`, `coordinates`, `unit`, `value`, `window`, `count_position`,
#'   `na.rm`, `overlap`, `empty` and `summary`. Default row counts use method
#'   `"count"`; other `FUN` calls use `"summary"` and record the function
#'   expression. Count units are `"features"`, coverage units are `"fraction"`.
#'   For weighted means and `FUN` summaries, `unit` is `NULL` and `value` names
#'   the input column; units are not inferred from column names. `window` is
#'   the requested width in bp; actual widths follow each row's boundaries.
#'   For `FUN`, `na.rm` records only an explicitly forwarded argument.
#' @examples
#' kar <- data.frame(Chr = "A", Start = 0, End = 2500)
#' snps <- data.frame(Chr = "A", Start = c(10, 20, 1500), Qual = c(30, 50, 99))
#' bin_genome(snps, kar, window = 1000)
#' bin_genome(snps, kar, window = 1000, value = "Qual", FUN = max)
#' intervals <- data.frame(Chr = "A", Start = 950, End = 1150)
#' counts <- bin_genome(intervals, kar, window = 1000, method = "count",
#'   count_position = "start")
#' attr(counts, "window_summary")
#' @export
bin_genome <- function(data, karyotype, window = 1e6, value = NULL, FUN = NULL,
                       empty = NULL, ..., method = NULL, group = NULL,
                       count_position = c("midpoint", "start", "end")) {
  position_supplied <- !missing(count_position)
  count_position <- match.arg(count_position)
  if (!is.null(method)) {
    method <- match.arg(method, c('count', 'coverage', 'weighted_mean'))
    if (method != 'count' && position_supplied) {
      stopf('`count_position` applies only to count or FUN summaries.')
    }
    return(bin_interval_windows(data, karyotype, window, value, FUN, empty,
      method, group, count_position, ...))
  }
  if (!is.null(group)) stopf('Grouped windows require an explicit `method`.')
  input <- window_input(data, karyotype, window)
  k <- input$karyotype
  if (!is.null(value) && !value %in% names(data)) {
    stopf("`value` names column `%s`, which is not in `data`.\n  Columns present: %s",
          value, paste0("`", names(data), "`", collapse = ", "))
  }
  method <- if (is.null(value) && is.null(FUN)) "count" else "summary"
  summary <- if (method == "count") NULL else if (is.null(FUN)) "mean" else
    paste(deparse(substitute(FUN)), collapse = " ")
  FUN <- FUN %||% if (is.null(value)) length else mean
  empty <- empty %||% if (is.null(value)) 0 else NA
  if (!"End" %in% names(data)) count_position <- "start"
  at <- window_position(input$start, input$end, count_position)
  v <- if (is.null(value)) rep(1, nrow(data)) else data[[value]]

  parts <- lapply(seq_len(nrow(k)), function(i) {
    e <- window_edges(k$.end[i] - k$.start[i], window)
    e$lower <- e$lower + k$.start[i]
    e$upper <- e$upper + k$.start[i]
    keep <- input$chr == k$.chr[i]
    out <- rep(empty, length(e$upper))
    if (any(keep)) {
      bin <- findInterval(at[keep], e$upper, left.open = TRUE) + 1L
      groups <- split(v[keep], factor(bin, levels = seq_along(e$upper)))
      occupied <- lengths(groups) > 0L
      out[occupied] <- vapply(groups[occupied],
        function(x) as.numeric(FUN(x, ...)), numeric(1))
    }
    data.frame(Chr = k$.chr[i], Start = e$lower + 1, End = e$upper, Value = out,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  na_argument <- match("na.rm", ...names())
  attr(out, "window_summary") <- window_metadata(method, window,
    value = value, count_position = count_position,
    na.rm = if (is.na(na_argument)) NULL else ...elt(na_argument),
    empty = empty, summary = summary)
  out
}

window_position <- function(start, end, position) {
  switch(position, midpoint = (start + end) / 2, start = start, end = end)
}

window_metadata <- function(method, window, group = NULL, value = NULL,
    count_position = NULL, na.rm = NULL, empty = NULL, summary = NULL) {
  list(method = method, group = group, coordinates = "1-based closed",
    unit = switch(method, count = "features", coverage = "fraction", NULL),
    value = value, window = window, count_position = count_position,
    na.rm = na.rm, overlap = switch(method, coverage = "union",
      weighted_mean = "independent observations", count_position),
    empty = empty, summary = summary)
}

window_input <- function(data, karyotype, window) {
  if (!is.data.frame(data)) {
    stopf("`data` must be a data frame, not %s.", class(data)[1])
  }
  miss <- setdiff(c("Chr", "Start"), names(data))
  if (length(miss)) {
    stopf("`data` is missing column%s %s.", if (length(miss) > 1) "s" else "",
          paste0("`", miss, "`", collapse = ", "))
  }
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window <= 0 || window != floor(window)) {
    stopf("`window` must be a positive integer.")
  }
  if (is.character(karyotype)) {
    karyotype <- utils::read.table(karyotype, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  }
  if (is.data.frame(karyotype) && !all(c("Chr", "End") %in% names(karyotype))) {
    stopf("`karyotype` needs columns Chr and End.")
  }
  if (is.data.frame(karyotype) && !"Start" %in% names(karyotype)) karyotype$Start <- 0
  k <- as_ideogram_data(karyotype)$karyotype
  if (any(k$.start != floor(k$.start) | k$.end != floor(k$.end))) {
    stopf("Window karyotypes require integer boundaries.")
  }
  chr <- validate_chr(data$Chr, "Chr")
  start <- validate_coordinate(data$Start, "Start")
  end <- if ("End" %in% names(data)) validate_coordinate(data$End, "End") else start
  i <- match(chr, k$.chr)
  if (anyNA(i)) stopf("Window data contain unknown chromosomes.")
  if (any(start != floor(start) | end != floor(end) | start > end |
          start <= k$.start[i] | end > k$.end[i])) {
    stopf("Window intervals require closed integer coordinates inside the source chromosome bounds.")
  }
  list(karyotype = k, chr = chr, start = start, end = end)
}
