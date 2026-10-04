# Explicit interval statistics for fixed chromosome windows.

bin_interval_windows <- function(data, karyotype, window, value, FUN, empty,
                                 method, group, count_position, ...) {
  input <- window_input(data, karyotype, window)
  k <- input$karyotype
  chr <- input$chr; start <- input$start; end <- input$end
  if (!is.null(FUN) || !is.null(empty)) stopf('Explicit methods define their own summary and empty values.')
  if (!is.null(group) && (!is.character(group) || anyNA(group) ||
      any(!group %in% names(data)) || anyDuplicated(group))) stopf('`group` must name distinct data columns.')
  reserved <- c('Chr', 'Start', 'End', 'Value', 'Width', 'N', 'N_valid', 'Rate')
  if (any(group %in% reserved)) stopf('Grouping columns conflict with window output columns.')
  options <- list(...)
  if (length(setdiff(names(options), 'na.rm')) || (length(options) && is.null(names(options)))) {
    stopf('Explicit window methods accept only `na.rm` in `...`.')
  }
  na.rm <- options$na.rm %||% FALSE
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) stopf('`na.rm` must be TRUE or FALSE.')
  if (method != 'count' && !'End' %in% names(data)) stopf('Interval statistics require an End column.')
  if (method == 'weighted_mean') {
    if (!is.character(value) || length(value) != 1L || !value %in% names(data) || !is.numeric(data[[value]])) {
      stopf('Weighted means require a numeric `value` column.')
    }
    values <- data[[value]]
    if (any(is.infinite(values))) stopf('Window values must be finite or missing.')
  } else {
    if (!is.null(value)) stopf('Count and coverage do not use a `value` column.')
    values <- rep(1, nrow(data))
  }
  if (method != 'count') count_position <- NULL else
    if (!'End' %in% names(data)) count_position <- 'start'
  at <- if (method == 'count') window_position(start, end, count_position) else NULL
  groups <- if (length(group)) unique(data[group]) else data.frame(.all = 1)
  parts <- list()
  for (gi in seq_len(nrow(groups))) {
    selected <- rep(TRUE, nrow(data))
    for (name in group) {
      v <- groups[[name]][gi]
      selected <- selected & if (is.na(v)) is.na(data[[name]]) else
        !is.na(data[[name]]) & data[[name]] == v
    }
    for (ci in seq_len(nrow(k))) {
      e <- window_edges(k$.end[ci] - k$.start[ci], window)
      lower <- e$lower + k$.start[ci] + 1
      upper <- e$upper + k$.start[ci]
      rows <- which(selected & chr == k$.chr[ci])
      result <- data.frame(Chr = k$.chr[ci], Start = lower, End = upper,
        Value = if (method == 'weighted_mean') NA_real_ else 0,
        Width = upper - lower + 1, N = 0L, N_valid = 0L, stringsAsFactors = FALSE)
      if (method == 'count') {
        index <- findInterval(at[rows], upper, left.open = TRUE) + 1L
        result$N <- tabulate(index, nbins = length(upper))
        result$N_valid <- result$N
        result$Value <- as.numeric(result$N)
        result$Rate <- result$Value / result$Width * 1e6
      } else result <- summarise_interval_windows(result, start[rows], end[rows],
        values[rows], method, na.rm)
      for (name in group) result[[name]] <- groups[[name]][gi]
      parts[[length(parts) + 1L]] <- result
    }
  }
  out <- if (length(parts)) do.call(rbind, parts) else data.frame(
    Chr = character(), Start = numeric(), End = numeric(), Value = numeric(),
    Width = numeric(), N = integer(), N_valid = integer())
  if (!length(parts)) {
    for (name in group) out[[name]] <- data[[name]][FALSE]
    if (method == 'count') out$Rate <- numeric()
  }
  rownames(out) <- NULL
  attr(out, 'window_summary') <- window_metadata(method, window, group, value,
    count_position, na.rm, empty = if (method == 'weighted_mean') NA_real_ else 0)
  out
}

summarise_interval_windows <- function(result, start, end, values, method, na.rm) {
  index <- order(start, end)
  start <- start[index]; end <- end[index]; values <- values[index]
  last <- findInterval(result$End, start)
  active <- integer(); previous <- 0L
  for (wi in seq_len(nrow(result))) {
    if (last[wi] > previous) active <- c(active, seq.int(previous + 1L, last[wi]))
    previous <- last[wi]
    active <- active[end[active] >= result$Start[wi]]
    if (!length(active)) next
    lo <- pmax(start[active], result$Start[wi])
    hi <- pmin(end[active], result$End[wi])
    valid <- !is.na(values[active])
    result$N[wi] <- length(active)
    result$N_valid[wi] <- sum(valid)
    if (method == 'coverage') {
      result$Value[wi] <- union_covered_bp(lo, hi) / result$Width[wi]
    } else if (any(valid) && (na.rm || all(valid))) {
      weights <- hi[valid] - lo[valid] + 1
      result$Value[wi] <- sum(values[active][valid] * (weights / sum(weights)))
    }
  }
  result
}

union_covered_bp <- function(start, end) {
  if (!length(start)) return(0)
  order <- order(start, end)
  start <- start[order]; end <- end[order]
  total <- 0; left <- start[1]; right <- end[1]
  if (length(start) > 1L) for (i in seq.int(2, length(start))) {
    if (start[i] <= right + 1) right <- max(right, end[i]) else {
      total <- total + right - left + 1
      left <- start[i]; right <- end[i]
    }
  }
  total + right - left + 1
}
