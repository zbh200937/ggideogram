# Explicit interval statistics for fixed chromosome windows.

encode_group_token <- function(col) {
  out <- as.character(col)
  if (is.numeric(col)) {
    is_nan <- is.nan(col)
    is_na <- is.na(col) & !is_nan
    out[is_nan] <- "\x01NaN"
    out[is_na] <- "\x01NA"
  } else {
    out[is.na(col)] <- "\x01NA"
  }
  out
}

bin_interval_windows <- function(data, karyotype, window, value, FUN, empty,
                                 method, group, count_position, ...) {
  input <- window_input(data, karyotype, window)
  k <- input$karyotype
  chr <- input$chr; start <- input$start; end <- input$end
  if (!is.null(FUN) || !is.null(empty)) stopf('Explicit methods define their own summary and empty values.')
  if (!is.null(group) && (!is.character(group) || anyNA(group) ||
      any(!group %in% names(data)) || anyDuplicated(group))) stopf('`group` must name distinct data columns.')
  if (length(group) && any(vapply(data[group], function(col) is.list(col) || !is.null(dim(col)), logical(1)))) {
    stopf('Grouping columns must be atomic vectors.')
  }
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

  if (length(group)) {
    g_tokens <- lapply(group, function(name) encode_group_token(data[[name]]))
    row_gkey <- do.call(paste, c(g_tokens, list(sep = "\x02")))
    uniq_gkey <- unique(row_gkey)
    first_g <- match(uniq_gkey, row_gkey)
    groups <- data[first_g, group, drop = FALSE]
    rownames(groups) <- NULL
    buckets <- split(seq_len(nrow(data)), paste(row_gkey, chr, sep = "\r"))
  } else {
    uniq_gkey <- ".all"
    groups <- data.frame(.all = 1)
    buckets <- split(seq_len(nrow(data)), paste(".all", chr, sep = "\r"))
  }

  # Precompute chromosome window edges once across all groups.
  n_chr <- nrow(k)
  chr_edges <- lapply(seq_len(n_chr), function(ci) {
    e <- window_edges(k$.end[ci] - k$.start[ci], window)
    lower <- e$lower + k$.start[ci] + 1
    upper <- e$upper + k$.start[ci]
    list(lower = lower, upper = upper, width = upper - lower + 1)
  })

  n_grp <- nrow(groups)
  parts <- vector("list", n_grp * n_chr)
  p_idx <- 0L
  for (gi in seq_len(n_grp)) {
    gk <- uniq_gkey[gi]
    for (ci in seq_len(n_chr)) {
      ce <- chr_edges[[ci]]
      lower <- ce$lower; upper <- ce$upper
      rows <- buckets[[paste(gk, k$.chr[ci], sep = "\r")]]
      result <- data.frame(Chr = k$.chr[ci], Start = lower, End = upper,
        Value = if (method == 'weighted_mean') NA_real_ else 0,
        Width = ce$width, N = 0L, N_valid = 0L, stringsAsFactors = FALSE)
      if (method == 'count') {
        if (length(rows)) {
          index <- findInterval(at[rows], upper, left.open = TRUE) + 1L
          result$N <- tabulate(index, nbins = length(upper))
          result$N_valid <- result$N
          result$Value <- as.numeric(result$N)
        }
        result$Rate <- result$Value / result$Width * 1e6
      } else if (length(rows)) {
        result <- summarise_interval_windows(result, start[rows], end[rows],
          values[rows], method, na.rm)
      }
      for (name in group) result[[name]] <- groups[[name]][gi]
      p_idx <- p_idx + 1L
      parts[[p_idx]] <- result
    }
  }
  out <- if (length(parts)) do.call(rbind, parts) else data.frame(
    Chr = character(), Start = numeric(), End = numeric(), Value = numeric(),
    Width = numeric(), N = integer(), N_valid = integer())
  if (!length(parts)) {
    if (method == 'count') out$Rate <- numeric()
    for (name in group) out[[name]] <- data[[name]][FALSE]
  }
  rownames(out) <- NULL
  attr(out, 'window_summary') <- window_metadata(method, window, group, value,
    count_position, na.rm, empty = if (method == 'weighted_mean') NA_real_ else 0)
  out
}

summarise_interval_windows <- function(result, start, end, values, method, na.rm) {
  index <- order(start, end)
  start <- start[index]; end <- end[index]; values <- values[index]
  w_start <- result$Start
  w_end <- result$End
  w_width <- result$Width
  n_w <- length(w_start)
  last <- findInterval(w_end, start)
  out_n <- integer(n_w)
  out_nv <- integer(n_w)
  out_val <- result$Value
  is_cov <- identical(method, 'coverage')
  active <- integer(); previous <- 0L
  for (wi in seq_len(n_w)) {
    if (last[wi] > previous) active <- c(active, seq.int(previous + 1L, last[wi]))
    previous <- last[wi]
    if (!length(active)) next
    active <- active[end[active] >= w_start[wi]]
    if (!length(active)) next
    lo <- pmax(start[active], w_start[wi])
    hi <- pmin(end[active], w_end[wi])
    v_act <- values[active]
    valid <- !is.na(v_act)
    out_n[wi] <- length(active)
    out_nv[wi] <- sum(valid)
    if (is_cov) {
      out_val[wi] <- union_covered_bp(lo, hi) / w_width[wi]
    } else if (any(valid) && (na.rm || all(valid))) {
      weights <- hi[valid] - lo[valid] + 1
      out_val[wi] <- sum(v_act[valid] * (weights / sum(weights)))
    }
  }
  result$N <- out_n
  result$N_valid <- out_nv
  result$Value <- out_val
  result
}

union_covered_bp <- function(start, end) {
  n <- length(start)
  if (!n) return(0)
  if (n == 1L) return(end[1L] - start[1L] + 1)
  ord <- order(start, end)
  start <- start[ord]; end <- end[ord]
  cm <- cummax(end)
  first_i <- which(c(TRUE, start[-1L] > cm[-n] + 1))
  last_i <- c(first_i[-1L] - 1L, n)
  sum(cm[last_i] - start[first_i] + 1)
}
