# Small shared helpers with no rendering or device assumptions.

`%||%` <- function(x, y) if (is.null(x)) y else x

stopf <- function(...) stop(sprintf(...), call. = FALSE)

# ggplot2 passes the expression used on the right-hand side of `+` through
# `...`.  Keeping that detail in one helper lets our S3 methods match the
# public ggplot_add(object, plot, ...) generic while retaining useful layer
# and attachment names in diagnostics.
ggplot_add_object_name <- function(..., fallback) {
  dots <- list(...)
  candidate <- dots$object_name
  if (is.null(candidate) && length(dots)) candidate <- dots[[1L]]
  if (is.null(candidate) || length(candidate) != 1L || is.na(candidate)) {
    return(fallback)
  }
  candidate <- as.character(candidate)
  if (!nzchar(candidate)) fallback else candidate
}

normalise_colour <- function(x) {
  x <- as.character(x)
  bare <- !is.na(x) & grepl("^([0-9A-Fa-f]{6}|[0-9A-Fa-f]{8})$", x)
  x[bare] <- paste0("#", x[bare])
  values <- unique(x[!is.na(x)])
  valid <- vapply(values, function(value) {
    !inherits(try(grDevices::col2rgb(value), silent = TRUE), "try-error")
  }, logical(1))
  if (any(!valid)) {
    stopf("Not a valid colour: %s.",
          paste0("`", values[!valid], "`", collapse = ", "))
  }
  unname(x)
}

NICE_AXIS_STEPS <- c(1, 2, 5, 10)

axis_breaks <- function(span, n = 6) {
  if (!is.finite(span) || span <= 0) return(0)
  raw <- span / max(n, 1)
  magnitude <- 10^floor(log10(raw))
  step <- magnitude * NICE_AXIS_STEPS[
    which.min(abs(NICE_AXIS_STEPS * magnitude - raw))
  ]
  seq(0, span, by = step)
}

axis_unit <- function(step, units) {
  table <- c(bp = 1, kb = 1e3, Mb = 1e6, Gb = 1e9)
  if (units != "auto") {
    return(list(div = unname(table[units]), suffix = units))
  }
  unit <- names(table)[max(which(step >= table | table == 1))]
  list(div = unname(table[unit]), suffix = unit)
}

axis_labels <- function(breaks, div, suffix, digits = NULL) {
  if (is.null(digits)) {
    step <- if (length(breaks) > 1L) min(diff(breaks)) else breaks[1]
    digits <- max(0, min(3, ceiling(-log10(max(step, 1) / div))))
  }
  paste0(
    formatC(breaks / div, format = "f", digits = digits),
    " ", suffix
  )
}

# Least-squares one-dimensional separation. Subtracting the required spacing
# turns the constraints into ordinary isotonic regression; clamping the fitted
# isotonic values preserves their order and enforces both bounds.
repel_1d <- function(position, distance, lo = -Inf, hi = Inf) {
  if (!length(position)) return(position)
  order <- order(position, seq_along(position))
  sorted <- position[order]
  count <- length(sorted)
  if ((count - 1) * distance > hi - lo) {
    adjusted <- seq(lo, hi, length.out = count)
  } else {
    offset <- (seq_len(count) - 1) * distance
    fitted <- stats::isoreg(sorted - offset)$yf
    lower <- lo
    upper <- hi - (count - 1) * distance
    fitted <- pmin(pmax(fitted, lower), upper)
    adjusted <- fitted + offset
  }
  result <- numeric(count)
  result[order] <- adjusted
  result
}
