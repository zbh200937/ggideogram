# Complete ggplot/grob components anchored to chromosome coordinates.

StatChrInset <- ggplot2::ggproto(
  "StatChrInset", ggplot2::Stat,
  required_aes = c("chr", "plot", "position|start", "position|end"),
  compute_panel = function(data, scales, placement = "beside",
                           side = "right", gap = 0.25, track = NULL,
                           na.rm = FALSE) {
    position <- if ("position" %in% names(data)) {
      data$position
    } else {
      (data$start + data$end) / 2
    }
    data$ideogram_chr <- as.character(data$chr)
    data$ideogram_position <- position
    if (placement == "beside") {
      data$ideogram_side <- side
      data$ideogram_gap <- gap
      if (!is.null(track)) {
        data$ideogram_marker_track <- track
        data$ideogram_inset_track_edge <- TRUE
      }
    }
    data$x <- 0
    data$y <- 0
    data
  }
)

GeomChrInset <- ggplot2::ggproto(
  "GeomChrInset", ggplot2::Geom,
  required_aes = c("x", "y", "plot"),
  default_aes = ggplot2::aes(),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord, width, height,
                        placement = "beside", hjust = NULL, vjust = NULL,
                        clip = "on", overlap = "error",
                        na.rm = FALSE) {
    coordinates <- coord$transform(data, panel_params)
    keep <- is.finite(coordinates$x) & is.finite(coordinates$y)
    coordinates <- coordinates[keep, , drop = FALSE]
    if (!nrow(coordinates)) return(grid::nullGrob())
    justification <- inset_justification(
      coord, coordinates, placement, hjust, vjust)
    check_inset_overlap(
      coord, panel_params, coordinates, width, height,
      justification$hjust, justification$vjust,
      placement = placement, action = overlap)

    children <- lapply(seq_len(nrow(coordinates)), function(index) {
      child <- coordinates$plot[[index]]
      grob <- if (inherits(child, "ggplot")) {
        ggplot2::ggplotGrob(child)
      } else {
        child
      }
      grid::grobTree(
        grob,
        vp = grid::viewport(
          x = grid::unit(coordinates$x[index], "npc"),
          y = grid::unit(coordinates$y[index], "npc"),
          width = width, height = height,
          just = c(justification$hjust[index],
                   justification$vjust[index]),
          clip = clip
        )
      )
    })
    do.call(grid::grobTree, children)
  }
)

#' Anchor complete ggplot or grob components to chromosomes
#'
#' `geom_chr_inset()` keeps each child plot intact: its coordinate system,
#' theme and guides are rebuilt inside a dedicated viewport and do not enter
#' the host plot's scales. The child is never rasterized or decomposed into
#' guessed semantic layers.
#'
#' Map `chr`, `plot`, and either `position` or both `start`/`end`. Interval
#' anchors use their midpoint. A `plot` column must be a list-column containing
#' complete ggplot objects or grid grobs.
#'
#' @param mapping Aesthetic mapping containing the chromosome keys and `plot`.
#' @param data Inset data frame.
#' @param stat Semantic inset Stat; normally unchanged.
#' @param position Standard ggplot2 position adjustment; normally identity.
#' @param ... Additional layer parameters.
#' @param width,height One positive scalar [grid::unit()] each. Use `"npc"`
#'   for normalized host-panel dimensions or physical units such as `"mm"`.
#' @param placement `"beside"` offsets the anchor transversely;
#'   `"center"` anchors on the chromosome axis.
#' @param side,gap Fallback beside placement, with `gap` in chromosome-width
#'   units.
#' @param track Optional declared track whose inner boundary owns the inset
#'   anchor and whose width reserves outward component space.
#' @param hjust,vjust Viewport justification at the anchor. `NULL` uses a safe
#'   placement-aware default: a beside inset attaches its inner edge to the
#'   anchor and extends only outwards; a centred inset uses `0.5`.
#' @param clip Clip the child to its viewport.
#' @param overlap What to do when normalized-panel (`"npc"`) inset rectangles
#'   overlap one another or a chromosome body: `"error"` (default), `"warn"`,
#'   or `"allow"`. Declare a sufficiently wide inset [track()] to resolve an
#'   error instead of hiding the collision.
#' @param na.rm Remove missing required coordinates silently.
#' @param show.legend Insets do not participate in host guides; retained for a
#'   ggplot2-style layer signature and fixed to `FALSE` by default.
#' @param inherit.aes Whether to inherit host aesthetics.
#'
#' @return A ggplot2 layer using `GeomChrInset`.
#' @examples
#' chromosomes <- data.frame(Chr = "A", Start = 0, End = 100)
#' child <- ggplot2::ggplot(
#'   data.frame(x = 1:3, y = c(1, 3, 2)),
#'   ggplot2::aes(x, y)
#' ) + ggplot2::geom_line()
#' inset <- data.frame(Chr = "A", Pos = 50)
#' inset$Plot <- list(child)
#' ggideogram(chromosomes) +
#'   geom_chr_inset(
#'     data = inset,
#'     ggplot2::aes(chr = Chr, position = Pos, plot = Plot),
#'     width = grid::unit(0.25, "npc"),
#'     height = grid::unit(0.25, "npc")
#'   )
#' @export
geom_chr_inset <- function(
    mapping = NULL,
    data = NULL,
    stat = StatChrInset,
    position = "identity",
    ...,
    width,
    height,
    placement = c("beside", "center"),
    side = c("right", "left"),
    gap = 0.25,
    track = NULL,
    hjust = NULL,
    vjust = NULL,
    clip = c("on", "off"),
    overlap = c("error", "warn", "allow"),
    na.rm = FALSE,
    show.legend = FALSE,
    inherit.aes = FALSE) {
  placement <- match.arg(placement)
  side <- match.arg(side)
  clip <- match.arg(clip)
  overlap <- match.arg(overlap)
  check_inset_unit(width, "width")
  check_inset_unit(height, "height")
  check_nonnegative_layout(gap, "gap")
  validate_optional_marker_track(track)
  check_inset_just(hjust, "hjust")
  check_inset_just(vjust, "vjust")
  validate_inset_mapping(data, mapping)

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomChrInset,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = c(list(
      width = width, height = height,
      placement = placement, side = side, gap = gap, track = track,
      hjust = hjust, vjust = vjust, clip = clip, overlap = overlap,
      na.rm = na.rm
    ), list(...))
  )
}

check_inset_unit <- function(x, what) {
  if (!grid::is.unit(x) || length(x) != 1L ||
      !is.finite(as.numeric(x)) || as.numeric(x) <= 0) {
    stopf("`%s` must be one positive `grid::unit()` value.", what)
  }
  invisible(x)
}

check_inset_just <- function(x, what) {
  if (is.null(x)) return(invisible(x))
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stopf("`%s` must be one finite number.", what)
  }
  invisible(x)
}

inset_justification <- function(coord, data, placement, hjust, vjust) {
  if (!is.null(hjust) && !is.null(vjust)) {
    return(list(
      hjust = rep(hjust, nrow(data)),
      vjust = rep(vjust, nrow(data))))
  }
  if (placement == "center") {
    return(list(
      hjust = rep(hjust %||% 0.5, nrow(data)),
      vjust = rep(vjust %||% 0.5, nrow(data))))
  }

  side <- as.character(data$ideogram_side)
  if ("ideogram_marker_track" %in% names(data)) {
    spec <- track_table_row(
      coord$layout, as.character(data$ideogram_marker_track[1]))
    side <- rep(spec$side, nrow(data))
  }
  if (any(!side %in% c("left", "right"))) {
    stopf("Beside inset tracks must be on the left or right side.")
  }
  if (coord$layout$orientation == "vertical") {
    automatic_hjust <- ifelse(side == "right", 0, 1)
    automatic_vjust <- rep(0.5, nrow(data))
  } else {
    automatic_hjust <- rep(0.5, nrow(data))
    automatic_vjust <- ifelse(side == "right", 1, 0)
  }
  list(
    hjust = if (is.null(hjust)) automatic_hjust else rep(hjust, nrow(data)),
    vjust = if (is.null(vjust)) automatic_vjust else rep(vjust, nrow(data))
  )
}

check_inset_overlap <- function(
    coord, panel_params, data, width, height, hjust, vjust,
    placement, action) {
  if (action == "allow" || grid::unitType(width) != "npc" ||
      grid::unitType(height) != "npc") {
    return(invisible(FALSE))
  }
  width <- as.numeric(width)
  height <- as.numeric(height)
  rectangle <- data.frame(
    left = data$x - hjust * width,
    right = data$x + (1 - hjust) * width,
    bottom = data$y - vjust * height,
    top = data$y + (1 - vjust) * height
  )
  collision <- character()
  if (nrow(rectangle) > 1L) {
    pairs <- utils::combn(seq_len(nrow(rectangle)), 2)
    overlap <- apply(pairs, 2, function(pair) {
      first <- rectangle[pair[1], ]
      second <- rectangle[pair[2], ]
      first$left < second$right && first$right > second$left &&
        first$bottom < second$top && first$top > second$bottom
    })
    if (any(overlap)) {
      pair <- pairs[, which(overlap)[1]]
      collision <- c(
        collision,
        sprintf("insets %d and %d overlap", pair[1], pair[2]))
    }
  }

  if (placement == "beside") {
    body_min <- coord$transform(
      data.frame(x = coord$layout$chrom$.x_min,
                 y = coord$layout$chrom$.y_min), panel_params)
    body_max <- coord$transform(
      data.frame(x = coord$layout$chrom$.x_max,
                 y = coord$layout$chrom$.y_max), panel_params)
    for (inset_index in seq_len(nrow(rectangle))) {
      hit <- rectangle$left[inset_index] < body_max$x &
        rectangle$right[inset_index] > body_min$x &
        rectangle$bottom[inset_index] < body_max$y &
        rectangle$top[inset_index] > body_min$y
      if (any(hit)) {
        collision <- c(
          collision,
          sprintf("inset %d overlaps chromosome %s", inset_index,
                  coord$layout$chrom$.chr[which(hit)[1]]))
        break
      }
    }
  }
  if (!length(collision)) return(invisible(FALSE))

  message <- paste0(
    "Inset viewport collision: ", paste(collision, collapse = "; "), ".\n",
    "  Declare a wider inset track, reduce the `npc` size, or set ",
    "`overlap = \"allow\"` intentionally.")
  if (action == "error") stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  invisible(TRUE)
}

validate_inset_mapping <- function(data, mapping) {
  if (is.null(data) || is.null(mapping)) return(invisible(TRUE))
  if (!is.data.frame(data)) {
    stopf("Inset `data` must be a data frame.")
  }
  base <- mapped_fields(data, mapping, c("chr", "plot"), what = "mapping")
  has_position <- "position" %in% names(mapping)
  has_interval <- all(c("start", "end") %in% names(mapping))
  if (!has_position && !has_interval) {
    stopf("Inset `mapping` needs `position` or both `start` and `end`.")
  }
  if (has_position && has_interval) {
    stopf("Inset `mapping` must use either `position` or `start`/`end`, not both.")
  }
  if (has_position) {
    mapped_fields(data, mapping, "position", what = "mapping")
  } else {
    interval <- mapped_fields(data, mapping, c("start", "end"),
                              what = "mapping")
    if (!is.numeric(interval$start) || !is.numeric(interval$end) ||
        anyNA(interval$start) || anyNA(interval$end) ||
        any(interval$start > interval$end)) {
      stopf("Inset intervals require finite numeric `start <= end`.")
    }
  }
  plots <- base$plot
  if (!is.list(plots) ||
      any(!vapply(plots, function(plot) {
        inherits(plot, "ggplot") || grid::is.grob(plot)
      }, logical(1)))) {
    stopf("Inset `plot` must be a list-column of ggplot objects or grobs.")
  }
  invisible(TRUE)
}

#' Convert an ideogram to a standard grob
#'
#' This is a stable convenience wrapper around [ggplot2::ggplotGrob()] and does
#' not create a second renderer or fixed page layout.
#'
#' @param plot A [ggideogram()] plot.
#'
#' @return A standard gtable/grob.
#' @export
as_ideogram_grob <- function(plot) {
  if (!inherits(plot, "ggplot") ||
      !inherits(plot$coordinates$layout, "ideogram_layout_v2")) {
    stopf("`plot` must be a `ggideogram()` ggplot object.")
  }
  ggplot2::ggplotGrob(plot)
}
