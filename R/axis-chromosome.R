# Chromosome ideograms as native ggplot2 discrete-axis guides.

#' Use chromosome ideograms as an x or y axis
#'
#' `scale_x_chromosome()` and `scale_y_chromosome()` are ordinary discrete
#' position scales whose guide draws one chromosome in the axis-label slot for
#' each visible chromosome break. The host layers remain untouched: a native
#' [ggplot2::geom_col()], [ggplot2::geom_line()] or
#' [ggplot2::geom_boxplot()] continues to use the plot's ordinary x/y data
#' coordinates.
#'
#' Chromosomes are matched by the raw scale break, not by the formatted display
#' label. This means `labels` may safely add prefixes or other presentation
#' text without changing biological alignment. The compact guide displays the
#' chromosome outline and centromere; use [ggideogram()] when full cytoband
#' detail or bp-position tracks are required.
#'
#' @param data A karyotype data frame or an [as_ideogram_data()] object.
#' @param mapping,centromere Passed to [as_ideogram_data()] for a data-frame
#'   input.
#' @param name,breaks,labels,limits,expand,position,drop,na.translate Standard
#'   discrete position-scale arguments. The default `limits = waiver()` uses
#'   the karyotype order and retains every chromosome. Supply `limits = NULL`
#'   for data-trained limits and standard `drop` behaviour.
#' @param guide `NULL` constructs [guide_chromosome_axis()]. Any explicitly
#'   supplied ggplot2 guide, including `"none"`, is passed through unchanged.
#' @param axis_length Maximum chromosome length in the guide. A one-dimensional
#'   absolute [grid::unit()] such as `unit(18, "mm")`; shorter chromosomes are
#'   drawn in their true relative length.
#' @param body_width Chromosome body width as an absolute [grid::unit()].
#' @param label_gap Clear physical gap between the chromosome and its displayed
#'   name.
#' @param show_labels Draw the scale's formatted labels beyond the chromosome
#'   bodies. Their typography comes from the normal `axis.text.*` theme
#'   element.
#' @param fill,colour,linewidth Chromosome body fill, outline colour and
#'   ggplot2-style line width.
#' @param curve_points Number of profile samples used for rounded caps and the
#'   centromere waist.
#' @param ... Additional arguments passed to [ggplot2::scale_x_discrete()] or
#'   [ggplot2::scale_y_discrete()].
#'
#' @return A discrete ggplot2 position scale.
#' @examples
#' data(human_karyotype, package = "ggideogram")
#' counts <- data.frame(
#'   Chr = human_karyotype$Chr,
#'   Count = seq_len(nrow(human_karyotype))
#' )
#'
#' ggplot2::ggplot(counts, ggplot2::aes(Chr, Count)) +
#'   ggplot2::geom_col() +
#'   scale_x_chromosome(human_karyotype)
#'
#' ggplot2::ggplot(counts, ggplot2::aes(Count, Chr)) +
#'   ggplot2::geom_col() +
#'   scale_y_chromosome(human_karyotype)
#' @export
scale_x_chromosome <- function(
    data,
    mapping = NULL,
    centromere = NULL,
    name = ggplot2::waiver(),
    breaks = ggplot2::waiver(),
    labels = ggplot2::waiver(),
    limits = ggplot2::waiver(),
    expand = ggplot2::waiver(),
    position = "bottom",
    drop = TRUE,
    na.translate = FALSE,
    guide = NULL,
    axis_length = grid::unit(18, "mm"),
    body_width = grid::unit(1.8, "mm"),
    label_gap = grid::unit(1, "mm"),
    show_labels = TRUE,
    fill = "#F7F7F7",
    colour = "#4D4D4D",
    linewidth = 0.4,
    curve_points = 32,
    ...) {
  chromosome_position_scale(
    axis = "x", data = data, mapping = mapping,
    centromere = centromere, name = name, breaks = breaks,
    labels = labels, limits = limits, expand = expand,
    position = position, drop = drop, na.translate = na.translate,
    guide = guide, axis_length = axis_length, body_width = body_width,
    label_gap = label_gap, show_labels = show_labels, fill = fill,
    colour = colour, linewidth = linewidth, curve_points = curve_points,
    ...
  )
}

#' @rdname scale_x_chromosome
#' @export
scale_y_chromosome <- function(
    data,
    mapping = NULL,
    centromere = NULL,
    name = ggplot2::waiver(),
    breaks = ggplot2::waiver(),
    labels = ggplot2::waiver(),
    limits = ggplot2::waiver(),
    expand = ggplot2::waiver(),
    position = "left",
    drop = TRUE,
    na.translate = FALSE,
    guide = NULL,
    axis_length = grid::unit(18, "mm"),
    body_width = grid::unit(1.8, "mm"),
    label_gap = grid::unit(1, "mm"),
    show_labels = TRUE,
    fill = "#F7F7F7",
    colour = "#4D4D4D",
    linewidth = 0.4,
    curve_points = 32,
    ...) {
  chromosome_position_scale(
    axis = "y", data = data, mapping = mapping,
    centromere = centromere, name = name, breaks = breaks,
    labels = labels, limits = limits, expand = expand,
    position = position, drop = drop, na.translate = na.translate,
    guide = guide, axis_length = axis_length, body_width = body_width,
    label_gap = label_gap, show_labels = show_labels, fill = fill,
    colour = colour, linewidth = linewidth, curve_points = curve_points,
    ...
  )
}

chromosome_position_scale <- function(
    axis, data, mapping, centromere, name, breaks, labels, limits,
    expand, position, drop, na.translate, guide, axis_length,
    body_width, label_gap, show_labels, fill, colour, linewidth,
    curve_points, ...) {
  semantic <- chromosome_axis_semantic(data, mapping, centromere)
  if (inherits(limits, "waiver")) {
    limits <- semantic$karyotype$.chr
  }
  if (is.null(guide)) {
    guide <- guide_chromosome_axis(
      semantic,
      axis_length = axis_length,
      body_width = body_width,
      label_gap = label_gap,
      show_labels = show_labels,
      fill = fill,
      colour = colour,
      linewidth = linewidth,
      curve_points = curve_points
    )
  }
  scale <- if (identical(axis, "x")) {
    ggplot2::scale_x_discrete
  } else {
    ggplot2::scale_y_discrete
  }
  scale(
    name = name, breaks = breaks, labels = labels, limits = limits,
    expand = expand, guide = guide, position = position,
    drop = drop, na.translate = na.translate, ...
  )
}

chromosome_axis_semantic <- function(data, mapping, centromere) {
  if (inherits(data, "ideogram_data")) {
    supplied <- list(mapping = mapping, centromere = centromere)
    supplied <- names(supplied)[!vapply(supplied, is.null, logical(1))]
    if (length(supplied)) {
      stopf("An `ideogram_data` input cannot also supply %s.",
            paste0("`", supplied, "`", collapse = ", "))
    }
    return(data)
  }
  as_ideogram_data(data, mapping = mapping, centromere = centromere)
}

#' Draw chromosome ideograms in a discrete-axis guide
#'
#' This guide is normally constructed by [scale_x_chromosome()] or
#' [scale_y_chromosome()]. It extends ggplot2's public [ggplot2::GuideAxis]
#' interface and occupies the normal axis-label slot; it does not inspect or
#' edit a built gtable.
#'
#' @param data A karyotype data frame or [as_ideogram_data()] object.
#' @param mapping,centromere Passed to [as_ideogram_data()].
#' @param axis_length,body_width,label_gap,show_labels,fill,colour,linewidth,curve_points
#'   See [scale_x_chromosome()].
#' @param title,theme,order,position Standard ggplot2 guide arguments.
#'
#' @return A ggplot2 guide object.
#' @export
guide_chromosome_axis <- function(
    data,
    mapping = NULL,
    centromere = NULL,
    axis_length = grid::unit(18, "mm"),
    body_width = grid::unit(1.8, "mm"),
    label_gap = grid::unit(1, "mm"),
    show_labels = TRUE,
    fill = "#F7F7F7",
    colour = "#4D4D4D",
    linewidth = 0.4,
    curve_points = 32,
    title = ggplot2::waiver(),
    theme = NULL,
    order = 0,
    position = ggplot2::waiver()) {
  semantic <- chromosome_axis_semantic(data, mapping, centromere)
  check_axis_absolute_unit(axis_length, "axis_length", positive = TRUE)
  check_axis_absolute_unit(body_width, "body_width", positive = TRUE)
  check_axis_absolute_unit(label_gap, "label_gap", positive = FALSE)
  if (!is.logical(show_labels) || length(show_labels) != 1L ||
      is.na(show_labels)) {
    stopf("`show_labels` must be TRUE or FALSE.")
  }
  check_nonnegative_layout(linewidth, "linewidth")
  curve_points <- validate_curve_points(curve_points)

  ggplot2::new_guide(
    title = title,
    theme = theme,
    name = "chromosome_axis",
    position = position,
    direction = NULL,
    angle = NULL,
    n.dodge = 1,
    minor.ticks = FALSE,
    cap = "none",
    order = order,
    check.overlap = FALSE,
    semantic = semantic,
    axis_length = axis_length,
    body_width = body_width,
    label_gap = label_gap,
    show_labels = show_labels,
    fill = fill,
    colour = colour,
    linewidth = linewidth,
    curve_points = curve_points,
    available_aes = c("x", "y"),
    super = GuideChromosomeAxis
  )
}

GuideChromosomeAxis <- ggplot2::ggproto(
  "GuideChromosomeAxis", ggplot2::GuideAxis,
  params = c(
    ggplot2::GuideAxis$params,
    list(
      semantic = NULL,
      axis_length = NULL,
      body_width = NULL,
      label_gap = NULL,
      show_labels = TRUE,
      fill = "#F7F7F7",
      colour = "#4D4D4D",
      linewidth = 0.4,
      curve_points = 32
    )
  ),
  extract_key = function(self, scale, aesthetic, ...) {
    key <- ggplot2::GuideAxis$extract_key(scale, aesthetic, ...)
    if (is.null(key) || !nrow(key)) return(key)
    breaks <- scale$get_breaks()
    if (length(breaks) != nrow(key)) {
      stopf("Chromosome-axis breaks could not be matched to the guide key.")
    }
    key$.chr <- as.character(breaks)
    key
  },
  build_labels = function(self, key, elements, params) {
    if (".type" %in% names(key)) {
      key <- key[key$.type == "major", , drop = FALSE]
    }
    if (!nrow(key)) return(list(grid::nullGrob()))
    chromosome_axis_strip(
      key = key,
      text_element = elements$text,
      params = params
    )
  }
)

chromosome_axis_strip <- function(key, text_element, params) {
  karyotype <- params$semantic$karyotype
  chromosome <- as.character(key$.chr)
  index <- match(chromosome, karyotype$.chr)
  if (anyNA(index)) {
    unknown <- unique(chromosome[is.na(index)])
    stopf("Chromosome-axis break%s not found in `data`: %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  karyotype <- karyotype[index, , drop = FALSE]
  position <- key[[params$aes]]
  labels <- as.character(key$.label %||% character())
  has_labels <- params$show_labels && length(labels) == length(position) &&
    any(!is.na(labels) & nzchar(labels))
  max_length_mm <- axis_unit_mm(params$axis_length, "axis_length")
  body_width_mm <- axis_unit_mm(params$body_width, "body_width")
  length_ratio <- (karyotype$.end - karyotype$.start) /
    max(params$semantic$karyotype$.end -
          params$semantic$karyotype$.start)
  body_length_mm <- max_length_mm * length_ratio

  body_grobs <- lapply(seq_len(nrow(karyotype)), function(i) {
    chromosome_axis_body_grob(
      karyotype[i, , drop = FALSE],
      axis_position = position[i],
      guide_position = params$position,
      body_length_mm = body_length_mm[i],
      body_width_mm = min(body_width_mm, body_length_mm[i]),
      fill = params$fill,
      colour = params$colour,
      linewidth = params$linewidth,
      curve_points = params$curve_points
    )
  })

  text <- chromosome_axis_text_grob(
    position = position,
    labels = labels,
    text_element = text_element,
    guide_position = params$position,
    show_labels = has_labels
  )
  text_width <- if (has_labels) grid::grobWidth(text) else
    grid::unit(0, "mm")
  text_height <- if (has_labels) grid::grobHeight(text) else
    grid::unit(0, "mm")
  label_gap <- if (has_labels) params$label_gap else grid::unit(0, "mm")
  is_vertical_axis <- params$position %in% c("left", "right")
  strip_width <- if (is_vertical_axis) {
    params$axis_length + label_gap + text_width
  } else {
    grid::unit(1, "npc")
  }
  strip_height <- if (is_vertical_axis) {
    grid::unit(1, "npc")
  } else {
    params$axis_length + label_gap + text_height
  }
  children <- c(body_grobs, if (has_labels) list(text) else list())
  strip <- grid::gTree(
    children = do.call(grid::gList, children),
    width = strip_width,
    height = strip_height,
    vp = grid::viewport(xscale = c(0, 1), yscale = c(0, 1)),
    cl = "chromosome_axis_strip"
  )
  list(strip)
}

chromosome_axis_text_grob <- function(
    position, labels, text_element, guide_position, show_labels) {
  if (!show_labels) return(grid::nullGrob())
  no_margin <- ggplot2::margin(0, 0, 0, 0)
  if (guide_position == "bottom") {
    ggplot2::element_grob(
      text_element, label = labels,
      x = grid::unit(position, "native"), y = grid::unit(0, "npc"),
      hjust = 0.5, vjust = 0, margin = no_margin)
  } else if (guide_position == "top") {
    ggplot2::element_grob(
      text_element, label = labels,
      x = grid::unit(position, "native"), y = grid::unit(1, "npc"),
      hjust = 0.5, vjust = 1, margin = no_margin)
  } else if (guide_position == "left") {
    ggplot2::element_grob(
      text_element, label = labels,
      x = grid::unit(0, "npc"), y = grid::unit(position, "native"),
      hjust = 0, vjust = 0.5, margin = no_margin)
  } else {
    ggplot2::element_grob(
      text_element, label = labels,
      x = grid::unit(1, "npc"), y = grid::unit(position, "native"),
      hjust = 1, vjust = 0.5, margin = no_margin)
  }
}

chromosome_axis_body_grob <- function(
    karyotype, axis_position, guide_position,
    body_length_mm, body_width_mm, fill, colour,
    linewidth, curve_points) {
  profile <- chromosome_axis_profile(
    karyotype, body_length_mm, body_width_mm, curve_points)
  vertical <- guide_position %in% c("bottom", "top")
  if (vertical) {
    x <- profile$cross
    y <- if (guide_position == "bottom") {
      1 - profile$long
    } else {
      profile$long
    }
    vp <- grid::viewport(
      x = grid::unit(axis_position, "native"),
      y = grid::unit(if (guide_position == "bottom") 1 else 0, "npc"),
      width = grid::unit(body_width_mm, "mm"),
      height = grid::unit(body_length_mm, "mm"),
      just = c("centre", if (guide_position == "bottom") "top" else "bottom")
    )
  } else {
    x <- if (guide_position == "left") 1 - profile$long else profile$long
    y <- profile$cross
    vp <- grid::viewport(
      x = grid::unit(if (guide_position == "left") 1 else 0, "npc"),
      y = grid::unit(axis_position, "native"),
      width = grid::unit(body_length_mm, "mm"),
      height = grid::unit(body_width_mm, "mm"),
      just = c(if (guide_position == "left") "right" else "left", "centre")
    )
  }
  grid::grobTree(
    grid::polygonGrob(
      x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
      gp = grid::gpar(
        fill = fill, col = colour,
        lwd = linewidth * ggplot2::.pt,
        linejoin = "round"
      )
    ),
    vp = vp
  )
}

chromosome_axis_profile <- function(
    karyotype, body_length_mm, body_width_mm, curve_points) {
  radius <- body_width_mm / 2
  cap_start <- seq(0, radius, length.out = curve_points)
  cap_end <- seq(body_length_mm - radius, body_length_mm,
                 length.out = curve_points)
  centromere <- numeric()
  has_centromere <- !is.na(karyotype$.centromere_start) &&
    !is.na(karyotype$.centromere_end) &&
    karyotype$.centromere_end > karyotype$.centromere_start
  if (has_centromere) {
    ce_start <- body_length_mm *
      (karyotype$.centromere_start - karyotype$.start) /
      (karyotype$.end - karyotype$.start)
    ce_end <- body_length_mm *
      (karyotype$.centromere_end - karyotype$.start) /
      (karyotype$.end - karyotype$.start)
    centromere <- c(ce_start, mean(c(ce_start, ce_end)), ce_end)
  }
  longitudinal <- sort(unique(c(0, body_length_mm, cap_start, cap_end,
                                centromere)))
  half_width <- rep(radius, length(longitudinal))
  first_cap <- longitudinal < radius
  half_width[first_cap] <- sqrt(pmax(
    0, radius^2 - (radius - longitudinal[first_cap])^2))
  last_cap <- longitudinal > body_length_mm - radius
  half_width[last_cap] <- sqrt(pmax(
    0, radius^2 -
      (longitudinal[last_cap] - (body_length_mm - radius))^2))
  if (has_centromere) {
    at_waist <- longitudinal >= ce_start & longitudinal <= ce_end
    half_width[at_waist] <- radius *
      abs(2 * (longitudinal[at_waist] - ce_start) /
            (ce_end - ce_start) - 1)
  }
  one_side <- data.frame(
    long = longitudinal / body_length_mm,
    cross = 0.5 + half_width / body_width_mm
  )
  other_side <- data.frame(
    long = rev(longitudinal / body_length_mm),
    cross = rev(0.5 - half_width / body_width_mm)
  )
  profile <- rbind(one_side, other_side)
  rbind(profile, profile[1, , drop = FALSE])
}

check_axis_absolute_unit <- function(x, what, positive) {
  if (!inherits(x, "unit") || length(x) != 1L) {
    stopf("`%s` must be one absolute grid unit.", what)
  }
  type <- grid::unitType(x)
  relative <- c("npc", "native", "null", "snpc", "strwidth", "strheight",
                "grobwidth", "grobheight")
  if (type %in% relative) {
    stopf("`%s` must use an absolute physical unit, not `%s`.", what, type)
  }
  value <- tryCatch(
    axis_unit_mm(x, what),
    error = function(error) NA_real_
  )
  invalid <- !is.finite(value) || if (positive) value <= 0 else value < 0
  if (invalid) {
    stopf("`%s` must be one %sfinite physical length.", what,
          if (positive) "positive " else "non-negative ")
  }
  invisible(x)
}

axis_unit_mm <- function(x, what) {
  value <- grid::convertWidth(x, "mm", valueOnly = TRUE)
  if (length(value) != 1L) {
    stopf("`%s` must resolve to one physical length.", what)
  }
  value
}

#' @method widthDetails chromosome_axis_strip
#' @importFrom grid widthDetails
#' @export
widthDetails.chromosome_axis_strip <- function(x) x$width

#' @method heightDetails chromosome_axis_strip
#' @importFrom grid heightDetails
#' @export
heightDetails.chromosome_axis_strip <- function(x) x$height
