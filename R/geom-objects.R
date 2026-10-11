# Public object-oriented grammar; source components can be replayed on relayout.

new_ideogram_object <- function(component) {
  structure(list(component = component), class = "ggideogram_object")
}

#' Draw a collection of chromosomes
#'
#' `geom_chr()` is the chromosome skeleton layer. On a plain ggplot it creates
#' the shared chromosome coordinate system from `data` or the plot data. On
#' an ideogram it adds the selected chromosome component to the same layout.
#' Select several chromosomes with `chr`, rather than making one function call
#' per chromosome. Names and bp axes use the same source coordinates.
#' When body and fill are added together to a plain ggplot, the plot data
#' supply the karyotype and the component data/mapping supply source intervals.
#'
#' @param mapping Karyotype mapping for a new plot, or interval mapping for
#'   `component = "fill"`.
#' @param data Karyotype data for a new skeleton, or source intervals for fill.
#'   If NULL, a new skeleton uses the plot's existing data.
#' @param chr NULL for all chromosomes, or a vector of chromosome keys.
#' @param component The component: `"body"`, `"name"`, `"axis"`, `"band"`,
#'   or `"fill"`. Several components can be supplied together.
#' @param names Optional logical to include chromosome names with the body.
#'   NULL uses the constructor default on a new plot.
#' @param axis Include the shared bp axis with the body.
#' @param track Track identifier for `component = "fill"`; inherited inside
#'   [geom_track()].
#' @param ... Component styling. A new plot also accepts [ggideogram()] layout
#'   arguments such as `orientation`. Names accept size/gap/position; axes
#'   accept side/breaks/units. Bodies accept fill/colour/linewidth/alpha/linetype.
#'   Unsupported parameters and visual mappings outside `component = "fill"`
#'   raise an error.
#' @return An additive chromosome component returning a standard ggplot.
#' @examples
#' k <- data.frame(Chr = c("A", "B"), Start = 0, End = c(100, 80))
#' ggplot2::ggplot(k) + geom_chr(chr = c("A", "B"),
#'   orientation = "horizontal", axis = TRUE)
#' @export
geom_chr <- function(mapping = NULL, data = NULL, chr = NULL,
    component = "body", names = NULL, axis = FALSE, track = NULL, ...) {
  component <- match.arg(component, c("body", "name", "axis", "band", "fill"), several.ok = TRUE)
  new_ideogram_object(structure(list(mapping = mapping, data = data, chr = chr,
    components = component, names = names, axis = axis, track = track, params = list(...)),
    class = "ggideogram_chr_component"))
}

#' Draw gene or transcript models in a chromosome track
#'
#' Map chr/start/end/type/strand and gene or transcript. Complete source models
#' are selected and clipped to a local view during preparation. Native segments
#' and rectangles draw backbones, exons, CDS and UTRs. In a [geom_track()] scope,
#' data, mapping and the track name can all be inherited.
#'
#' @param mapping,data Annotation mappings and source data.
#' @param track Optional declared track identifier.
#' @param mode Draw one lane per gene or transcript.
#' @param ... Gene styling: optional lane_order (complete identifier vector),
#'   block_height = 0.35 (fraction of one lane), native grid::arrow() or NULL,
#'   labels = TRUE, label_size = 3 (mm), label_gap = 0.3 (em),
#'   label_colour = "#202020", label_family = NULL (inherits `base_family`),
#'   label_fontface = "plain" (use "italic" for gene symbols), and standard
#'   aesthetics. Intron arrows follow strand; arrow_min_bp, arrow_spacing_bp
#'   and arrow_margin_bp control their genomic spacing.
#' @return An additive gene-model component.
#' @export
geom_genemodel <- function(mapping = NULL, data = NULL, track = NULL,
    mode = c("gene", "transcript"), ...) {
  fun <- if (match.arg(mode) == "gene") geom_chr_gene else geom_chr_transcript
  new_ideogram_object(fun(mapping, data, track, ...))
}

#' Draw native points, text or intervals at chromosome loci
#'
#' Ordinary point/text layers retain their Geom and physical styling. Map chr
#' and position, or chr/start/end for 1-based closed intervals. Interval display
#' boundaries are start - 1 and end, retaining single-base features.
#' `geom = "link"` draws a native
#' leader from the true genomic anchor to its displayed marker position.
#' Text layout is a property of this same locus annotation: `position =
#' "identity"` keeps its position, `"spread"` packs text in genomic order,
#' and `"repel"` uses the optional ggrepel for free text-box repulsion.
#' Spreading aligns the near text edges on one straight baseline or circular
#' arc, with readable facing and native segment leaders to the source loci.
#'
#' @param mapping,data Source locus mappings and data.
#' @param geom A point/text/segment layer or constructor; also `"point"`,
#'   `"text"`, `"interval"`, or `"link"`.
#' @param track Optional declared track identifier, inherited within a track.
#' @param track_position Relative transverse position in a track: 0 is the
#'   near edge, 1 the far edge, and 0.5 (default) the middle. This is annotation
#'   placement, independent of a track's scientific values.
#' @param side,gap Beside placement in chromosome-body-width units.
#' @param position A standard Position, including [position_chr_repel()], or
#'   `"spread"` / `"repel"` for text.
#' @param label_width Space reserved for spread/repel text without an explicit
#'   track, in chromosome-body widths. NULL selects 32 for linear/outer text
#'   and 14 for inner circular text. An explicit track owns its width.
#' @param ... Native geom parameters and aesthetics. Ordered text also accepts
#'   padding = 0.6 mm, label_gap = 0 em, connection_height = 5 mm,
#'   segment_colour = "#626A73", and segment_linewidth = 0.22 mm. Its default
#'   gap is 0.2 body widths. Text defaults to size = 3 mm. Free repulsion
#'   accepts columns, seed, max.iter, box.padding and overflow.
#'   Explicit arguments override an existing layer's styling and position;
#'   omitted arguments retain the layer's settings. Native layers must use
#'   `stat = "identity"`; summarise the source data before adding loci.
#' @return An additive locus component using native point/text/segment geoms
#'   or a compatible text extension.
#' @export
geom_locus <- function(mapping = NULL, data = NULL, geom = ggplot2::geom_point,
    track = NULL, side = c("right", "left", "inner", "outer"), gap = 0.25,
    position = "identity", track_position = 0.5, label_width = NULL, ...) {
  check_track_position(track_position)
  side <- canonical_chr_side(match.arg(side))
  if (is.character(position) && length(position) == 1L && position %in% c("spread", "repel")) {
    if (missing(gap)) gap <- 0.2
    return(new_ideogram_object(locus_text_component(mapping, data, geom,
      track, side, gap, label_width, position, list(...))))
  }
  if (is.character(geom)) {
    geom <- match.arg(geom, c("point", "text", "interval", "link"))
    fun <- switch(geom, point = geom_chr_marker, text = geom_chr_text,
      interval = geom_chr_interval, link = geom_chr_link)
    params <- list(...)
    if (!is.null(params$stat) && !identical(params$stat, "identity") &&
        !inherits(params$stat, "StatIdentity")) {
      stopf("Locus layers require stat = 'identity'; summarise source data before adding loci.")
    }
    params$stat <- NULL
    args <- c(list(mapping = mapping, data = data, track = track, side = side,
      gap = gap, position = position, track_position = track_position), params)
    if (geom == "text" && is.null(args$size) && is.null(mapping$size)) args$size <- ideogram_text_size("label")
    if (geom == "interval" && !is.null(track)) args$placement <- "beside"
    return(new_ideogram_object(do.call(fun, args)))
  }
  prototype <- locus_native_layer(geom, list(...),
    if (missing(position)) NULL else position)
  if (!inherits(prototype, "Layer")) stopf("`geom` must be a point, text, or segment layer.")
  if (!inherits(prototype$stat, "StatIdentity")) {
    stopf("Locus layers require stat = 'identity'; summarise source data before adding loci.")
  }
  mapping <- scope_mapping(mapping, prototype$mapping)
  data <- scope_layer_data(prototype$data, data)
  if (!is.null(mapping$x) && is.null(mapping$position)) mapping$position <- mapping$x
  if (inherits(prototype$geom, "GeomSegment") && all(c("start", "end") %in% names(mapping))) {
    base <- geom_chr_interval(mapping, data, track = track, side = side,
      gap = gap, placement = if (is.null(track)) "overlay" else "beside", track_position = track_position)
  } else if (inherits(prototype$geom, "GeomPoint") || inherits(prototype$geom, "GeomText")) {
    base <- geom_chr_marker(mapping, data, track = track, side = side, gap = gap,
      track_position = track_position)
  } else stopf("`geom` must be a point/text layer, or an interval segment layer.")
  base$geom <- prototype$geom
  base$position <- if (inherits(base$stat, "StatChrInterval")) {
    circular_segment_position(prototype$position)
  } else locus_position(prototype$position)
  base$geom_params <- prototype$geom_params
  base$aes_params <- prototype$aes_params
  if (inherits(prototype$geom, "GeomText") && is.null(base$aes_params$size) && is.null(mapping$size)) {
    base$aes_params$size <- ideogram_text_size("label")
  }
  base$show.legend <- prototype$show.legend
  new_ideogram_object(base)
}

locus_native_layer <- function(geom, params, position = NULL) {
  if (is.function(geom)) {
    return(do.call(geom, c(list(position = position %||% "identity"), params)))
  }
  if (!inherits(geom, "Layer") || (!length(params) && is.null(position))) return(geom)
  settings <- c("stat", "show.legend", "inherit.aes", "check.aes", "check.param")
  args <- list(geom = geom$geom, stat = geom$stat, mapping = geom$mapping,
    data = if (inherits(geom$data, "waiver")) NULL else geom$data,
    position = position %||% geom$position,
    show.legend = geom$show.legend, inherit.aes = geom$inherit.aes)
  args <- utils::modifyList(args, params[intersect(names(params), settings)])
  defaults <- c(geom$stat_params, geom$geom_params, geom$aes_params)
  defaults <- defaults[!duplicated(names(defaults), fromLast = TRUE)]
  args$params <- utils::modifyList(defaults, params[!names(params) %in% settings])
  do.call(ggplot2::layer, args)
}

locus_text_component <- function(mapping, data, geom, track, side, gap,
    width, position, params) {
  fun <- if (position == "spread") geom_chr_labels else geom_chr_text_repel
  controls <- setdiff(union(names(formals(geom_chr_labels)),
    names(formals(geom_chr_text_repel))),
    c("mapping", "data", "...", "stat", "position", "parse", "na.rm", "show.legend", "inherit.aes"))
  if (is.character(geom)) {
    if (!identical(geom, "text")) stopf("`position = '%s'` requires a text geom.", position)
    geom <- ggplot2::geom_text
  }
  prototype <- locus_native_layer(geom, params[!names(params) %in% controls], "identity")
  if (!inherits(prototype, "Layer") || !inherits(prototype$geom, "GeomText")) {
    stopf("`position = '%s'` requires a text geom.", position)
  }
  if (!inherits(prototype$stat, "StatIdentity")) stopf("Locus text layout requires stat = 'identity'.")
  mapping <- scope_mapping(mapping, prototype$mapping)
  data <- scope_layer_data(prototype$data, data)
  if (!is.null(mapping$x) && is.null(mapping$position)) mapping$position <- mapping$x
  args <- utils::modifyList(utils::modifyList(prototype$geom_params, prototype$aes_params), params)
  native <- args[intersect(names(args), prototype$geom$parameters())]
  args <- args[!names(args) %in% setdiff(names(native), c("parse", "na.rm"))]
  args$show.legend <- args$show.legend %||% prototype$show.legend
  args$stat <- "identity"
  object <- do.call(fun, c(list(mapping = mapping, data = data, track = track,
    side = side, gap = gap, width = width), args))
  object$text_geom <- prototype$geom
  object$text_params <- native
  object
}

#' Connect chromosome loci or genomic intervals
#'
#' Endpoints can be mapped directly, or resolved by their IDs from [chr_nodes()]
#' using `aes(from = ..., to = ...)`. Same-chromosome local relations and
#' reversed interval relations default to smooth arcs. Set curvature = 0 for
#' straight connections. Curves use native segments; interval bands use native
#' polygons and retain their original genomic bounds.
#'
#' @param mapping,data Edge mappings and source edge data.
#'   Direct point mappings are chr1/position1/chr2/position2; interval mappings
#'   are chr1/start1/end1/chr2/start2/end2, with optional orientation (+/-).
#'   ID mappings are from/to and use the supplied nodes table.
#' @param type Point connections or interval bands.
#' @param nodes Optional source node table from [chr_nodes()].
#' @param ... Endpoint side1/side2, gap, bend/curvature, curve_points,
#'   orientation for intervals, and native colour/fill/alpha/linewidth/arrow.
#' @return An additive native connection or band component.
#' @export
geom_chrlink <- function(mapping = NULL, data = NULL,
    type = c("point", "interval"), nodes = NULL, ...) {
  fun <- if (match.arg(type) == "point") geom_chr_connection else geom_chr_synteny
  new_ideogram_object(fun(mapping, data, nodes = nodes, ...))
}

#' Anchor a complete ggplot or grob to a chromosome locus
#'
#' Complete plots preserve their own coordinate system, theme and guides.
#' For plots that should share genomic coordinates, use [geom_track()] with
#' ordinary layers instead.
#'
#' @param plot A complete ggplot or grob; NULL when mapping an existing plot
#'   list-column.
#' @param mapping,data Chromosome anchor mapping and data.
#' @param ... width and height are required positive scalar grid::unit()
#'   values, either physical units such as "mm" or host-panel proportions
#'   "npc". Optional track references an already declared beside region;
#'   its width is independent of the child viewport. placement defaults
#'   to "beside", side to "right" (outer for circles), gap to 0.25 body widths,
#'   clip to "on", overlap to "error". Beside plots attach their inner edge
#'   to the true anchor. Map position or 1-based closed start/end (midpoint
#'   `(start - 1 + end) / 2`) and, without plot, a plot
#'   list-column. hjust/vjust can override placement-aware justification.
#' @return An additive positioned-plot component.
#' @export
geom_locus_inset <- function(plot = NULL, mapping = NULL, data = NULL, ...) {
  if (inherits(plot, "uneval") && is.null(mapping)) { mapping <- plot; plot <- NULL }
  if (!is.null(plot)) {
    if (!inherits(plot, "ggplot") && !inherits(plot, "grob")) stopf("`plot` must be a complete ggplot or grob.")
    if (!is.data.frame(data)) stopf("An inset requires data-frame chromosome anchors.")
    column <- ".ideogram_child"
    while (column %in% names(data)) column <- paste0(column, "_")
    data[[column]] <- rep(list(plot), nrow(data))
    mapping$plot <- rlang::new_quosure(rlang::sym(column))
  }
  new_ideogram_object(geom_chr_inset(mapping, data, ...))
}

#' @method ggplot_add ggideogram_object
#' @export
ggplot_add.ggideogram_object <- function(object, plot, ...) {
  initial <- !inherits(plot$coordinates$layout, "ideogram_layout_v2")
  before <- length(plot$layers)
  result <- ggplot2::ggplot_add(object$component, plot,
    ggplot_add_object_name(..., fallback = "ideogram_object"))
  if (initial) return(result)
  layout <- ideogram_plot_layout(result)
  id <- (layout$recipe_count %||% 0L) + 1L
  layout$recipe_count <- id
  if (length(result$layers) > before) for (i in seq.int(before + 1L, length(result$layers))) {
    result$layers[[i]]$ideogram_recipe_id <- id
    result$layers[[i]]$ideogram_recipe <- object$component
    result$layers[[i]] <- record_ideogram_aesthetics(result$layers[[i]])
  }
  update_plot_ideogram_layout(result, layout)
}

#' @method ggplot_add ggideogram_chr_component
#' @export
ggplot_add.ggideogram_chr_component <- function(object, plot, ...) {
  if (!"fill" %in% object$components) {
    visual <- setdiff(names(object$mapping),
      c("chr", "start", "end", "genome", "assembly", "homolog", "label"))
    if (length(visual)) stopf("Chromosome components do not support mapped aesthetics: %s. Use fixed styling or component = 'fill'.",
      paste(visual, collapse = ", "))
  }
  if (!inherits(plot$coordinates$layout, "ideogram_layout_v2")) {
    has_fill <- "fill" %in% object$components
    data <- if (has_fill) plot$data else object$data %||% plot$data
    mapping <- if (has_fill) plot$mapping else object$mapping %||% plot$mapping
    if (!length(mapping)) mapping <- NULL
    args <- c(list(data = data, mapping = mapping, chr = object$chr,
      axis = object$axis, show_names = object$names %||% TRUE), object$params)
    result <- do.call(ggideogram, args)
    if (has_fill) {
      fill <- object
      fill$components <- "fill"; fill$names <- FALSE; fill$axis <- FALSE
      constructor_args <- setdiff(names(formals(ggideogram)),
        c("fill", "colour", "linewidth"))
      fill$params <- fill$params[!names(fill$params) %in% constructor_args]
      result <- result + new_ideogram_object(fill)
    }
    for (layer in plot$layers) result <- result + layer
    result$scales <- plot$scales; result$labels <- plot$labels
    result$guides <- plot$guides; result$facet <- plot$facet
    result$theme <- result$theme + plot$theme
    return(result)
  }
  original <- ideogram_plot_layout(plot)
  selected <- original
  chr <- object$chr
  if (isTRUE(chr)) chr <- NULL
  if (!is.null(object$data) && is.null(chr) && !"fill" %in% object$components) {
    chr <- as_ideogram_data(object$data, object$mapping)$karyotype$.chr
  }
  if (!is.null(chr)) {
    if (any(!chr %in% selected$chrom$.chr)) stopf("Unknown selected chromosome: %s.", format_chr_rows(setdiff(chr, selected$chrom$.chr)))
    selected$chrom <- selected$chrom[match(chr, selected$chrom$.chr), , drop = FALSE]
    selected$index <- stats::setNames(seq_len(nrow(selected$chrom)), selected$chrom$.chr)
  }
  plot <- update_plot_ideogram_layout(plot, selected)
  components <- object$components
  if (isTRUE(object$names) && !"name" %in% components) components <- c(components, "name")
  if (isTRUE(object$axis) && !"axis" %in% components) components <- c(components, "axis")
  accepted <- unique(unlist(lapply(components, chromosome_component_parameters)))
  unused <- setdiff(names(object$params), accepted)
  if (length(unused)) stopf("Unsupported chromosome component parameter%s: %s.",
    if (length(unused) > 1L) "s" else "", paste(unused, collapse = ", "))
  for (kind in components) {
    if (kind == "fill") {
      params <- object$params
      if (length(components) > 1L) {
        accepted <- c(names(formals(geom_chr_fill)), ggplot2::GeomRect$aesthetics(),
          ggplot2::GeomRect$parameters(extra = TRUE), "key_glyph")
        params <- params[intersect(names(params), accepted)]
      }
      plot <- plot + do.call(geom_chr_fill, c(list(mapping = object$mapping,
        data = object$data, track = object$track), params))
      next
    }
    fun <- switch(kind, body = geom_chromosome, name = geom_chr_name,
      axis = geom_chr_axis, band = geom_chr_cytoband)
    accepted <- intersect(names(object$params), names(formals(fun)))
    parameters <- object$params[accepted]
    if (kind == "axis" && !is.null(chr)) parameters$chr <- chr
    plot <- plot + do.call(fun, parameters)
  }
  updated <- ideogram_plot_layout(plot)
  original$bounds <- updated$bounds
  if (!is.null(updated$circular_name_axis)) original$circular_name_axis <- updated$circular_name_axis
  update_plot_ideogram_layout(plot, original)
}

chromosome_component_parameters <- function(kind) {
  fun <- switch(kind, body = geom_chromosome, name = geom_chr_name,
    axis = geom_chr_axis, band = geom_chr_cytoband, fill = geom_chr_fill)
  parameters <- setdiff(names(formals(fun)), "...")
  if (kind == "fill") parameters <- union(parameters,
    c(ggplot2::GeomRect$aesthetics(), ggplot2::GeomRect$parameters(extra = TRUE), "key_glyph"))
  parameters
}
