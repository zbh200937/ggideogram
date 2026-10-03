#' Declare source nodes for chromosome connections
#'
#' Nodes have unique identifiers and retain every input column. A point node
#' maps `position` in original bp. An interval node maps one-based closed
#' `start` and `end`; its point anchor is `(start - 1 + end) / 2`, while the
#' complete interval remains available to [geom_chrlink()] with
#' `type = "interval"`. Different
#' identifiers may occupy the same position.
#'
#' The full source table can be used in local views. Source chromosome names
#' and bounds are checked when a connection is added to an ideogram, before
#' selecting or clipping visible pairs. Use [chr_key()] explicitly for nodes
#' belonging to multiple genomes or assemblies.
#'
#' @param data A source node data frame.
#' @param mapping An [ggplot2::aes()] mapping with `id`, `chr`, and either
#'   `position` or both `start` and `end`.
#' @return A data frame inheriting from `chr_nodes`, with original columns and
#'   canonical `.node_id`, `.node_chr`, `.node_position` fields. Interval nodes
#'   also have `.node_start` and `.node_end` in original closed coordinates.
#' @examples
#' genes <- data.frame(Gene = c("g1", "g2"), Chr = "1",
#'                     Start = c(11, 61), End = c(20, 80))
#' nodes <- chr_nodes(genes, ggplot2::aes(id = Gene, chr = Chr,
#'                                      start = Start, end = End))
#' edges <- data.frame(From = "g1", To = "g2")
#' ggideogram(data.frame(Chr = "1", Start = 0, End = 100)) +
#'   geom_chrlink(data = edges, nodes = nodes,
#'                ggplot2::aes(from = From, to = To))
#' @export
chr_nodes <- function(data, mapping) {
  if (!is.data.frame(data)) stopf('`data` must be a node data frame.')
  interval <- any(c('start', 'end') %in% names(mapping))
  if (interval && 'position' %in% names(mapping)) {
    stopf('Node mapping must use either `position` or `start` and `end`.')
  }
  fields <- mapped_fields(data, mapping,
    c('id', 'chr', if (interval) c('start', 'end') else 'position'), 'mapping')
  reserved <- c('.node_id', '.node_chr', '.node_position', '.node_start', '.node_end')
  if (any(reserved %in% names(data))) stopf('Source data contain reserved `.node_*` fields.')
  id <- validate_node_id(fields$id, 'mapping$id')
  if (anyDuplicated(id)) {
    stopf('Node identifiers must be unique; repeated: %s.',
      format_chr_rows(unique(id[duplicated(id)])))
  }
  data$.node_id <- id
  data$.node_chr <- validate_chr(fields$chr, 'mapping$chr')
  if (interval) {
    start <- validate_coordinate(fields$start, 'mapping$start')
    end <- validate_coordinate(fields$end, 'mapping$end')
    if (any(start < 1 | end < start | start != floor(start) | end != floor(end))) {
      stopf('Interval nodes require integer one-based closed coordinates with 1 <= start <= end.')
    }
    data$.node_start <- start
    data$.node_end <- end
    data$.node_position <- (start - 1 + end) / 2
  } else {
    data$.node_position <- validate_coordinate(fields$position, 'mapping$position')
    if (any(data$.node_position < 0)) stopf('Node positions must be non-negative original bp.')
  }
  class(data) <- c('chr_nodes', class(data))
  data
}

validate_node_id <- function(x, what) {
  x <- as.character(x)
  if (anyNA(x) || any(!nzchar(x))) stopf('`%s` must contain non-missing, non-empty node identifiers.', what)
  x
}

validate_chr_nodes_source <- function(nodes, layout) {
  if (anyDuplicated(nodes$.node_id)) stopf('Node identifiers must be unique.')
  source <- layout$data$karyotype
  i <- match(nodes$.node_chr, source$.chr)
  if (anyNA(i)) {
    stopf('Nodes contain unknown source chromosomes: %s.',
      format_chr_rows(unique(nodes$.node_chr[is.na(i)])))
  }
  interval <- all(c('.node_start', '.node_end') %in% names(nodes))
  low <- if (interval) nodes$.node_start - 1 else nodes$.node_position
  high <- if (interval) nodes$.node_end else nodes$.node_position
  bad <- low < source$.start[i] | high > source$.end[i]
  if (any(bad)) stopf('Node coordinates are outside source bounds for %s.',
    format_chr_rows(nodes$.node_id[bad]))
  invisible(nodes)
}

resolve_chr_pair_nodes <- function(object, layout) {
  nodes <- object$nodes
  if (is.null(nodes)) return(object)
  if (!inherits(nodes, 'chr_nodes')) stopf('`nodes` must be created with `chr_nodes()`.')
  validate_chr_nodes_source(nodes, layout)
  interval <- all(c('.node_start', '.node_end') %in% names(nodes))
  if (object$type == 'interval' && !interval) {
    stopf('Synteny requires interval nodes mapping `start` and `end`.')
  }
  d <- object$data
  if (!is.data.frame(d)) stopf('Chromosome pairs require data-frame `data`.')
  fields <- mapped_fields(d, object$mapping, c('from', 'to'), 'mapping')
  pair_aes <- c('chr1', 'position1', 'start1', 'end1', 'chr2', 'position2', 'start2', 'end2')
  if (any(pair_aes %in% names(object$mapping))) {
    stopf('With `nodes`, map `from` and `to` instead of explicit endpoint coordinates.')
  }
  source_names <- setdiff(names(nodes), c('.node_id', '.node_chr', '.node_position', '.node_start', '.node_end'))
  for (j in 1:2) {
    field <- c('from', 'to')[j]
    id <- validate_node_id(fields[[field]], field)
    at <- match(id, nodes$.node_id)
    if (anyNA(at)) stopf('Unknown node identifiers in `%s`: %s.', field,
      format_chr_rows(unique(id[is.na(at)])))
    extra <- stats::setNames(lapply(source_names, function(name) nodes[[name]][at]),
      paste0('.node', j, '_', source_names))
    canonical <- c('id', 'chr', 'position', if (interval) c('start', 'end'))
    extra <- c(extra, stats::setNames(lapply(canonical, function(name) nodes[[paste0('.node_', name)]][at]),
      paste0('.node_', canonical, j)))
    if (any(names(extra) %in% names(d))) stopf('Edge data contain reserved node endpoint fields.')
    d[names(extra)] <- extra
  }
  endpoint <- if (object$type == 'point') c('chr', 'position') else c('chr', 'start', 'end')
  mapping <- track_mapping_without(object$mapping, c('from', 'to'))
  for (j in 1:2) for (field in endpoint) {
    mapping[[paste0(field, j)]] <- rlang::new_quosure(
      rlang::sym(paste0('.node_', field, j)), env = rlang::empty_env())
  }
  object$data <- d
  object$mapping <- mapping
  object
}
