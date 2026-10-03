#' Identify chromosomes across genomes and assemblies
#'
#' Keys preserve the distinction between equal chromosome names in different
#' genomes. Map the returned vector to `chr` in annotations, tracks, paired
#' connections or an ordinary chromosome position scale.
#'
#' @param genome,chr Character vectors of genome and chromosome identifiers.
#' @param assembly Optional assembly identifiers, included in the key.
#'   Length-one inputs are recycled; other lengths must agree.
#' @return A character vector of unambiguous composite keys.
#' @export
chr_key <- function(genome, chr, assembly = NULL) {
  fields <- list(genome = genome, chr = chr)
  if (!is.null(assembly)) fields$assembly <- assembly
  sizes <- lengths(fields)
  n <- max(sizes)
  if (!n) return(character())
  if (any(!sizes %in% c(1L, n))) {
    stopf("Chromosome key inputs must have equal lengths or length one.")
  }
  fields <- lapply(names(fields), function(name) {
    value <- enc2utf8(validate_chr(fields[[name]], name))
    value <- rep_len(value, n)
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  paste0(if (is.null(assembly)) "gc:" else "gca:",
         do.call(paste0, fields))
}

canonical_genome_fields <- function(karyotype, mapping) {
  optional <- intersect(c("genome", "assembly", "homolog", "label"), names(mapping))
  fields <- mapped_fields(karyotype, mapping, optional, "mapping")
  karyotype$.chr_name <- karyotype$.chr
  for (field in optional) {
    target <- switch(field, homolog = ".homolog_group", label = ".display_name",
                     paste0(".", field))
    karyotype[[target]] <- validate_chr(fields[[field]], paste0("mapping$", field))
  }
  if ("assembly" %in% optional && !"genome" %in% optional) {
    stopf("Mapping `assembly` requires an explicit `genome` mapping.")
  }
  if ("genome" %in% optional) {
    karyotype$.chr <- chr_key(karyotype$.genome, karyotype$.chr_name,
                             karyotype$.assembly)
  }
  if (!"label" %in% optional) {
    karyotype$.display_name <- if ("genome" %in% optional) {
      paste(karyotype$.genome, karyotype$.chr_name, sep = "\n")
    } else karyotype$.chr_name
  }
  karyotype
}

layout_level_order <- function(values, order, what) {
  levels <- unique(as.character(values))
  if (is.null(order)) return(levels)
  order <- validate_chr(order, what)
  if (anyDuplicated(order) || !setequal(order, levels)) {
    stopf("`%s` must list every observed group exactly once.", what)
  }
  order
}

chromosome_layout_slots <- function(k, ncol, order_by, genome_order, homolog_order) {
  n <- nrow(k)
  if (order_by == "input") {
    if (!is.null(genome_order) || !is.null(homolog_order)) {
      stopf("Explicit group orders require `order_by = 'genome'` or 'homolog'.")
    }
    cols <- validate_layout_ncol(ncol, n)
    return(list(index = seq_len(n), row = (seq_len(n) - 1L) %/% cols,
                col = (seq_len(n) - 1L) %% cols, ncol = cols))
  }
  if (is.null(k$.genome)) stopf("Grouped layout requires a `genome` mapping.")
  genomes <- layout_level_order(k$.genome, genome_order, "genome_order")
  group_key <- chr_key(k$.genome, "group", k$.assembly)
  group_levels <- unlist(lapply(genomes, function(genome) {
    unique(group_key[k$.genome == genome])
  }), use.names = FALSE)
  group <- match(group_key, group_levels) - 1L
  has_homolog <- !is.null(k$.homolog_group)
  if (order_by == "homolog" && !has_homolog) {
    stopf("`order_by = 'homolog'` requires an explicit `homolog` mapping.")
  }
  if (!has_homolog && !is.null(homolog_order)) {
    stopf("`homolog_order` requires a `homolog` mapping.")
  }
  if (has_homolog) {
    homologs <- layout_level_order(k$.homolog_group, homolog_order, "homolog_order")
    homolog <- match(k$.homolog_group, homologs) - 1L
    if (order_by == "homolog") {
      cols <- validate_layout_ncol(ncol, length(group_levels))
      row <- homolog * ceiling(length(group_levels) / cols) + group %/% cols
      col <- group %% cols
    } else {
      cols <- validate_layout_ncol(ncol, length(homologs))
      row <- group * ceiling(length(homologs) / cols) + homolog %/% cols
      col <- homolog %% cols
    }
  } else {
    counts <- tabulate(group + 1L, nbins = length(group_levels))
    cols <- validate_layout_ncol(ncol, max(counts))
    starts <- c(0, utils::head(cumsum(ceiling(counts / cols)), -1L))
    rank <- stats::ave(seq_len(n), group, FUN = seq_along) - 1L
    row <- starts[group + 1L] + rank %/% cols
    col <- rank %% cols
  }
  if (anyDuplicated(paste(row, col))) {
    stopf("Each genome/assembly must have at most one chromosome per homolog group.")
  }
  index <- order(row, col)
  list(index = index, row = row[index], col = col[index], ncol = cols)
}
