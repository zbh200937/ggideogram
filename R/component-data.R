# Explicit semantic data attachments for ideogram plots.

#' Attach chromosome metadata to an ideogram
#'
#' `chr_data()` declares an explicit chromosome-key connection. Unique
#' metadata columns are joined into the semantic karyotype and layout tables;
#' the original attachment is also kept in a named registry for third-party
#' extensions. Duplicate chromosome keys are rejected unless
#' `relationship = "many-to-many"`, in which case the data is registered but
#' not flattened into one chromosome row.
#'
#' @param data A data frame to attach.
#' @param by A named character vector mapping canonical `chr` to the source
#'   column, for example `c(chr = "Chr")`.
#' @param relationship Expected key relationship. `"many-to-one"` permits a
#'   subset of unique chromosome keys; `"one-to-one"` additionally requires
#'   every chromosome exactly once; `"many-to-many"` must be explicit and is
#'   stored without flattening.
#' @param name Optional stable attachment name. If `NULL`, the expression used
#'   to add the component is used.
#'
#' @return An additive ideogram data component.
#' @export
chr_data <- function(
    data,
    by = c(chr = "Chr"),
    relationship = c("many-to-one", "one-to-one", "many-to-many"),
    name = NULL) {
  relationship <- match.arg(relationship)
  validate_attachment_input(data, by, required = "chr", name = name)
  new_ideogram_data_component("chr", data, by, relationship, name)
}

#' Attach locus data to an ideogram
#'
#' `locus_data()` validates chromosome and genomic-position keys and stores a
#' canonical `.chr`/`.position` table in the plot's attachment registry. It
#' does not guess how an arbitrary downstream geom should use the records;
#' extensions can retrieve them with [get_ideogram_attachment()] and use the
#' public projection functions.
#'
#' @param data A data frame to attach.
#' @param by A named character vector mapping canonical `chr` and `position` to
#'   source columns, for example `c(chr = "Chr", position = "Pos")`.
#' @param relationship Expected key relationship. Repeated chromosome-position
#'   pairs require the explicit `"many-to-many"` choice.
#' @param name Optional stable attachment name.
#'
#' @return An additive ideogram data component.
#' @export
locus_data <- function(
    data,
    by = c(chr = "Chr", position = "Pos"),
    relationship = c("one-to-one", "many-to-many"),
    name = NULL) {
  relationship <- match.arg(relationship)
  validate_attachment_input(
    data, by, required = c("chr", "position"), name = name)
  new_ideogram_data_component("locus", data, by, relationship, name)
}

validate_attachment_input <- function(data, by, required, name) {
  if (!is.data.frame(data)) {
    stopf("Attachment `data` must be a data frame, not %s.", class(data)[1])
  }
  if (!nrow(data)) stopf("Attachment `data` has no rows.")
  if (!is.character(by) || is.null(names(by)) ||
      !identical(sort(names(by)), sort(required)) ||
      anyNA(by) || any(!nzchar(by))) {
    stopf("`by` must map exactly %s to source column names.",
          paste0("`", required, "`", collapse = " and "))
  }
  missing <- setdiff(unname(by), names(data))
  if (length(missing)) {
    stopf("Attachment `data` has no key column%s %s.",
          if (length(missing) > 1L) "s" else "",
          paste0("`", missing, "`", collapse = ", "))
  }
  if (!is.null(name) &&
      (!is.character(name) || length(name) != 1L ||
       is.na(name) || !nzchar(name))) {
    stopf("Attachment `name` must be NULL or one non-empty string.")
  }
  invisible(TRUE)
}

new_ideogram_data_component <- function(kind, data, by, relationship, name) {
  structure(
    list(
      kind = kind, data = data, by = by,
      relationship = relationship, name = name
    ),
    class = "ggideogram_data_component"
  )
}

#' @method ggplot_add ggideogram_data_component
#' @importFrom ggplot2 ggplot_add
#' @export
ggplot_add.ggideogram_data_component <- function(object, plot, ...) {
  layout <- ideogram_plot_layout(plot)
  object_name <- ggplot_add_object_name(..., fallback = "attachment")
  attachment_name <- object$name %||% object_name
  if (!is.character(attachment_name) || length(attachment_name) != 1L ||
      is.na(attachment_name) || !nzchar(attachment_name)) {
    stopf("Could not derive a stable attachment name; supply `name`.")
  }
  layout$attachments <- layout$attachments %||%
    list(chr = list(), locus = list())
  if (!is.null(layout$attachments[[object$kind]][[attachment_name]])) {
    stopf("A `%s` attachment named `%s` already exists.",
          object$kind, attachment_name)
  }

  layout <- if (object$kind == "chr") {
    attach_chr_component(layout, object, attachment_name)
  } else {
    attach_locus_component(layout, object, attachment_name)
  }
  update_plot_ideogram_layout(plot, layout)
}

attach_chr_component <- function(layout, object, name) {
  chr <- as.character(object$data[[object$by[["chr"]]]])
  validate_attachment_chr(layout, chr)
  duplicated_key <- duplicated(chr) | duplicated(chr, fromLast = TRUE)
  if (any(duplicated_key) && object$relationship != "many-to-many") {
    stopf(paste0(
      "Chromosome attachment `%s` has duplicate key%s %s.\n",
      "  Use `relationship = \"many-to-many\"` only when duplication is ",
      "intentional."), name, if (length(unique(chr[duplicated_key])) > 1L) "s" else "",
      format_chr_rows(unique(chr[duplicated_key])))
  }
  if (object$relationship == "one-to-one" &&
      !setequal(chr, layout$chrom$.chr)) {
    missing <- setdiff(layout$chrom$.chr, chr)
    stopf("One-to-one chromosome attachment `%s` is missing %s.",
          name, format_chr_rows(missing))
  }

  attached <- object$data
  attached$.chr <- chr
  layout$attachments$chr[[name]] <- attached
  if (object$relationship == "many-to-many") return(layout)

  source_columns <- setdiff(names(object$data), object$by[["chr"]])
  collisions <- intersect(
    source_columns,
    union(names(layout$data$karyotype), names(layout$chrom)))
  if (length(collisions)) {
    stopf(paste0(
      "Chromosome attachment `%s` would overwrite column%s %s.\n",
      "  Rename the metadata columns before attaching."),
      name, if (length(collisions) > 1L) "s" else "",
      paste0("`", collisions, "`", collapse = ", "))
  }
  match_index <- match(layout$chrom$.chr, chr)
  for (column in source_columns) {
    value <- object$data[[column]][match_index]
    layout$data$karyotype[[column]] <- value
    layout$chrom[[column]] <- value
  }
  layout
}

attach_locus_component <- function(layout, object, name) {
  chr <- as.character(object$data[[object$by[["chr"]]]])
  position <- object$data[[object$by[["position"]]]]
  validate_attachment_chr(layout, chr)
  if (!is.numeric(position) || anyNA(position) ||
      any(!is.finite(position))) {
    stopf("Locus attachment positions must be finite numeric coordinates.")
  }
  project_positions_checked(layout, chr, position)
  key <- paste(chr, format(position, digits = 17, scientific = FALSE),
               sep = "\r")
  duplicated_key <- duplicated(key) | duplicated(key, fromLast = TRUE)
  if (any(duplicated_key) && object$relationship != "many-to-many") {
    detail <- unique(sprintf("%s:%s", chr[duplicated_key],
                             format(position[duplicated_key], trim = TRUE)))
    stopf(paste0(
      "Locus attachment `%s` has duplicate key%s %s.\n",
      "  Use `relationship = \"many-to-many\"` only when duplication is ",
      "intentional."), name, if (length(detail) > 1L) "s" else "",
      paste0("`", detail, "`", collapse = ", "))
  }
  attached <- object$data
  attached$.chr <- chr
  attached$.position <- as.numeric(position)
  layout$attachments$locus[[name]] <- attached
  layout
}

validate_attachment_chr <- function(layout, chr) {
  if (anyNA(chr) || any(!nzchar(chr))) {
    stopf("Attachment chromosome keys must be non-missing and non-empty.")
  }
  unknown <- setdiff(unique(chr), layout$chrom$.chr)
  if (length(unknown)) {
    stopf("Attachment contains unknown chromosome%s %s.",
          if (length(unknown) > 1L) "s" else "",
          format_chr_rows(unknown))
  }
  invisible(chr)
}

#' Retrieve data attached to an ideogram plot
#'
#' This is the stable read boundary for extensions; callers do not need to
#' inspect the ggplot object or coordinate internals.
#'
#' @param plot A [ggideogram()] plot.
#' @param name Attachment name supplied to [chr_data()] or [locus_data()].
#' @param type Attachment registry, `"chr"` or `"locus"`.
#'
#' @return The canonical attached data frame.
#' @export
get_ideogram_attachment <- function(plot, name,
                                     type = c("chr", "locus")) {
  type <- match.arg(type)
  layout <- ideogram_plot_layout(plot)
  data <- layout$attachments[[type]][[name]]
  if (is.null(data)) {
    available <- names(layout$attachments[[type]])
    stopf("No `%s` attachment named `%s`; available: %s.",
          type, name,
          if (length(available)) paste(available, collapse = ", ") else "none")
  }
  data
}
