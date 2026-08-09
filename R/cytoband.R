# Giemsa-stained cytogenetic band data and palettes.

CYTOBAND_SCHEMES <- list(
  circos = c(
    gneg = "#FFFFFF", gpos = "#000000",
    gpos25 = "#C8C8C8", gpos33 = "#D2D2D2", gpos50 = "#C8C8C8",
    gpos66 = "#A0A0A0", gpos75 = "#828282", gpos100 = "#000000",
    gvar = "#DCDCDC", acen = "#D92F27", stalk = "#647FA4"
  ),
  biovizbase = c(
    gneg = "#FFFFFF", gpos = "#000000",
    gvar = "#000000", acen = "#8B2323", stalk = "#CD3333"
  ),
  only.centromeres = c(
    gneg = "#C8C8C8", gpos = "#C8C8C8",
    gvar = "#C8C8C8", acen = "#D92F27", stalk = "#C8C8C8"
  )
)

#' Colours for Giemsa stains
#'
#' Maps cytoband stain names to established palettes. Untabulated `gposNN`
#' names are interpolated between `gneg` and `gpos`.
#'
#' @param stain Character vector such as `"gneg"`, `"gpos50"`, or `"acen"`.
#' @param scheme Palette family.
#' @param palette Optional named colour overrides.
#' @param bleach Proportion of white mixed into each colour, from 0 to 1.
#'
#' @return A hexadecimal colour vector.
#' @export
cytoband_colours <- function(
    stain,
    scheme = c("circos", "biovizbase", "only.centromeres"),
    palette = NULL,
    bleach = 0) {
  colours <- CYTOBAND_SCHEMES[[match.arg(scheme)]]
  if (!is.null(palette)) {
    palette <- unlist(palette)
    if (is.null(names(palette)) || !all(nzchar(names(palette)))) {
      stopf("`palette` must be named, e.g. `c(gpos = \"#08306b\")`.")
    }
    colours[names(palette)] <- normalise_colour(palette)
    if (any(c("gneg", "gpos") %in% names(palette))) {
      tabulated <- grep("^gpos[0-9]+$", names(colours), value = TRUE)
      colours <- colours[
        setdiff(names(colours), setdiff(tabulated, names(palette)))
      ]
    }
  }
  if (!is.numeric(bleach) || length(bleach) != 1L || is.na(bleach) ||
      bleach < 0 || bleach > 1) {
    stopf("`bleach` must be a single number between 0 and 1.")
  }

  stain <- tolower(trimws(as.character(stain)))
  result <- unname(colours[match(stain, names(colours))])
  percentage <- suppressWarnings(as.numeric(sub("^gpos", "", stain)))
  interpolate <- which(is.na(result) & !is.na(percentage))
  if (length(interpolate)) {
    ramp <- grDevices::colorRamp(c(colours[["gneg"]], colours[["gpos"]]))
    rgb <- ramp(pmin(pmax(percentage[interpolate] / 100, 0), 1))
    result[interpolate] <- grDevices::rgb(
      rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255
    )
  }

  unknown <- unique(stain[is.na(result)])
  if (length(unknown)) {
    stopf(paste0(
      "Unknown Giemsa stain%s: %s.\n",
      "  Known names include %s, plus any gposNN. Add custom names with ",
      "`palette`."),
      if (length(unknown) > 1L) "s" else "",
      paste0("`", unknown, "`", collapse = ", "),
      paste(utils::head(sort(names(colours)), 8), collapse = ", ")
    )
  }
  if (bleach > 0) {
    rgb <- grDevices::col2rgb(result) * (1 - bleach) + 255 * bleach
    result <- grDevices::rgb(
      rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255
    )
  }
  result
}

#' Read a UCSC cytoband table
#'
#' @param file A path, URL, or data frame containing the five positional UCSC
#'   columns: chromosome, start, end, band name and stain.
#' @param chr Optional chromosome selection and order.
#'
#' @return A data frame with `Chr`, `Start`, `End`, `Name`, and `Stain`.
#' @export
read_cytoband <- function(file, chr = NULL) {
  bands <- if (is.data.frame(file)) {
    file
  } else {
    utils::read.table(
      file, sep = "\t", header = FALSE, quote = "",
      comment.char = "#", stringsAsFactors = FALSE
    )
  }
  if (ncol(bands) < 5L) {
    stopf("`file` needs 5 cytoBand columns; got %d.", ncol(bands))
  }
  result <- data.frame(
    Chr = as.character(bands[[1]]),
    Start = as.numeric(bands[[2]]) + 1,
    End = as.numeric(bands[[3]]),
    Name = as.character(bands[[4]]),
    Stain = as.character(bands[[5]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(chr)) {
    chr <- as.character(chr)
    missing <- setdiff(chr, result$Chr)
    if (length(missing)) {
      stopf("Not in the cytoband table: %s.", format_chr_rows(missing))
    }
    result <- result[result$Chr %in% chr, , drop = FALSE]
    result <- result[order(match(result$Chr, chr), result$Start), , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}

cytoband_alias <- function(cytoband) {
  if (is.data.frame(cytoband) && !"Stain" %in% names(cytoband) &&
      "gieStain" %in% names(cytoband)) {
    cytoband$Stain <- cytoband$gieStain
  }
  cytoband
}

is_cytoband <- function(x) {
  is.data.frame(x) && any(c("Stain", "gieStain") %in% names(x))
}

#' Derive a karyotype from cytobands
#'
#' Chromosome ends come from the final band and centromeres from `acen` bands.
#'
#' @param cytoband A cytoband data frame.
#' @param chr Optional chromosome selection and order.
#'
#' @return A karyotype data frame accepted by [ggideogram()].
#' @export
cytoband_karyotype <- function(cytoband, chr = NULL) {
  cytoband <- cytoband_alias(cytoband)
  if (!is.data.frame(cytoband)) {
    stopf("`cytoband` must be a data frame.")
  }
  required <- c("Chr", "Start", "End", "Stain")
  missing <- setdiff(required, names(cytoband))
  if (length(missing)) {
    stopf("`cytoband` is missing column%s %s.",
          if (length(missing) > 1L) "s" else "",
          paste0("`", missing, "`", collapse = ", "))
  }
  cytoband$Chr <- as.character(cytoband$Chr)
  keep <- if (is.null(chr)) unique(cytoband$Chr) else as.character(chr)
  unknown <- setdiff(keep, cytoband$Chr)
  if (length(unknown)) {
    stopf("Not in the cytoband table: %s.", format_chr_rows(unknown))
  }

  result <- data.frame(
    Chr = keep,
    Start = 0,
    End = vapply(
      keep,
      function(id) max(cytoband$End[cytoband$Chr == id]),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  centromere <- cytoband[
    tolower(trimws(cytoband$Stain)) == "acen", , drop = FALSE
  ]
  if (nrow(centromere)) {
    result$CE_start <- vapply(keep, function(id) {
      value <- centromere$Start[centromere$Chr == id]
      if (length(value)) min(value) else 0
    }, numeric(1))
    result$CE_end <- vapply(keep, function(id) {
      value <- centromere$End[centromere$Chr == id]
      if (length(value)) max(value) else 0
    }, numeric(1))
  }
  rownames(result) <- NULL
  result
}
