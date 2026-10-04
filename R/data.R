# Datasets, carried over unchanged from RIdeogram so that its examples and any
# comparison scripts keep working against this package.

#' Human karyotype
#'
#' The 24 human chromosomes with centromere positions, in the five-column form
#' accepted by [ggideogram()].
#'
#' @format A data frame with 24 rows and 5 columns:
#' \describe{
#'   \item{Chr}{Chromosome name.}
#'   \item{Start, End}{Chromosome extent, in base pairs.}
#'   \item{CE_start, CE_end}{Centromere extent, in base pairs.}
#' }
#' @source RIdeogram: \url{https://github.com/TickingClock1992/RIdeogram}
"human_karyotype"

#' Liriodendron karyotype
#'
#' The 19 *Liriodendron chinense* chromosomes, in the three-column form. Without
#' `CE_start`/`CE_end` the chromosomes are drawn as plain capsules.
#'
#' @format A data frame with 19 rows and 3 columns: `Chr`, `Start`, `End`.
#' @source RIdeogram: \url{https://github.com/TickingClock1992/RIdeogram}
"liriodendron_karyotype"

#' Human gene counts by window
#'
#' Gene counts in 1 Mb windows across [human_karyotype]. Suitable for a
#' [geom_track()] with native ggplot2::geom_line() or ggplot2::geom_col() after
#' mapping each window midpoint to x and its count to y.
#' Terminal windows can be shorter than 1 Mb. `Value` is a count, not a
#' width-normalized density; genes per Mb are `Value / ((End - Start + 1) / 1e6)`.
#' The original annotation release is not recorded by the upstream dataset.
#'
#' @format A data frame with 3102 rows and 4 columns: `Chr`, `Start`, `End`,
#'   `Value`.
#' @source GENCODE annotations, distributed with RIdeogram:
#'   \url{https://github.com/TickingClock1992/RIdeogram}.
"gene_density"

#' Human LTR retrotransposon counts by window
#'
#' LTR counts in the same 1 Mb windows as [gene_density], so the two can be
#' drawn as two tracks of the same figure.
#' `Value` is the count in the actual closed window; terminal windows can be
#' shorter than 1 Mb. Divide by `(End - Start + 1) / 1e6` for counts per Mb.
#'
#' @format A data frame with 3102 rows and 4 columns: `Chr`, `Start`, `End`,
#'   `Value`.
#' @source RIdeogram: \url{https://github.com/TickingClock1992/RIdeogram}
"LTR_density"

#' 500 random non-coding RNAs
#'
#' A random sample of 500 tRNA, rRNA and miRNA annotations from GENCODE,
#' inherited unchanged from RIdeogram. Counts of this sample describe sampled
#' records, not the complete RNA annotation catalogue or RNA expression.
#' The original GENCODE release and sampling seed are not recorded upstream.
#' Use the interval midpoint as the
#' `position` aesthetic in [geom_locus()]. `color` retains the source
#' package's bare hexadecimal strings, so add `#` before passing them to a
#' standard ggplot2 scale.
#'
#' @format A data frame with 500 rows and 6 columns:
#' \describe{
#'   \item{Type}{Marker class; one legend row per distinct value.}
#'   \item{Shape}{Source shape name (`circle`, `box`, or `triangle`), which maps
#'     naturally to standard filled ggplot2 shapes 21, 22, and 24.}
#'   \item{Chr}{Chromosome name, matched against the karyotype.}
#'   \item{Start, End}{Marker extent; the anchor is the midpoint.}
#'   \item{color}{Fill colour.}
#' }
#' @source GENCODE annotations, randomly sampled and distributed with RIdeogram:
#'   \url{https://github.com/TickingClock1992/RIdeogram}.
"Random_RNAs_500"
