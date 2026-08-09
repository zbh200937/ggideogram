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

#' Human gene density
#'
#' Gene counts in 1 Mb windows across [human_karyotype]. Suitable for a
#' [geom_track_line()] or [geom_track_col()] layer after mapping each window
#' midpoint to `position`.
#'
#' @format A data frame with 3102 rows and 4 columns: `Chr`, `Start`, `End`,
#'   `Value`.
#' @source RIdeogram: \url{https://github.com/TickingClock1992/RIdeogram}
"gene_density"

#' Human LTR retrotransposon density
#'
#' LTR counts in the same 1 Mb windows as [gene_density], so the two can be
#' drawn as two tracks of the same figure.
#'
#' @format A data frame with 3102 rows and 4 columns: `Chr`, `Start`, `End`,
#'   `Value`.
#' @source RIdeogram: \url{https://github.com/TickingClock1992/RIdeogram}
"LTR_density"

#' 500 random non-coding RNAs
#'
#' Marker positions inherited from RIdeogram. Use the interval midpoint as the
#' `position` aesthetic in [geom_chr_marker()]. `color` retains the source
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
#' @source RIdeogram: \url{https://github.com/TickingClock1992/RIdeogram}
"Random_RNAs_500"
