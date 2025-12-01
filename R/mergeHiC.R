#' @title Merge two Hi-C replicates
#'
#' @description
#' This function preprocesses two Hi-C replicate files of the same resolution,
#' filters them for specified chromosomes, and merges them into a
#' single data frame.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr filter transmute
#' @importFrom stats na.omit
#' @importFrom rlang .data
#'
#' @param hic_rep1 Path to the text file (CSV, TSV, or TXT) of the first Hi-C
#' replicate. The file should be a data frame with the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chromosomeI} \tab character; chromosome name of
#'   fragment I. \cr
#'   column 2: \tab \code{fragmentI} \tab integer; coordinate of the start
#'   point of fragment I. \cr
#'   column 3: \tab \code{chromosomeJ} \tab character; chromosome name of
#'   fragment J. \cr
#'   column 4: \tab \code{fragmentJ} \tab integer; coordinate of the start
#'   point of fragment J. \cr
#'   column 5: \tab \code{pvalue} \tab numeric; p-value measuring the
#'   significance of the interaction between the two fragments. \cr
#' }
#'
#' @param hic_rep2 Path to the text file (CSV, TSV, or TXT) of the second Hi-C
#' replicate. The file should have the same columns as \code{hic_rep1}.
#' @param chrI A character string specifying the chromosome name of fragment I
#' to be analyzed.
#' @param chrJ A character string specifying the chromosome name of fragment J
#' to be analyzed.
#'
#' @return A data frame containing the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chromosomeI} \tab character; chromosome name of
#'   fragment I. \cr
#'   column 2: \tab \code{chromosomeJ} \tab character; chromosome name of
#'   fragment J. \cr
#'   column 3: \tab \code{fragmentI} \tab integer; coordinate of the start
#'   point of fragment I. \cr
#'   column 4: \tab \code{fragmentJ} \tab integer; coordinate of the start
#'   point of fragment J. \cr
#'   column 5: \tab \code{pvalue1} \tab numeric; p-value for the first
#'   replicate. \cr
#'   column 6: \tab \code{pvalue2} \tab numeric; p-value for the second
#'   replicate. \cr
#'   column 7: \tab \code{obs1} \tab numeric; -log(\code{pvalue1}). \cr
#'   column 8: \tab \code{obs2} \tab numeric; -log(\code{pvalue2}). \cr
#'   column 9: \tab \code{dist} \tab numeric; genomic distance between
#'   fragment I and fragment J (kb). \cr
#' }
#'
#' @examples
#' \dontrun{
#' mergeHiC("replicate1.txt", "replicate2.txt", chrI = "chr18", chrJ = "chr18")
#' }
#'
#' @export
mergeHiC <- function(hic_rep1, hic_rep2, chrI, chrJ){

  rep1 <- data.table::fread(hic_rep1)
  rep2 <- data.table::fread(hic_rep2)

  required_cols <-
    c("chromosomeI", "fragmentI", "chromosomeJ", "fragmentJ", "pvalue")
  if (!all(required_cols %in% colnames(rep1))) {
    stop("Replicate 1 does not have full columns as required.")
  }
  if (!all(required_cols %in% colnames(rep2))) {
    stop("Replicate 2 does not have full columns as required.")
  }

  rep1 <- rep1 %>%
    stats::na.omit() %>%
    dplyr::filter(.data$chromosomeI == chrI &
                    .data$chromosomeJ == chrJ)
  rep2 <- rep2 %>%
    stats::na.omit() %>%
    dplyr::filter(.data$chromosomeI == chrI &
                    .data$chromosomeJ == chrJ)

  hic_data <- merge(rep1, rep2, by = c("fragmentI", "fragmentJ")) %>%
    dplyr::mutate(
      obs1 = jitter(-log(pmax(.data$pvalue.x, .Machine$double.eps)), factor = 1e-4),
      obs2 = jitter(-log(pmax(.data$pvalue.y, .Machine$double.eps)), factor = 1e-4),
      dist = abs(.data$fragmentI - .data$fragmentJ) / 1000
    ) %>%
    dplyr::transmute(
      chromosomeI = .data$chromosomeI.x,
      chromosomeJ = .data$chromosomeJ.x,
      fragmentI = .data$fragmentI,
      fragmentJ = .data$fragmentJ,
      pvalue1 = .data$pvalue.x,
      pvalue2 = .data$pvalue.y,
      obs1 = .data$obs1,
      obs2 = .data$obs2,
      dist = .data$dist
    )

  if (nrow(hic_data) == 0) {
    stop("No data available to perform the reproducibility analysis.")
  }

  return(hic_data)
}

