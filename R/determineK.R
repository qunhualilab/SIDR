#' @title Decide the number of strata for Hi-C reproducibility analysis
#'
#' @description
#' Determines the number of distance-defined strata using a data-driven procedure
#' based on rank correlations between Hi-C replicates across genomic distance.
#'
#' @param hic_data A data frame returned by \code{mergeHiC}, containing merged
#' Hi-C interaction data.
#' @param max_ns An integer specifying the maximum number of strata to consider.
#' The default is 10.
#' @param corr_threshold A numeric threshold for the minimum difference in
#' adjacent stratum-specific rank correlations. The default is 0.1.
#'
#' @return An integer giving the selected number of strata.
#'
#' @examples
#' set.seed(1)    # uses random jitter
#' rep1 <- system.file("extdata", "rep1.txt", package = "SIDR")
#' rep2 <- system.file("extdata", "rep2.txt", package = "SIDR")
#'
#' # Merge Hi-C replicates
#' hic_data <- mergeHiC(rep1, rep2, chrI = "chr18", chrJ = "chr18")
#'
#' # Decide the number of strata
#' K <- determineK(hic_data, max_ns = 10, corr_threshold = 0.1)
#'
#' @export
determineK <- function(hic_data, max_ns = 10, corr_threshold = 0.1) {
  corr_diff <- Inf
  ns <- 1

  while (ns < max_ns && all(corr_diff >= corr_threshold)) {
    ns <- ns + 1
    ind <- stratifyData(hic_data, ns = ns)
    corr_ind <- rep(NA, ns)

    for (k in seq_len(ns)) {
      rk1 <- rank(hic_data$obs1[ind[[k]]], ties.method = "random")
      rk2 <- rank(hic_data$obs2[ind[[k]]], ties.method = "random")
      corr_ind[k] <- cor(rk1, rk2)
    }

    corr_diff <- round(abs(diff(corr_ind)), 2)
  }

  K <- ns - 1

  if (K == 1) {
    message("Please use an unstratified model instead.")
  } else {
    return(K)
  }
}
