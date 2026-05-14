#' @title Stratify Hi-C data into distance-based strata
#'
#' @description
#' Divides Hi-C interactions into distance-based strata to
#' account for distance-dependent changes in reproducibility.
#' It returns a list of row indices for each stratum.
#'
#' @param hic_data A data frame returned by \code{mergeHiC}, containing merged
#' Hi-C interaction data.
#' @param ns An integer specifying the number of strata into which Hi-C
#' interactions are stratified (must be >= 2).
#'
#' @return A list of length \code{ns}, where each element contains the row
#' indices of that stratum.
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
#' # Stratify data
#' ind <- stratifyData(hic_data, ns = K)
#'
#' @export
stratifyData <- function(hic_data, ns){

  if (!all(c("obs1", "obs2", "dist") %in% names(hic_data))) {
    stop("The input data must contain columns: obs1, obs2, and dist.")
  }
  if (ns < 2) stop("ns must be at least 2.")

  r1 <- rank(hic_data$obs1, ties.method = "random")
  r2 <- rank(hic_data$obs2, ties.method = "random")

  # Compute correlation by distance
  corr <- tapply(seq_len(nrow(hic_data)), hic_data$dist, function(idx) {
    if (length(idx) > 1) cor(r1[idx], r2[idx])
    else NA
  })

  # Determine thresholds by correlation intervals
  corr <- corr[!is.na(corr)]
  corr_min <- min(corr)
  corr_max <- max(corr)
  breaks <- seq(corr_max, corr_min, length.out = ns + 1)

  dthres <- sapply(2:ns, function(k) {
    idx <- which(corr <= breaks[k])
    as.numeric(names(idx)[1])
  })

  # Create index list
  indx <- cut(hic_data$dist, breaks = c(-Inf, dthres, Inf), labels = 1:ns)
  ind <- NULL
  for (k in seq_len(ns)){
    ind[[k]] <- (1:dim(hic_data)[1])[indx == k]
  }

  return(ind)
}





