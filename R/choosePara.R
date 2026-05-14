#' @title Decide initial parameters for Hi-C reproducibility analysis
#'
#' @description
#' Compute initial parameters for the Gaussian mixture copula model used
#' in stratified IDR analysis.
#'
#' @importFrom stats sd cor
#' @param hic_data A data frame returned by \code{mergeHiC}, containing merged
#' Hi-C interaction data.
#' @param ns An integer specifying the number of strata.
#' @param ind A list returned by \code{stratifyData}, containing row indices
#' of each stratum.
#' @param prp A numeric value between 0 and 1 specifying the proportion
#' threshold for selecting top-ranked observations (highly likely to be
#' reproducible) to compute correlation.
#' @param rho A numeric value between 0 and 1 specifying the correlation for the
#' bottom-ranked signals (highly likely to be irreproducible).
#'
#' @return A list containing:
#' \describe{
#'   \item{mixps}{A numeric vector of reproducible proportions.}
#'   \item{mus}{A numeric vector of means for the reproducible component in each stratum.}
#'   \item{rho}{Correlation coefficient for the irreproducible component.}
#'   \item{sigma}{Standard deviation (set to 1) for the reproducible component.}
#'   \item{omega}{A numeric vector of correlation scale parameters for the reproducible component.}
#' }
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
#' # Select initial parameters
#' init <- choosePara(hic_data, ns = K, ind = ind, prp = 0.7, rho = 0.15)
#'
#' @export
choosePara <- function(hic_data, ns, ind, prp = 0.7, rho = 0.15) {

  if (!is.numeric(prp) || prp <= 0 || prp >= 1)
    stop("prp must be between 0 and 1.")

  if (!is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("rho must be between 0 and 1.")

  if (!is.list(ind) || length(ind) != ns)
    stop("ind must be a list of length ns.")

  mixps <- (ns:1) / (ns + 1)
  mus <- 1.5 - 0.1 * (0:(ns - 1))

  if (any(mus <= 0)) {
    stop("All values in 'mus' must be positive.")
  }

  y1 <- hic_data$obs1
  y2 <- hic_data$obs2

  omega <- numeric(ns)

  for (k in seq_len(ns)) {
    indk <- ind[[k]]
    nk <- length(indk)
    yk1 <- y1[indk]
    yk2 <- y2[indk]

    # Compute correlation among top ranks
    r1 <- rank(yk1, ties.method = "random")
    r2 <- rank(yk2, ties.method = "random")
    idx <- which(r1 >= prp * max(r1) & r2 > prp * max(r2))
    cov_val <- cor(r1[idx], r2[idx])

    if (cov_val < rho) stop("Please choose a smaller rho.")

    # Initialize omega using logit transformation
    omega[k] <- round(sqrt(logit(cov_val) - logit(rho)), 2)
  }

  return(list(mixps = mixps, mus = mus, sigma = 1, rho = rho, omega = omega))
}


