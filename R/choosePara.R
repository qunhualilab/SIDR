#' @title Decide initial parameters for Hi-C reproducibility analysis.
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
#'   \item{mus}{A numeric vector of mean of the reproducible group in each stratum.}
#'   \item{rho}{Correlation coefficient for the irreproducible group.}
#'   \item{sigma}{Standard deviation (set to 1).}
#'   \item{omega}{A numeric vector of omega parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Merge Hi-C replicates
#' data <- mergeHiC("rep1.txt", "rep2.txt", chrI = "chr18", chrJ = "chr18")
#'
#' # Stratify data
#' stratified_indices <- stratifyData(hic_data, ns = 3)
#'
#' # Decide initial parameters
#' params <- choosePara(data, ns = 3, ind = ind, prp = 0.6, rho = 0.15)
#' }
#'
#' @export
choosePara <- function(hic_data, ns, ind, prp = 0.6, rho = 0.15) {

  if (!ns %in% c(3, 4, 5))
    stop("ns must be 3, 4, or 5.")

  if (!is.numeric(prp) || prp <= 0 || prp >= 1)
    stop("prp must be between 0 and 1.")

  if (!is.numeric(rho) || rho <= 0 || rho >= 1)
    stop("rho must be between 0 and 1.")

  if (!is.list(ind) || length(ind) != ns)
    stop("ind must be a list of length ns.")

  if (ns == 3) {
    mixps <- c(0.70, 0.50, 0.30)
    mus   <- c(0.80, 0.70, 0.60)
  } else if (ns == 4) {
    mixps <- c(0.70, 0.60, 0.30, 0.10)
    mus   <- c(0.80, 0.70, 0.60, 0.50)
  } else if (ns == 5) {
    mixps <- c(0.75, 0.60, 0.45, 0.30, 0.15)
    mus   <- c(0.80, 0.70, 0.60, 0.50, 0.40)
  } else {
    stop("ns must be 3, 4, or 5.")
  }

  y1 <- hic_data$obs1
  y2 <- hic_data$obs2
  p0 <- 1 - mixps

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

