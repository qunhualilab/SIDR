#' @title Fit the SIDR Model
#'
#' @description
#' This function estimates parameters of the stratified Gaussian copula mixture
#' model using the Nelder–Mead optimization method, and computes the local and
#' global irreproducible discovery rate (idr and IDR) for Hi-C interactions
#' based on the estimated parameters. The computation is inspired by
#' the IDR methodology implemented in the
#' \href{https://bioconductor.org/packages/release/bioc/html/idr.html}{IDR package}
#' (Li et al., 2011).
#'
#' @importFrom stats optim
#' @param hic_data A data frame returned by \code{mergeHiC}, containing merged
#' Hi-C interaction data. Alternatively, it can be a data frame
#' containing well-paired observations for the two replicates (obs1 and obs2).
#' @param ns An integer specifying the number of strata.
#' @param ind A list returned by \code{stratifyData}, where each element contains
#' the indices of the observations belonging to that stratum.
#' @param thet.ini A numeric vector of length \code{3 * ns + 2}.
#' \enumerate{
#'   \item The first \code{ns} elements: logit-transformed initial mixture
#'   proportions of the reproducible component (\eqn{\pi_k});
#'   \item The next \code{ns} elements: initial means of the reproducible
#'   component (\eqn{\mu_k});
#'   \item The next element: initial standard deviation of the reproducible
#'   component (\eqn{\sigma});
#'   \item The next element: logit-transformed initial correlation coefficient
#'   of the irreproducible component (\eqn{\rho});
#'   \item The final \code{ns} elements: initial correlation scaling parameters
#'   of the reproducible component (\eqn{\omega_k}).
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{est}{A numeric vector of estimated model parameters.}
#'   \item{idr}{A list of local irreproducible discovery rates (idr) per stratum.}
#'   \item{IDR}{A numeric vector of global irreproducible discovery rate (IDR).}
#' }
#'
#' @seealso
#' \code{\link{choosePara}} for initializing parameters, and
#' \code{\link{stratifyData}} for defining distance-based strata.
#'
#' @references
#' Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011).
#' Measuring reproducibility of high-throughput experiments.
#' \emph{Annals of Applied Statistics}, 5(3), 1752–1779.
#' \url{https://bioconductor.org/packages/release/bioc/html/idr.html}
#'
#' @examples
#' \dontrun{
#' # Merge Hi-C replicates
#' data <- mergeHiC("rep1.txt", "rep2.txt", chrI = "chr18", chrJ = "chr18")
#'
#' # Stratify data
#' ind <- stratifyData(data, ns = 3)
#'
#' # Decide initial parameters
#' init <- choosePara(data, ns = 3, ind = ind)
#'
#' # Fit the stratified IDR model
#' thet.ini <- c(logit(init$mixps), init$mus, init$sigma, logit(init$rho), init$omega)
#' res <- fit_IDR_stratified(data, ns = 3, ind = ind, thet.ini = thet.ini)
#' }
#'
#' @export
fit_IDR_stratified <- function(hic_data, ns, ind, thet.ini){

  # Extract observations
  y1 <- hic_data$obs1
  y2 <- hic_data$obs2

  # Optimize log-likelihood
  fit <- optim(
    par = thet.ini,
    fn = function(thets) -gmcmLik(y1, y2, ns, ind, thets),
    method = "Nelder-Mead",
    hessian = TRUE
  )

  est <- fit$par
  est.gmcm <- c(
    logit.inv(est[1:ns]),             # mixture proportions
    abs(est[(ns + 1):(2 * ns)]),      # means
    est[(2 * ns + 1)],                # sigma
    logit.inv(est[(2 * ns + 2)]),     # correlation (rho)
    est[(2 * ns + 3):(3 * ns + 2)]    # omega
  )

  mixp.gmcm  <- est.gmcm[1:ns]
  mu.gmcm    <- est.gmcm[(ns + 1):(2 * ns)]
  sigma.gmcm <- est.gmcm[2 * ns + 1]
  rho.gmcm   <- est.gmcm[2 * ns + 2]
  omega.gmcm <- est.gmcm[(2 * ns + 3):(3 * ns + 2)]

  # Compute idr per stratum
  idr.gmcm <- get.idr.gmcm(
    y1, y2, ns, ind, mixp.gmcm,
    mu.gmcm, sigma.gmcm, rho.gmcm, omega.gmcm
  )

  idr.all.gmcm <- numeric(nrow(hic_data))
  for(k in seq_len(ns)){
    idr.all.gmcm[ind[[k]]] <- idr.gmcm$idr[[k]]
  }

  # Compute IDR
  IDR.all.gmcm <- get.IDR(idr.all.gmcm)

  return(list(
    est = est.gmcm,
    idr = idr.gmcm$idr,
    IDR = IDR.all.gmcm
  ))
}

