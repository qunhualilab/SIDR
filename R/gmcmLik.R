#' @title Compute log-likelihood of the SIDR model
#'
#' @description
#' This function computes log-likelihood of all Hi-C observations under the
#' stratified Gaussian copula mixture model. It serves as the objective function
#' optimized in \code{\link{fit_IDR_stratified}} to estimate model parameters.
#'
#' @importFrom stats dnorm
#'
#' @param y1 A numeric vector of Hi-C interaction scores
#' (-log p-values) from replicate 1.
#' @param y2 A numeric vector of the same length as \code{y1} from
#' replicate 2.
#' @param ns An integer specifying the number of strata.
#' @param ind A list of length \code{ns}, where each element contains
#' the row indices corresponding to that stratum.
#' @param thets A numeric vector of length \code{3 * ns + 2}.
#' \enumerate{
#'   \item The first \code{ns} elements: mixture proportions of the reproducible
#'   component (\eqn{\pi_k});
#'   \item The next \code{ns} elements: means of the reproducible
#'   component (\eqn{\mu_k});
#'   \item The next element: standard deviation of the reproducible
#'   component. (\eqn{\sigma});
#'   \item The next element: correlation coefficient of the irreproducible
#'   component (\eqn{\rho});
#'   \item The final \code{ns} elements: correlation scaling parameters of the
#'   reproducible component (\eqn{\omega_k}).
#' }
#'
#' @return A single numeric value specifying the total log-likelihood
#' of all observations.
#'
#' @seealso \code{\link{fit_IDR_stratified}} for parameter estimation
#' and \code{\link{get.idr.gmcm}} for local IDR computation.
#'
#' @export
gmcmLik <- function(y1, y2, ns, ind, thets){

  mixps  <- logit.inv(thets[1:ns])
  mus    <- abs(thets[(ns + 1):(2 * ns)])
  sigma  <- thets[(2 * ns + 1)]
  Sigma  <- array(0, c(2, 2, ns))
  Sigma1 <- array(0, c(2, 2, ns))
  rho    <- logit.inv(thets[(2 * ns + 2)])
  omega  <- thets[(2 * ns + 3):(3 * ns + 2)]

  # Construct covariance matrices
  for (k in seq_len(ns)){
    Sigma[,,k]  <- rbind(c(1, rho), c(rho, 1))
    adj_rho <- logit.inv(logit(rho) + omega[k]^2 + 1e-4)
    Sigma1[,,k] <- rbind(c(sigma^2, sigma^2 * adj_rho),
                         c(sigma^2 * adj_rho, sigma^2))
  }

  # Pseudo data
  n <- length(y1)
  y1.rank <- y2.rank <- numeric(n)

  # Marginal CDFs estimate by empirical CDFs
  for (k in seq_len(ns)) {
    indk <- ind[[k]]
    nk <- length(indk)
    y1.rank[indk] <- rank(y1[indk], ties.method = "random") / (nk + 1)
    y2.rank[indk] <- rank(y2[indk], ties.method = "random") / (nk + 1)
  }

  y1.pseudo <- get.pseudo(y1.rank, ns, ind, mixps, mus, sigma)
  y2.pseudo <- get.pseudo(y2.rank, ns, ind, mixps, mus, sigma)
  uu <- cbind(y1.pseudo, y2.pseudo)

  pdf2d <- pdf1 <- pdf2 <- rep(0.1, n)

  for(k in seq_len(ns)){
    indk <- ind[[k]]
    pdf2d[indk] <- (1 - mixps[k]) * dmvnorm(uu[indk, ], mean = c(0, 0),
                                            sigma = Sigma[,,k]) +
      mixps[k] * dmvnorm(uu[indk, ], mean = c(mus[k], mus[k]),
                         sigma = Sigma1[,,k])
    pdf1[indk] <- (1 - mixps[k]) * dnorm(y1.pseudo[indk], mean = 0, sd = 1) +
      mixps[k] * dnorm(y1.pseudo[indk], mean = mus[k], sd = sigma)
    pdf2[indk] <- (1 - mixps[k]) * dnorm(y2.pseudo[indk], mean = 0, sd = 1) +
      mixps[k] * dnorm(y2.pseudo[indk], mean = mus[k], sd = sigma)
  }

  lik <- sum(log(pdf2d)) - sum(log(pdf1) + log(pdf2))
  return(lik)
}

