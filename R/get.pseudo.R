#' @title Compute pseudo-data per stratum
#'
#' @description
#' This function computes pseudo-observations corresponding to given CDF values
#' within each stratum. It is adapted from the \pkg{idr} package (Li et al., 2011),
#' and serves as a stratified extension of \code{\link{get.pseudo.mix}}.
#'
#' @param u A numeric vector of CDF values (between 0 and 1),
#' typically obtained from empirical ranks.
#' @param ns An integer specifying the number of strata.
#' @param ind A list of length \code{ns}, where each element contains
#' the row indices corresponding to that stratum.
#' @param mixps Numeric vector of length \code{ns},
#' mixture proportions (\eqn{\pi_k}) for each stratum.
#' @param mus Numeric vector of length \code{ns}, means (\eqn{\mu_k}) of
#' the reproducible component for each stratum.
#' @param sigma A numeric value, the standard deviation (\eqn{\sigma}) of
#' the reproducible component.
#'
#' @return A numeric vector of pseudo-observations,
#' with the same length and order as \code{u}.
#'
#' @references
#' Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011).
#' Measuring reproducibility of high-throughput experiments.
#' *The Annals of Applied Statistics*, 5(3), 1752–1779.
#' See also the \pkg{idr} package: <https://cran.r-project.org/package=idr>
#'
#' @seealso \code{\link{get.pseudo.mix}} for the mixture-based inverse transformation.
#'
#' @export
get.pseudo <- function(u, ns, ind, mixps, mus, sigma){
  u.qnt <- numeric(length(u))
  for(k in seq_len(ns)){
    indk <- ind[[k]]
    u.qnt[indk] <- get.pseudo.mix(u[indk], mixps[k], mus[k], sigma)
  }

  invisible(u.qnt)
}

#' @title The inverse of marginal distribution by package idr
#'
#' @description
#' This function computes pseudo-observations corresponding to given
#' CDF values (\code{z}) under a two-component Gaussian mixture model.
#' It is adapted from the \pkg{idr} package (Li et al., 2011).
#'
#' @param z A numeric vector of CDF values (between 0 and 1).
#' @param mixp Proportion of the reproducible component.
#' @param mu Mean of the reproducible component.
#' @param sigma Standard deviation of the reproducible component.
#' @references
#' Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011).
#' Measuring reproducibility of high-throughput experiments.
#' *The Annals of Applied Statistics*, 5(3), 1752–1779.
#' See also the \pkg{idr} package: <https://cran.r-project.org/package=idr>
#'
#' @export
get.pseudo.mix <- function (z, mixp, mu, sigma){
  nw <- 1000
  w <- seq(min(-3, mu - 3 * sigma), max(mu + 3 * sigma, 3), length = nw)
  w.cdf <- (1 - mixp) * pnorm(w, mean = 0, sd = 1) +
    mixp * pnorm(w, mean = mu, sd = sigma)

  quan.z <- rep(NA, length(z))
  for (ii in seq_len(nw)) {
    index <- which(z >= w.cdf[ii] & z < w.cdf[ii + 1])
    quan.z[index] <- (z[index] - w.cdf[ii]) *
      (w[ii + 1] - w[ii])/(w.cdf[ii + 1] - w.cdf[ii]) + w[ii]
  }
  index <- which(z < w.cdf[1])
  if (length(index) > 0)
    quan.z[index] <- w[1]
  index <- which(z > w.cdf[nw])
  if (length(index) > 0)
    quan.z[index] <- w[nw]
  invisible(quan.z)
}

