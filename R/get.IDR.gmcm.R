#' @title Compute idr for the SIDR Model
#'
#' @description
#' This function computes the local irreproducible discovery rate (idr)
#' for each stratum under the fitted stratified Gaussian copula mixture model.
#'
#' @importFrom stats ecdf
#' @param y1 A numeric vector of Hi-C interaction scores
#' (-log p-values) from replicate 1.
#' @param y2 A numeric vector of the same length as \code{y1} from
#' replicate 2.
#' @param ns An integer specifying the number of strata.
#' @param ind A list of length \code{ns}, where each element contains
#' the row indices corresponding to that stratum.
#' @param mixps A numeric vector of length \code{ns}, the estimated mixture
#' proportions for the reproducible component in each stratum.
#' @param mus A numeric vector of length \code{ns}, the estimated means
#' for the reproducible component for each stratum.
#' @param sigma A numeric value specifying the estimated standard deviation
#' for the reproducible component.
#' @param rho A numeric value specifying the estimated correlation for the
#' irreproducible component.
#' @param omega A numeric vector of length \code{ns}, the estimated correlation
#' scale parameters for the reproducible component.
#'
#' @return idr A list of length \code{ns}, each element containing the local
#' IDR values for observations in that stratum.
#'
#' @export
get.idr.gmcm <- function(y1, y2, ns, ind, mixps, mus, sigma, rho, omega){

  # Construct covariance matrices
  Sigma  <- array(0, c(2, 2, ns))
  Sigma1 <- array(0, c(2, 2, ns))
  for (k in seq_len(ns)){
    Sigma[,,k]  <- rbind(c(1, rho), c(rho, 1))
    adj_rho <- logit.inv(logit(rho) + omega[k]^2 + 1e-4)
    Sigma1[,,k] <- rbind(
      c(sigma^2, sigma^2 * adj_rho),
      c(sigma^2 * adj_rho, sigma^2)
    )
  }

  idr <-  NULL
  # Compute local IDR per stratum
  for(k in seq_len(ns)){
    indk <- ind[[k]]
    yk1 <- y1[indk]
    yk2 <- y2[indk]

    yk1.cdf.func <- ecdf(yk1)
    yk2.cdf.func <- ecdf(yk2)
    afactor <- length(yk1)/(length(yk1) + 1)
    yk1.cdf <- yk1.cdf.func(yk1) * afactor
    yk2.cdf <- yk2.cdf.func(yk2) * afactor

    yk1.pseudo <- get.pseudo.mix(yk1.cdf, mixps[k], mus[k], sigma)
    yk2.pseudo <- get.pseudo.mix(yk2.cdf, mixps[k], mus[k], sigma)
    yk.pseudo <- cbind(yk1.pseudo, yk2.pseudo)

    pdf0 <- (1 - mixps[k]) * dmvnorm(yk.pseudo, mean = c(0, 0), sigma = Sigma[,,k])
    pdf1 <- mixps[k] * dmvnorm(yk.pseudo, mean = c(mus[k], mus[k]), sigma = Sigma1[,,k])

    idr[[k]] <- pdf0 / (pdf0 + pdf1)
  }

  return(list(idr = idr))
}

#' @title Compute global IDR
#'
#' @description
#' Computes the global irreproducible discovery rate (IDR) for all
#' observations based on local IDR values.
#'
#' @param idr A numeric vector of local IDR values.
#'
#' @return A numeric vector of global IDR values for all observations.
#'
#' @export
get.IDR <- function(idr){
  o <- order(idr)
  idr.ordered <- idr[o]
  n.sel <- seq_along(idr)
  IDR.o <- cumsum(idr.ordered) / n.sel
  IDR <- numeric(length(idr))
  IDR[o] <- IDR.o
  return(IDR)
}

