#' @title Simulate stratified observations under the SIDR model
#'
#' @description
#' This function simulates paired observations under a stratified Gaussian
#' mixture copula model.Each stratum has its own reproducible mixture proportion,
#' mean, and correlation scaling.
#'
#' @param n Integer. Total number of observations to simulate.
#' @param ns Integer. Number of strata.
#' @param ind A list of length \code{ns}, where each element contains indices
#' corresponding to one stratum.
#' @param mixps Numeric vector of length \code{ns}, mixture proportions of the
#' reproducible component (\eqn{\pi_k}).
#' @param mus Numeric vector of length \code{ns}, means of the reproducible
#' component (\eqn{\mu_k}).
#' @param sigma Numeric. Standard deviation of the reproducible
#' component (\eqn{\sigma}).
#' @param rho Numeric. Correlation of the irreproducible component.
#' @param omega Numeric vector of length \code{ns}, correlation scaling
#' parameters (\eqn{\omega_k}) for the reproducible component.
#'
#' @return
#' A data frame with \code{n} rows and three columns:
#' \describe{
#'   \item{obs1}{Simulated observation from replicate 1.}
#'   \item{obs2}{Simulated observation from replicate 2.}
#'   \item{genuine}{Binary indicator (1 = reproducible / genuine signal, 0 = noise).}
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm
#' @importFrom mvtnorm dmvnorm
#' @export
#'
#' @examples
#' set.seed(456)
#' n <- 10000
#' ns <- 3
#' x <- rexp(n, 1/4)
#' indx <-  1 * (x < 1.62) + 2 * (1.62 <= x & x < 3.89) + 3 * (x >= 3.89)
#' ind <- NULL
#' for (k in seq_len(ns)){
#'   ind[[k]] <- (1:n)[indx == k]
#' }
#' mixps <- c(0.8, 0.4, 0.2)
#' mus <- c(1.1, 1.0, 0.9)
#' sigma <- 1
#' rho <- 0.2
#' omega <- c(1.16, 0.75, 0.56)
#' sim_data <- simulateSIDR(n, ns, ind, mixps, mus, sigma, rho, omega)
#' thet.ini <- c(logit(mixps), mus, sigma, logit(rho), omega)
#' result <- fit_IDR_stratified(sim_data, ns, ind, thet.ini)
#'
#' @export
simulateSIDR <- function(n, ns, ind, mixps, mus, sigma, rho, omega){
  # Validate input lengths
  if (length(mixps) != ns || length(mus) != ns || length(omega) != ns) {
    stop("mixps, mus, and omega must each be of length ns.")
  }

  # Covariance matrices for each stratum
  Sigma  <- array(0, c(2, 2, ns))
  Sigma1 <- array(0, c(2, 2, ns))

  for (k in seq_len(ns)){
    Sigma[,,k]  <- rbind(c(1, rho), c(rho, 1))
    adj_rho <- logit.inv(logit(rho) + omega[k]^2 + 1e-4)
    Sigma1[,,k] <- rbind(
      c(sigma^2, sigma^2 * adj_rho),
      c(sigma^2 * adj_rho, sigma^2))
  }

  y <- matrix(0, n, ncol = 2)
  genuine <- rep(0, n)

  # Simulate within each stratum
  for (k in seq_len(ns)){
    indk <- ind[[k]]
    nk <- length(indk)
    Kk <- sample(0:1, size = nk, replace = TRUE, prob = c(1-mixps[k], mixps[k]))  # 1 = genuine
    z <- (Kk == 0) * mvrnorm(nk, c(0, 0), Sigma[,,k]) +
      (Kk == 1) * mvrnorm(nk, c(mus[k], mus[k]), Sigma1[,,k])
    genuine[indk] <- Kk
    # Convert to -log(p) scale
    z <- -log(1 - pnorm(z, mean = 0, sd = 1))
    y[indk, ] <- z
  }

  dat <- data.frame(obs1 = y[,1], obs2 = y[,2], genuine = genuine)
  return(dat)
}
