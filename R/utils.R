#' @title Logit, Inverse Logit Transformations,
#' and Derivative of Inverse Logit
#'
#' @description
#' Utility functions to perform the logit and inverse logit transformations,
#' as well as compute the derivative of the inverse logit function.
#' These are used internally in the SIDR model to ensure that
#' parameters such as correlation and mixture proportions remain within
#' the valid interval (0, 1).
#'
#' @param p A numeric vector of probabilities (for \code{logit}),
#' with all elements strictly between 0 and 1.
#' @param x A numeric vector of real values (for \code{logit.inv})
#' corresponding to the logit-transformed scale.
#'
#' @return
#' \itemize{
#'   \item \code{logit(p)} returns \eqn{\log(p / (1 - p))}.
#'   \item \code{logit.inv(x)} returns \eqn{1 / (1 + \exp(-x))}.
#'   \item \code{logit.inv.dev1(x)} returns the first derivative of
#'   the inverse logit, i.e., \eqn{\exp(x) / (1 + \exp(x))^2}.
#' }
#'
#' @seealso
#' Used in parameter transformations within
#' \code{\link{gmcmLik}}, \code{\link{get.idr.gmcm}},
#' and \code{\link{choosePara}}.
#'
#' @export
logit <- function(p) {
  if (any(p <= 0 | p >= 1)) {
    stop("All input values to logit() must be strictly between 0 and 1.")
  }
  log(p / (1 - p))
}

#' @rdname logit
#' @export
logit.inv <- function(x) {
  1 / (1 + exp(-x))
}

#' @rdname logit
#' @export
logit.inv.dev1 <- function(x) {
  exp(x) / (1 + exp(x))^2
}
