#' Simulated stratified Gaussian copula data
#'
#' A small dataset of simulated paired scores from the stratified Gaussian
#' mixture copula model, used for demonstration and testing purposes.
#'
#' @format A data frame with 5,000 rows and 3 columns:
#' \describe{
#'   \item{obs1}{Numeric; simulated score from replicate 1 (-log p-value).}
#'   \item{obs2}{Numeric; simulated score from replicate 2 (-log p-value).}
#'   \item{genuine}{Logical; TRUE if the observation belongs to the reproducible component.}
#' }
#' @source Generated for demonstration.
"sim_data"


#' Example Hi-C dataset: Chromosome 18, IMR90, 50kb resolution
#'
#' A subset of Hi-C interaction data used to illustrate
#' the stratified IDR analysis workflow. P-values were generated
#' by HiC-DC+.
#'
#' @format A data frame with 29,501 rows and 9 columns:
#' \describe{
#'   \item{chromosomeI}{Character; chromosome name of fragment I.}
#'   \item{chromosomeJ}{Character; chromosome name of fragment J.}
#'   \item{fragmentI}{Integer; coordinate of fragment I start.}
#'   \item{fragmentJ}{Integer; coordinate of fragment J start.}
#'   \item{pvalue1}{Numeric; p-value for replicate 1.}
#'   \item{pvalue2}{Numeric; p-value for replicate 2.}
#'   \item{obs1}{Numeric; -log(pvalue1).}
#'   \item{obs2}{Numeric; -log(pvalue2).}
#'   \item{dist}{Numeric; genomic distance between fragments (kb).}
#' }
#' @source Generated for demonstration.
"hic_example"




