test_that("fit_IDR_stratified returns valid fitted results", {
  skip_on_cran()
  set.seed(1)

  rep1 <- system.file("extdata", "rep1.txt", package = "SIDR")
  rep2 <- system.file("extdata", "rep2.txt", package = "SIDR")

  hic_data <- mergeHiC(rep1, rep2,
                       chrI = "chr18",
                       chrJ = "chr18")

  K <- determineK(hic_data, max_ns = 10, corr_threshold = 0.1)

  skip_if(K < 2)

  ind <- stratifyData(hic_data, ns = K)

  init <- choosePara(hic_data, ns = K, ind = ind,
                     prp = 0.7, rho = 0.15)

  thet.ini <- c(
    logit(init$mixps),
    init$mus,
    init$sigma,
    logit(init$rho),
    init$omega
  )

  fit <- fit_IDR_stratified(hic_data, ns = K, ind = ind, thet.ini = thet.ini)

  expect_type(fit, "list")
  expect_named(fit, c("est", "idr", "IDR"))

  expect_length(fit$est, 3 * K + 2)
  expect_length(fit$idr, K)
  expect_length(fit$IDR, nrow(hic_data))

  expect_true(all(is.finite(fit$est)))
  expect_true(all(fit$IDR >= 0 & fit$IDR <= 1))
})
