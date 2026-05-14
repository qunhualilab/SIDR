test_that("choosePara returns valid initial parameters", {
  set.seed(1)

  rep1 <- system.file("extdata", "rep1.txt", package = "SIDR")
  rep2 <- system.file("extdata", "rep2.txt", package = "SIDR")

  hic_data <- mergeHiC(rep1, rep2,
                       chrI = "chr18",
                       chrJ = "chr18")

  K <- determineK(hic_data, max_ns = 10, corr_threshold = 0.1)

  skip_if(K < 2)

  ind <- stratifyData(hic_data, ns = K)

  init <- choosePara(hic_data, ns = K, ind = ind, prp = 0.7, rho = 0.15)

  expect_type(init, "list")
  expect_named(init, c("mixps", "mus", "sigma", "rho", "omega"))

  expect_length(init$mixps, K)
  expect_length(init$mus, K)
  expect_length(init$omega, K)

  expect_true(all(init$mixps > 0 & init$mixps < 1))
  expect_true(all(init$mus > 0))

})
