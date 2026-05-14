test_that("determineK returns a valid number of strata", {

  set.seed(1)

  rep1 <- system.file("extdata", "rep1.txt", package = "SIDR")
  rep2 <- system.file("extdata", "rep2.txt", package = "SIDR")

  hic_data <- mergeHiC(rep1, rep2,
                       chrI = "chr18",
                       chrJ = "chr18")

  K <- determineK(hic_data,
                  max_ns = 10,
                  corr_threshold = 0.1)

  expect_true(is.numeric(K))
  expect_true(K > 1)
})
