test_that("mergeHiC returns a valid merged Hi-C data frame", {
  set.seed(1)

  rep1 <- system.file("extdata", "rep1.txt", package = "SIDR")
  rep2 <- system.file("extdata", "rep2.txt", package = "SIDR")

  hic_data <- mergeHiC(rep1, rep2,
                       chrI = "chr18",
                       chrJ = "chr18")

  expect_s3_class(hic_data, "data.frame")
  expect_gt(nrow(hic_data), 0)
})
