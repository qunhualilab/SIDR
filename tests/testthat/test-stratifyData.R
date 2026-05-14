test_that("stratifyData returns a valid list of strata", {
  set.seed(1)

  rep1 <- system.file("extdata", "rep1.txt", package = "SIDR")
  rep2 <- system.file("extdata", "rep2.txt", package = "SIDR")

  hic_data <- mergeHiC(rep1, rep2,
                       chrI = "chr18",
                       chrJ = "chr18")

  K <- determineK(hic_data,
                  max_ns = 10,
                  corr_threshold = 0.1)

  skip_if(K < 2)

  ind <- stratifyData(hic_data, ns = K)

  expect_type(ind, "list")
  expect_length(ind, K)

  all_indices <- unlist(ind)

  expect_true(all(all_indices >= 1))
  expect_true(all(all_indices <= nrow(hic_data)))
  expect_equal(sort(all_indices), seq_len(nrow(hic_data)))
})
