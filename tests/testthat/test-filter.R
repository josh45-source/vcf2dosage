test_that("filter_maf removes low-MAF variants", {
  dosage <- matrix(c(
    0, 0, 0, 0, 0, # MAF = 0 (monomorphic)
    0, 1, 0, 1, 0, # MAF = 0.2
    0, 1, 1, 1, 2 # MAF = 0.5
  ), nrow = 3, ncol = 5, byrow = TRUE)

  result <- filter_maf(dosage, maf = 0.1)
  expect_equal(nrow(result$dosage), 2) # removed monomorphic
  expect_equal(length(result$maf_values), 3)
})

test_that("filter_maf keeps all when threshold is 0", {
  dosage <- matrix(c(0, 0, 0, 1, 1, 2), nrow = 2, ncol = 3, byrow = TRUE)
  result <- filter_maf(dosage, maf = 0)
  expect_equal(nrow(result$dosage), 2)
})

test_that("filter_missingness removes high-missing variants", {
  dosage <- matrix(c(
    NA, NA, NA, NA, 1, # 80% missing
    0,  1,  NA, 1,  0, # 20% missing
    0,  1,  1,  1,  0 # 0% missing
  ), nrow = 3, ncol = 5, byrow = TRUE)

  result <- filter_missingness(dosage, max_missing = 0.5)
  expect_equal(nrow(result$dosage), 2) # removed first row
})

test_that("filter_samples removes high-missing samples", {
  dosage <- matrix(
    c(
      0, NA, 1,
      1, NA, 0,
      0, NA, 1,
      1, NA, 0
    ),
    nrow = 4, ncol = 3, byrow = TRUE,
    dimnames = list(NULL, c("s1", "s2", "s3"))
  )

  result <- filter_samples(dosage, max_missing = 0.5)
  expect_equal(ncol(result$dosage), 2) # removed s2
  expect_false("s2" %in% colnames(result$dosage))
})

test_that("variant_summary returns correct structure", {
  dosage <- matrix(c(0, 1, 2, NA, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
  result <- variant_summary(dosage)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("maf", "miss_rate", "het_rate") %in% names(result)))
  expect_equal(nrow(result), 2)
})
