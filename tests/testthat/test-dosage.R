# ---- Tests for dosage conversion ----

test_that("gt_to_dosage handles diploid genotypes", {
  gt <- matrix(c("0/0", "0/1", "1/1", "0/0"), nrow = 2, ncol = 2, byrow = TRUE)
  dosage <- gt_to_dosage(gt)

  expect_equal(dosage[1, 1], 0)
  expect_equal(dosage[1, 2], 1)
  expect_equal(dosage[2, 1], 2)
  expect_equal(dosage[2, 2], 0)
})

test_that("gt_to_dosage handles phased genotypes", {
  gt <- matrix(c("0|0", "0|1", "1|1"), nrow = 3, ncol = 1)
  dosage <- gt_to_dosage(gt)

  expect_equal(dosage[1, 1], 0)
  expect_equal(dosage[2, 1], 1)
  expect_equal(dosage[3, 1], 2)
})

test_that("gt_to_dosage handles tetraploid genotypes", {
  gt <- matrix(c("0/0/0/0", "0/0/0/1", "0/0/1/1", "0/1/1/1", "1/1/1/1"),
    nrow = 5, ncol = 1
  )
  dosage <- gt_to_dosage(gt)

  expect_equal(dosage[1, 1], 0)
  expect_equal(dosage[2, 1], 1)
  expect_equal(dosage[3, 1], 2)
  expect_equal(dosage[4, 1], 3)
  expect_equal(dosage[5, 1], 4)
})

test_that("gt_to_dosage handles missing genotypes", {
  gt <- matrix(c("0/0", "./.", ".", "0/1"), nrow = 4, ncol = 1)
  dosage <- gt_to_dosage(gt)

  expect_equal(dosage[1, 1], 0)
  expect_true(is.na(dosage[2, 1]))
  expect_true(is.na(dosage[3, 1]))
  expect_equal(dosage[4, 1], 1)
})

test_that("gt_to_dosage handles partially missing genotypes", {
  gt <- matrix(c("0/.", "./1"), nrow = 2, ncol = 1)
  dosage <- gt_to_dosage(gt)

  expect_true(is.na(dosage[1, 1]))
  expect_true(is.na(dosage[2, 1]))
})

test_that("gt_to_dosage preserves dimensions and names", {
  gt <- matrix(c("0/0", "0/1", "1/1", "0/0"),
    nrow = 2, ncol = 2,
    dimnames = list(c("snp1", "snp2"), c("sampleA", "sampleB"))
  )
  dosage <- gt_to_dosage(gt)

  expect_equal(dim(dosage), c(2, 2))
  expect_equal(rownames(dosage), c("snp1", "snp2"))
  expect_equal(colnames(dosage), c("sampleA", "sampleB"))
})

test_that("detect_ploidy identifies diploid", {
  gt <- matrix(c("0/0", "0/1", "1/1", "0/1"),
    nrow = 2, ncol = 2,
    dimnames = list(NULL, c("s1", "s2"))
  )
  ploidy <- detect_ploidy(gt)
  expect_equal(as.integer(ploidy), c(2L, 2L))
})

test_that("detect_ploidy identifies tetraploid", {
  gt <- matrix(c("0/0/0/0", "0/0/0/1", "0/0/1/1", "0/1/1/1"),
    nrow = 2, ncol = 2,
    dimnames = list(NULL, c("s1", "s2"))
  )
  ploidy <- detect_ploidy(gt)
  expect_equal(as.integer(ploidy), c(4L, 4L))
})
