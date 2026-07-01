test_that("create_example_vcf creates a valid file", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))

  create_example_vcf(tmp, type = "diploid")
  expect_true(file.exists(tmp))

  lines <- readLines(tmp)
  expect_true(any(startsWith(lines, "##fileformat")))
  expect_true(any(startsWith(lines, "#CHROM")))
})

test_that("read_vcf reads a diploid VCF", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "diploid", n_variants = 5, n_samples = 3)

  vcf <- read_vcf(tmp)
  expect_s3_class(vcf, "vcf_data")
  expect_equal(vcf$n_variants, 5)
  expect_equal(vcf$n_samples, 3)
  expect_equal(length(vcf$samples), 3)
})

test_that("read_vcf reads a tetraploid VCF", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "tetraploid", n_variants = 5, n_samples = 3)

  vcf <- read_vcf(tmp)
  expect_s3_class(vcf, "vcf_data")
  expect_equal(vcf$n_variants, 5)
})

test_that("extract_gt extracts genotype field", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "diploid", n_variants = 5, n_samples = 3)

  vcf <- read_vcf(tmp)
  gt <- extract_gt(vcf)

  expect_true(is.matrix(gt))
  expect_equal(dim(gt), c(5, 3))
  # Should contain "/" separators
  expect_true(any(grepl("/", gt)))
})

test_that("print.vcf_data doesn't error", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "diploid")

  vcf <- read_vcf(tmp)
  expect_no_error(print(vcf))
})

test_that("vcf2dosage full pipeline works on diploid", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "diploid", n_variants = 10, n_samples = 5)

  result <- vcf2dosage(tmp, ploidy = 2, min_maf = 0, max_missing_variant = 1)
  expect_s3_class(result, "vcf2dosage_result")
  expect_true(is.matrix(result$dosage))
  expect_true(is.matrix(result$dosage_t))
  expect_s3_class(result$map, "tbl_df")
  expect_true(ncol(result$dosage_t) == nrow(result$dosage))
})

test_that("vcf2dosage full pipeline works on tetraploid", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "tetraploid", n_variants = 10, n_samples = 5)

  result <- vcf2dosage(tmp, ploidy = 4, min_maf = 0, max_missing_variant = 1)
  expect_s3_class(result, "vcf2dosage_result")
  expect_equal(result$ploidy, 4L)
})

test_that("vcf2dosage auto-detects ploidy", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "tetraploid", n_variants = 10, n_samples = 5)

  result <- vcf2dosage(tmp, ploidy = NULL, min_maf = 0, max_missing_variant = 1)
  expect_equal(result$ploidy, 4L)
})

test_that("print.vcf2dosage_result doesn't error", {
  tmp <- tempfile(fileext = ".vcf")
  on.exit(unlink(tmp))
  create_example_vcf(tmp, type = "diploid", n_variants = 5, n_samples = 3)

  result <- vcf2dosage(tmp, min_maf = 0, max_missing_variant = 1)
  expect_no_error(print(result))
})
