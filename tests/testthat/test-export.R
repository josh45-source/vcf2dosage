make_test_dosage <- function() {
  matrix(c(0, 1, 2, 0, 1, 0),
    nrow = 2, ncol = 3, byrow = TRUE,
    dimnames = list(c("chr1_100", "chr1_200"), c("s1", "s2", "s3"))
  )
}

test_that("export_rrblup transposes and centers", {
  dosage <- make_test_dosage()
  result <- export_rrblup(dosage)

  expect_equal(dim(result), c(3, 2)) # samples x markers
  expect_equal(rownames(result), c("s1", "s2", "s3"))
  # Centered: 0,1,2 -> -1,0,1
  expect_equal(result[1, 1], -1)
  expect_equal(result[2, 1], 0)
  expect_equal(result[3, 1], 1)
})

test_that("export_rrblup writes to file", {
  dosage <- make_test_dosage()
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))

  result <- export_rrblup(dosage, file = tmp)
  expect_true(file.exists(tmp))
})

test_that("export_gwaspoly creates correct structure", {
  dosage <- make_test_dosage()
  variant_info <- tibble::tibble(
    CHROM = c("chr1", "chr1"),
    POS = c(100L, 200L),
    ID = c("snp1", "snp2"),
    REF = c("A", "C"),
    ALT = c("T", "G")
  )

  result <- export_gwaspoly(dosage, variant_info)
  expect_true("Marker" %in% names(result))
  expect_true("Chrom" %in% names(result))
  expect_true("Position" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("export_bglr transposes and imputes", {
  dosage <- matrix(c(0, NA, 2, 0, 1, 0),
    nrow = 2, ncol = 3, byrow = TRUE,
    dimnames = list(c("snp1", "snp2"), c("s1", "s2", "s3"))
  )

  result <- export_bglr(dosage, impute = TRUE)
  expect_equal(dim(result), c(3, 2))
  expect_false(any(is.na(result)))
})

test_that("export_sommer basic transpose works", {
  dosage <- make_test_dosage()
  result <- export_sommer(dosage)

  expect_equal(dim(result), c(3, 2)) # samples x markers
})

test_that("export_sommer relationship matrix", {
  dosage <- make_test_dosage()
  result <- export_sommer(dosage, as_relationship = TRUE)

  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
})

test_that("export_marker_map creates correct tibble", {
  variant_info <- tibble::tibble(
    CHROM = c("chr1", "chr2"),
    POS = c(100L, 200L),
    ID = c("snp1", "snp2")
  )

  result <- export_marker_map(variant_info)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(c("marker", "chrom", "pos") %in% names(result)))
})
