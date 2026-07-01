# ---- Dosage Matrix Conversion ----


#' Extract Genotype Calls from VCF Data
#'
#' Extracts the GT (genotype) field from the full genotype strings in
#' a VCF object. Handles both simple GT-only formats and complex
#' multi-field formats (e.g. GT:DP:GQ).
#'
#' @param vcf A `vcf_data` object from [read_vcf()].
#'
#' @return A character matrix (variants x samples) containing only
#'   the GT field (e.g. "0/0", "0/1/1/1", "0|0|1|1").
#'
#' @export
extract_gt <- function(vcf) {
  validate_vcf(vcf)

  # Find GT position in FORMAT
  format_fields <- strsplit(vcf$format, ":")

  gt_matrix <- vcf$gt

  # Extract GT field (first field in most VCFs, but check)
  result <- matrix(NA_character_, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
  dimnames(result) <- dimnames(gt_matrix)

  for (i in seq_len(nrow(gt_matrix))) {
    fields <- format_fields[[i]]
    gt_pos <- match("GT", fields)

    if (is.na(gt_pos)) {
      cli_alert_warning("No GT field at variant {i}. Skipping.")
      next
    }

    if (gt_pos == 1L && !any(grepl(":", gt_matrix[i, 1]))) {
      # Simple case: GT only, no colon-separated fields
      result[i, ] <- gt_matrix[i, ]
    } else {
      # Extract GT from colon-separated fields
      result[i, ] <- vapply(
        strsplit(gt_matrix[i, ], ":"),
        function(x) if (length(x) >= gt_pos) x[gt_pos] else NA_character_,
        character(1)
      )
    }
  }

  result
}


#' Convert Genotype Calls to Dosage Matrix
#'
#' Converts a character matrix of VCF genotype calls (e.g. "0/0", "0/1",
#' "1/1" for diploids, or "0/0/1/1" for tetraploids) into a numeric
#' dosage matrix counting the number of alternate alleles.
#'
#' Handles both phased (`|`) and unphased (`/`) separators, missing
#' genotypes (`./.`, `.`), and arbitrary ploidy levels.
#'
#' @param gt_matrix A character matrix of genotype calls (variants x samples),
#'   as returned by [extract_gt()].
#' @param missing_value Numeric. Value to use for missing genotypes.
#'   Default `NA_real_`.
#'
#' @return A numeric matrix (variants x samples) of alternate allele dosages.
#'   For a diploid: 0, 1, or 2. For a tetraploid: 0, 1, 2, 3, or 4.
#'
#' @examples
#' \dontrun{
#' vcf <- read_vcf("genotypes.vcf.gz")
#' gt <- extract_gt(vcf)
#' dosage <- gt_to_dosage(gt)
#' }
#'
#' @export
gt_to_dosage <- function(gt_matrix, missing_value = NA_real_) {
  n_variants <- nrow(gt_matrix)
  n_samples <- ncol(gt_matrix)

  dosage <- matrix(missing_value, nrow = n_variants, ncol = n_samples)
  dimnames(dosage) <- dimnames(gt_matrix)

  for (i in seq_len(n_variants)) {
    for (j in seq_len(n_samples)) {
      geno <- gt_matrix[i, j]
      if (is.na(geno) || geno == "." || geno == "./." || geno == ".|.") {
        dosage[i, j] <- missing_value
        next
      }

      # Split by / or |
      alleles <- strsplit(geno, "[/|]")[[1]]

      # Check for missing alleles
      if (any(alleles == ".")) {
        dosage[i, j] <- missing_value
        next
      }

      # Count alternate alleles (anything != "0")
      allele_nums <- suppressWarnings(as.integer(alleles))
      if (any(is.na(allele_nums))) {
        dosage[i, j] <- missing_value
        next
      }

      dosage[i, j] <- sum(allele_nums > 0L)
    }
  }

  dosage
}


#' Detect Ploidy from Genotype Calls
#'
#' Examines a genotype matrix to determine the ploidy level of each sample
#' based on the number of alleles per genotype call.
#'
#' @param gt_matrix A character matrix of genotype calls (variants x samples).
#'
#' @return A named integer vector of ploidy per sample. Returns the most
#'   common ploidy observed across non-missing genotypes.
#'
#' @examples
#' \dontrun{
#' vcf <- read_vcf("tetraploid.vcf.gz")
#' gt <- extract_gt(vcf)
#' detect_ploidy(gt)
#' # sample1 sample2 sample3
#' #       4       4       4
#' }
#'
#' @export
detect_ploidy <- function(gt_matrix) {
  n_samples <- ncol(gt_matrix)
  ploidy <- setNames(integer(n_samples), colnames(gt_matrix))

  for (j in seq_len(n_samples)) {
    # Sample some non-missing genotypes
    genos <- gt_matrix[, j]
    genos <- genos[!is.na(genos) & genos != "." & genos != "./." & genos != ".|."]

    if (length(genos) == 0) {
      ploidy[j] <- NA_integer_
      next
    }

    # Count alleles in each genotype
    n_alleles <- vapply(head(genos, 100), function(g) {
      length(strsplit(g, "[/|]")[[1]])
    }, integer(1))

    # Most common count = ploidy
    ploidy[j] <- as.integer(names(sort(table(n_alleles), decreasing = TRUE))[1])
  }

  ploidy
}
