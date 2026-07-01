# ---- Internal Utilities ----


#' Validate a VCF Data Object
#'
#' @param vcf Object to validate.
#' @return Invisible NULL. Throws error if invalid.
#' @keywords internal
validate_vcf <- function(vcf) {
  if (!inherits(vcf, "vcf_data")) {
    cli_abort(c(
      "{.arg vcf} must be a {.cls vcf_data} object.",
      "i" = "Create one with {.fn read_vcf}."
    ))
  }
  invisible(NULL)
}


#' Create a Minimal Example VCF File
#'
#' Writes a small VCF file with known genotypes for testing.
#' Includes diploid, tetraploid, and mixed-ploidy examples.
#'
#' @param file Character. Path to write the VCF file.
#' @param type Character. Type of example: `"diploid"`, `"tetraploid"`,
#'   or `"mixed"`. Default `"diploid"`.
#' @param n_variants Integer. Number of variants. Default 10.
#' @param n_samples Integer. Number of samples. Default 5.
#'
#' @return The file path (invisibly).
#'
#' @examples
#' \dontrun{
#' tmp <- tempfile(fileext = ".vcf")
#' create_example_vcf(tmp, type = "diploid")
#' vcf <- read_vcf(tmp)
#' }
#'
#' @export
create_example_vcf <- function(file,
                               type = "diploid",
                               n_variants = 10L,
                               n_samples = 5L) {
  type <- match.arg(type, c("diploid", "tetraploid", "mixed"))

  sample_names <- paste0("Sample", seq_len(n_samples))

  # Header
  header <- c(
    "##fileformat=VCFv4.3",
    "##source=vcf2dosage_example",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", sample_names
    ), collapse = "\t")
  )

  # Generate genotype data
  set.seed(42)
  lines <- character(n_variants)

  for (i in seq_len(n_variants)) {
    chrom <- paste0("chr", ((i - 1) %% 3) + 1)
    pos <- i * 1000
    id <- paste0("snp", i)
    ref <- sample(c("A", "C", "G", "T"), 1)
    alt <- sample(setdiff(c("A", "C", "G", "T"), ref), 1)

    genos <- character(n_samples)
    for (j in seq_len(n_samples)) {
      if (type == "diploid") {
        alleles <- sample(0:1, 2, replace = TRUE, prob = c(0.7, 0.3))
        genos[j] <- paste(alleles, collapse = "/")
      } else if (type == "tetraploid") {
        alleles <- sample(0:1, 4, replace = TRUE, prob = c(0.7, 0.3))
        genos[j] <- paste(alleles, collapse = "/")
      } else {
        # Mixed: some diploid, some tetraploid
        if (j <= 2) {
          alleles <- sample(0:1, 2, replace = TRUE, prob = c(0.7, 0.3))
        } else {
          alleles <- sample(0:1, 4, replace = TRUE, prob = c(0.7, 0.3))
        }
        genos[j] <- paste(alleles, collapse = "/")
      }
      # Add depth info
      dp <- sample(10:50, 1)
      genos[j] <- paste(genos[j], dp, sep = ":")
    }

    # Add one missing genotype randomly
    if (i == 3) genos[2] <- "./.:0"
    if (type == "tetraploid" && i == 3) genos[2] <- "./././.:0"

    lines[i] <- paste(c(
      chrom, pos, id, ref, alt, ".", "PASS", ".", "GT:DP",
      genos
    ), collapse = "\t")
  }

  writeLines(c(header, lines), file)
  cli_alert_success("Created example VCF ({type}): {n_variants} variants x {n_samples} samples.")
  invisible(file)
}
