# ---- Main Wrapper ----


#' Convert a VCF File to a Dosage Matrix
#'
#' The main entry point for vcf2dosage. Reads a VCF file, extracts genotype
#' calls, converts to dosage, filters by MAF and missingness, and returns
#' a ready-to-use dosage matrix.
#'
#' @param file Character. Path to a VCF file (`.vcf` or `.vcf.gz`).
#' @param ploidy Integer or NULL. Expected ploidy level. If NULL,
#'   auto-detected from genotype calls. Default NULL.
#' @param min_maf Numeric. Minimum minor allele frequency. Variants below
#'   this threshold are removed. Default 0.05.
#' @param max_missing_variant Numeric. Maximum proportion of missing
#'   genotypes per variant (0-1). Default 0.2.
#' @param max_missing_sample Numeric. Maximum proportion of missing
#'   genotypes per sample (0-1). Default 0.5.
#' @param n_max Integer or Inf. Maximum number of variants to read.
#'   Default Inf.
#'
#' @return An S3 object of class `"vcf2dosage_result"` containing:
#' \describe{
#'   \item{dosage}{Numeric matrix (variants x samples) of allele dosages.}
#'   \item{dosage_t}{Numeric matrix (samples x markers) — transposed,
#'     ready for most GS packages.}
#'   \item{map}{Tibble of variant positions (marker, chrom, pos).}
#'   \item{samples}{Character vector of retained sample names.}
#'   \item{ploidy}{Detected or specified ploidy level.}
#'   \item{filter_summary}{Tibble summarizing how many variants/samples
#'     were removed at each step.}
#'   \item{vcf}{The original parsed VCF data.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- vcf2dosage("my_genotypes.vcf.gz", ploidy = 2)
#' result
#'
#' # Use with rrBLUP
#' Z <- export_rrblup(result$dosage)
#' # rrBLUP::mixed.solve(y = pheno, Z = Z)
#'
#' # Tetraploid potato
#' result <- vcf2dosage("potato.vcf.gz", ploidy = 4, min_maf = 0.01)
#' }
#'
#' @export
vcf2dosage <- function(file,
                       ploidy = NULL,
                       min_maf = 0.05,
                       max_missing_variant = 0.2,
                       max_missing_sample = 0.5,
                       n_max = Inf) {
  cli_h1("vcf2dosage: VCF to Dosage Matrix Conversion")

  # ---- 1. Read VCF ----
  cli_h2("Step 1: Reading VCF")
  vcf <- read_vcf(file, n_max = n_max)

  n_start_variants <- vcf$n_variants
  n_start_samples <- vcf$n_samples

  # ---- 2. Extract genotypes ----
  cli_h2("Step 2: Extracting genotype calls")
  gt <- extract_gt(vcf)

  # ---- 3. Detect ploidy ----
  if (is.null(ploidy)) {
    cli_h2("Step 3: Detecting ploidy")
    ploidy_vec <- detect_ploidy(gt)
    ploidy <- as.integer(names(sort(table(ploidy_vec), decreasing = TRUE))[1])
    cli_alert_info("Detected ploidy: {ploidy} (from {length(unique(ploidy_vec))} unique level(s))")

    if (length(unique(ploidy_vec[!is.na(ploidy_vec)])) > 1) {
      cli_alert_warning("Mixed ploidy detected! Using mode = {ploidy}.")
    }
  } else {
    cli_alert_info("Using specified ploidy: {ploidy}")
  }

  # ---- 4. Convert to dosage ----
  cli_h2("Step 4: Converting to dosage")
  dosage <- gt_to_dosage(gt)
  rownames(dosage) <- paste0(vcf$fix$CHROM, "_", vcf$fix$POS)
  cli_alert_success("Dosage matrix: {nrow(dosage)} variants x {ncol(dosage)} samples.")

  # ---- 5. Filter samples ----
  cli_h2("Step 5: Filtering")
  sample_result <- filter_samples(dosage, max_missing = max_missing_sample)
  dosage <- sample_result$dosage

  # Filter variants by missingness
  miss_result <- filter_missingness(dosage, max_missing = max_missing_variant)
  dosage <- miss_result$dosage

  # Filter by MAF
  maf_result <- filter_maf(dosage, maf = min_maf, ploidy = ploidy)
  dosage <- maf_result$dosage

  # ---- 6. Build marker map ----
  # Align variant info with filtered dosage
  kept_variants <- which(miss_result$kept)[maf_result$kept]
  kept_samples <- which(sample_result$kept)

  fix_filtered <- vcf$fix[kept_variants, ]

  map <- export_marker_map(fix_filtered)

  # Update rownames
  rownames(dosage) <- map$marker

  # ---- 7. Build filter summary ----
  filter_summary <- tibble(
    step = c("Input", "Sample filter", "Missingness filter", "MAF filter"),
    n_variants = c(
      n_start_variants,
      n_start_variants,
      sum(miss_result$kept),
      nrow(dosage)
    ),
    n_samples = c(
      n_start_samples,
      ncol(dosage),
      ncol(dosage),
      ncol(dosage)
    )
  )

  # ---- 8. Result ----
  cli_h2("Done")
  cli_alert_success(
    "Final: {nrow(dosage)} variants x {ncol(dosage)} samples (ploidy = {ploidy})"
  )

  structure(
    list(
      dosage         = dosage,
      dosage_t       = t(dosage),
      map            = map,
      samples        = colnames(dosage),
      ploidy         = ploidy,
      filter_summary = filter_summary,
      vcf            = vcf
    ),
    class = "vcf2dosage_result"
  )
}


#' Print a vcf2dosage Result
#'
#' @param x A `vcf2dosage_result` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.vcf2dosage_result <- function(x, ...) {
  cli_h3("vcf2dosage Result")
  cli_ul(c(
    "Variants: {nrow(x$dosage)}",
    "Samples:  {ncol(x$dosage)}",
    "Ploidy:   {x$ploidy}",
    "Chroms:   {length(unique(x$map$chrom))}",
    "Source:   {x$vcf$source_file}"
  ))
  cat("\nFilter summary:\n")
  print(x$filter_summary)
  invisible(x)
}
