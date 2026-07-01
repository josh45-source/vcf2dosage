# ---- Variant Filtering ----


#' Filter Variants by Minor Allele Frequency
#'
#' Removes variants where the minor allele frequency is below a threshold.
#' Calculates MAF from the dosage matrix.
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param maf Numeric. Minimum minor allele frequency. Default 0.05.
#' @param ploidy Integer. Ploidy level for MAF calculation. Default 2.
#'
#' @return A list with `dosage` (filtered matrix), `kept` (logical vector),
#'   and `maf_values` (numeric vector of MAF per variant).
#'
#' @export
filter_maf <- function(dosage, maf = 0.05, ploidy = 2L) {
  n_variants <- nrow(dosage)

  maf_values <- apply(dosage, 1, function(row) {
    non_na <- row[!is.na(row)]
    if (length(non_na) == 0) {
      return(0)
    }
    p <- mean(non_na) / ploidy
    min(p, 1 - p)
  })

  kept <- maf_values >= maf
  n_removed <- sum(!kept)

  cli_alert_info("MAF filter (>= {maf}): removed {n_removed}/{n_variants} variants.")

  list(
    dosage     = dosage[kept, , drop = FALSE],
    kept       = kept,
    maf_values = maf_values
  )
}


#' Filter Variants by Missingness
#'
#' Removes variants where the proportion of missing genotypes exceeds
#' a threshold.
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param max_missing Numeric. Maximum proportion of missing genotypes
#'   per variant (0-1). Default 0.2 (20%).
#'
#' @return A list with `dosage` (filtered matrix), `kept` (logical vector),
#'   and `miss_rates` (numeric vector of missingness per variant).
#'
#' @export
filter_missingness <- function(dosage, max_missing = 0.2) {
  n_variants <- nrow(dosage)
  n_samples <- ncol(dosage)

  miss_rates <- apply(dosage, 1, function(row) {
    sum(is.na(row)) / length(row)
  })

  kept <- miss_rates <= max_missing
  n_removed <- sum(!kept)

  cli_alert_info(
    "Missingness filter (<= {max_missing * 100}%): removed {n_removed}/{n_variants} variants."
  )

  list(
    dosage     = dosage[kept, , drop = FALSE],
    kept       = kept,
    miss_rates = miss_rates
  )
}


#' Filter Samples by Missingness
#'
#' Removes samples where the proportion of missing genotypes exceeds
#' a threshold.
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param max_missing Numeric. Maximum proportion of missing genotypes
#'   per sample (0-1). Default 0.3 (30%).
#'
#' @return A list with `dosage` (filtered matrix), `kept` (logical vector),
#'   and `miss_rates` (numeric vector of missingness per sample).
#'
#' @export
filter_samples <- function(dosage, max_missing = 0.3) {
  n_samples <- ncol(dosage)

  miss_rates <- apply(dosage, 2, function(col) {
    sum(is.na(col)) / length(col)
  })

  kept <- miss_rates <= max_missing
  n_removed <- sum(!kept)

  if (n_removed > 0) {
    removed_names <- colnames(dosage)[!kept]
    cli_alert_warning(
      "Sample filter: removed {n_removed}/{n_samples} sample(s) with >{max_missing * 100}% missing."
    )
  } else {
    cli_alert_info("Sample filter: all {n_samples} samples passed.")
  }

  list(
    dosage     = dosage[, kept, drop = FALSE],
    kept       = kept,
    miss_rates = miss_rates
  )
}


#' Summary of Variant Quality Metrics
#'
#' Calculates per-variant quality metrics: MAF, missingness rate,
#' heterozygosity rate.
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param variant_info Tibble or NULL. Variant information (CHROM, POS, etc.)
#'   to join with metrics.
#' @param ploidy Integer. Ploidy level. Default 2.
#'
#' @return A tibble with per-variant quality metrics.
#'
#' @export
variant_summary <- function(dosage, variant_info = NULL, ploidy = 2L) {
  n_variants <- nrow(dosage)

  summary_data <- tibble(
    variant_index = seq_len(n_variants),
    maf = apply(dosage, 1, function(row) {
      non_na <- row[!is.na(row)]
      if (length(non_na) == 0) {
        return(NA_real_)
      }
      p <- mean(non_na) / ploidy
      min(p, 1 - p)
    }),
    miss_rate = apply(dosage, 1, function(row) {
      sum(is.na(row)) / length(row)
    }),
    het_rate = apply(dosage, 1, function(row) {
      non_na <- row[!is.na(row)]
      if (length(non_na) == 0) {
        return(NA_real_)
      }
      # Heterozygous = dosage not 0 and not equal to ploidy
      sum(non_na > 0 & non_na < ploidy) / length(non_na)
    }),
    mean_dosage = apply(dosage, 1, mean, na.rm = TRUE),
    n_missing = apply(dosage, 1, function(row) sum(is.na(row)))
  )

  # Add variant info if provided
  if (!is.null(variant_info) && nrow(variant_info) == n_variants) {
    summary_data <- bind_cols(
      variant_info[, intersect(
        c("CHROM", "POS", "ID", "REF", "ALT"),
        names(variant_info)
      )],
      summary_data
    )
  }

  summary_data
}
