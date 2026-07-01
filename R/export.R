# ---- Export Functions ----
# Export dosage matrices in formats compatible with GS packages.


#' Export Dosage Matrix for rrBLUP
#'
#' Exports the dosage matrix in the format expected by `rrBLUP::mixed.solve()`
#' and `rrBLUP::A.mat()`: a numeric matrix with samples as rows and markers
#' as columns, coded as -1, 0, 1 (centered on heterozygote).
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples)
#'   from [vcf2dosage()].
#' @param file Character or NULL. If provided, writes to a CSV file.
#'
#' @return A numeric matrix (samples x markers) coded as -1, 0, 1.
#'   Returned invisibly if `file` is provided.
#'
#' @export
export_rrblup <- function(dosage, file = NULL) {
  # Transpose: rrBLUP wants samples as rows
  mat <- t(dosage)

  # Center: 0,1,2 -> -1,0,1
  mat <- mat - 1

  if (!is.null(file)) {
    write.csv(mat, file = file, quote = FALSE)
    cli_alert_success("Exported rrBLUP format to {.path {file}}")
    return(invisible(mat))
  }

  mat
}


#' Export Dosage Matrix for GWASpoly
#'
#' Exports the dosage matrix in the format expected by `GWASpoly::read.GWASpoly()`:
#' a data frame with columns for marker info (Marker, Chrom, Position) followed
#' by one column per sample containing integer dosage values.
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param variant_info Tibble. Variant information with at least CHROM and POS
#'   columns (from `vcf$fix`).
#' @param file Character or NULL. If provided, writes to a CSV file.
#'
#' @return A data frame in GWASpoly format. Returned invisibly if `file`
#'   is provided.
#'
#' @export
export_gwaspoly <- function(dosage, variant_info, file = NULL) {
  if (nrow(variant_info) != nrow(dosage)) {
    cli_abort("variant_info must have the same number of rows as dosage matrix.")
  }

  # Build marker IDs
  marker_ids <- if ("ID" %in% names(variant_info) &&
    !all(variant_info$ID == "." | is.na(variant_info$ID))) {
    variant_info$ID
  } else {
    paste0(variant_info$CHROM, "_", variant_info$POS)
  }

  result <- data.frame(
    Marker = marker_ids,
    Chrom = variant_info$CHROM,
    Position = variant_info$POS,
    dosage,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!is.null(file)) {
    write.csv(result, file = file, row.names = FALSE, quote = FALSE)
    cli_alert_success("Exported GWASpoly format to {.path {file}}")
    return(invisible(result))
  }

  result
}


#' Export Dosage Matrix for BGLR
#'
#' Exports the dosage matrix in the format expected by `BGLR::BGLR()`:
#' a numeric matrix with samples as rows and markers as columns.
#' Missing values are imputed with column means (BGLR requires
#' complete data).
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param impute Logical. Whether to impute missing values with marker
#'   means. Default TRUE.
#' @param file Character or NULL. If provided, writes to a CSV file.
#'
#' @return A numeric matrix (samples x markers).
#'
#' @export
export_bglr <- function(dosage, impute = TRUE, file = NULL) {
  mat <- t(dosage)

  if (impute && any(is.na(mat))) {
    n_missing <- sum(is.na(mat))
    # Impute with column (marker) means
    col_means <- colMeans(mat, na.rm = TRUE)
    for (j in seq_len(ncol(mat))) {
      na_idx <- is.na(mat[, j])
      if (any(na_idx)) {
        mat[na_idx, j] <- col_means[j]
      }
    }
    cli_alert_info("Imputed {n_missing} missing value(s) with marker means for BGLR.")
  }

  if (!is.null(file)) {
    write.csv(mat, file = file, quote = FALSE)
    cli_alert_success("Exported BGLR format to {.path {file}}")
    return(invisible(mat))
  }

  mat
}


#' Export Dosage Matrix for sommer
#'
#' Exports the dosage matrix for use with `sommer::mmer()`: a numeric
#' matrix with samples as rows and markers as columns. Optionally
#' computes the additive relationship matrix (A = ZZ'/p).
#'
#' @param dosage Numeric matrix. Dosage matrix (variants x samples).
#' @param as_relationship Logical. If TRUE, returns the genomic
#'   relationship matrix instead of the raw dosage. Default FALSE.
#' @param file Character or NULL. If provided, writes to a CSV file.
#'
#' @return A numeric matrix (samples x markers), or a relationship
#'   matrix if `as_relationship = TRUE`.
#'
#' @export
export_sommer <- function(dosage, as_relationship = FALSE, file = NULL) {
  mat <- t(dosage)

  if (as_relationship) {
    # Center the matrix
    p <- colMeans(mat, na.rm = TRUE) / 2
    Z <- scale(mat, center = 2 * p, scale = FALSE)
    Z[is.na(Z)] <- 0
    denom <- 2 * sum(p * (1 - p))
    mat <- (Z %*% t(Z)) / denom
    cli_alert_info("Computed genomic relationship matrix ({nrow(mat)} x {ncol(mat)}).")
  }

  if (!is.null(file)) {
    write.csv(mat, file = file, quote = FALSE)
    cli_alert_success("Exported sommer format to {.path {file}}")
    return(invisible(mat))
  }

  mat
}


#' Export Marker Map
#'
#' Exports a marker map tibble with columns for marker name, chromosome,
#' and position — compatible with most GS and GWAS packages.
#'
#' @param variant_info Tibble. Variant information from `vcf$fix`.
#' @param file Character or NULL. If provided, writes to a CSV file.
#'
#' @return A tibble with columns: `marker`, `chrom`, `pos`.
#'
#' @export
export_marker_map <- function(variant_info, file = NULL) {
  marker_ids <- if ("ID" %in% names(variant_info) &&
    !all(variant_info$ID == "." | is.na(variant_info$ID))) {
    variant_info$ID
  } else {
    paste0(variant_info$CHROM, "_", variant_info$POS)
  }

  map <- tibble(
    marker = marker_ids,
    chrom  = variant_info$CHROM,
    pos    = as.integer(variant_info$POS)
  )

  if (!is.null(file)) {
    write.csv(map, file = file, row.names = FALSE, quote = FALSE)
    cli_alert_success("Exported marker map to {.path {file}}")
    return(invisible(map))
  }

  map
}
