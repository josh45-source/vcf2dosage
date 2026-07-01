# ---- VCF File Reading ----


#' Read a VCF File
#'
#' Reads a VCF file and returns a structured list containing the header
#' metadata, variant information (CHROM, POS, ID, REF, ALT, QUAL, FILTER,
#' INFO, FORMAT), and genotype matrix.
#'
#' Uses `data.table::fread()` for fast reading of large files.
#'
#' @param file Character. Path to a VCF file (`.vcf` or `.vcf.gz`).
#' @param n_max Integer or Inf. Maximum number of variants to read.
#'   Default `Inf` (read all).
#'
#' @return An S3 object of class `"vcf_data"` with components:
#' \describe{
#'   \item{meta}{Character vector of header lines (lines starting with `##`).}
#'   \item{fix}{Tibble of fixed variant fields: CHROM, POS, ID, REF, ALT,
#'     QUAL, FILTER, INFO.}
#'   \item{gt}{Character matrix of genotype fields (samples x variants
#'     after transposing), with the FORMAT column parsed.}
#'   \item{samples}{Character vector of sample names.}
#'   \item{format}{Character vector of FORMAT field entries per variant.}
#'   \item{n_variants}{Integer. Number of variants.}
#'   \item{n_samples}{Integer. Number of samples.}
#' }
#'
#' @examples
#' \dontrun{
#' vcf <- read_vcf("my_genotypes.vcf.gz")
#' vcf
#' }
#'
#' @export
read_vcf <- function(file, n_max = Inf) {
  if (!file.exists(file)) {
    cli_abort("File not found: {.path {file}}")
  }

  cli_alert_info("Reading VCF: {.path {basename(file)}}")

  # Read header lines
  con <- if (grepl("\\.gz$", file)) gzcon(file(file, "rb")) else file(file, "r")
  on.exit(close(con), add = TRUE)

  meta <- character()
  header_line <- NULL
  line_count <- 0L

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    line_count <- line_count + 1L

    if (startsWith(line, "##")) {
      meta <- c(meta, line)
    } else if (startsWith(line, "#CHROM")) {
      header_line <- line
      break
    } else {
      break
    }
  }
  close(con)
  on.exit(NULL)

  if (is.null(header_line)) {
    cli_abort("No valid VCF header found (missing #CHROM line).")
  }

  # Parse column names from header
  col_names <- strsplit(sub("^#", "", header_line), "\t")[[1]]

  # Read data with data.table::fread for speed
  nrows_arg <- if (is.infinite(n_max)) -1L else as.integer(n_max)

  dt <- fread(
    file,
    skip = line_count,
    header = FALSE,
    sep = "\t",
    nrows = nrows_arg,
    col.names = col_names,
    colClasses = "character",
    showProgress = FALSE
  )

  # Split into fixed fields and genotype data
  fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  format_col <- "FORMAT"
  sample_names <- setdiff(col_names, c(fixed_cols, format_col))

  fix <- as_tibble(dt[, fixed_cols, with = FALSE])
  fix$POS <- as.integer(fix$POS)

  # Extract FORMAT and genotype matrix
  format_vec <- dt[[format_col]]
  gt_matrix <- as.matrix(dt[, sample_names, with = FALSE])

  n_variants <- nrow(gt_matrix)
  n_samples <- length(sample_names)

  cli_alert_success("Read {n_variants} variants x {n_samples} samples.")

  structure(
    list(
      meta = meta,
      fix = fix,
      gt = gt_matrix,
      samples = sample_names,
      format = format_vec,
      n_variants = n_variants,
      n_samples = n_samples,
      source_file = basename(file)
    ),
    class = "vcf_data"
  )
}


#' Print a VCF Data Object
#'
#' @param x A `vcf_data` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.vcf_data <- function(x, ...) {
  cli_h3("VCF Data")
  cli_ul(c(
    "Source:    {x$source_file}",
    "Variants: {x$n_variants}",
    "Samples:  {x$n_samples}",
    "Chroms:   {length(unique(x$fix$CHROM))}"
  ))
  invisible(x)
}
