# Convert a VCF File to a Dosage Matrix

The main entry point for vcf2dosage. Reads a VCF file, extracts genotype
calls, converts to dosage, filters by MAF and missingness, and returns a
ready-to-use dosage matrix.

## Usage

``` r
vcf2dosage(
  file,
  ploidy = NULL,
  min_maf = 0.05,
  max_missing_variant = 0.2,
  max_missing_sample = 0.5,
  n_max = Inf
)
```

## Arguments

- file:

  Character. Path to a VCF file (`.vcf` or `.vcf.gz`).

- ploidy:

  Integer or NULL. Expected ploidy level. If NULL, auto-detected from
  genotype calls. Default NULL.

- min_maf:

  Numeric. Minimum minor allele frequency. Variants below this threshold
  are removed. Default 0.05.

- max_missing_variant:

  Numeric. Maximum proportion of missing genotypes per variant (0-1).
  Default 0.2.

- max_missing_sample:

  Numeric. Maximum proportion of missing genotypes per sample (0-1).
  Default 0.5.

- n_max:

  Integer or Inf. Maximum number of variants to read. Default Inf.

## Value

An S3 object of class `"vcf2dosage_result"` containing:

- dosage:

  Numeric matrix (variants x samples) of allele dosages.

- dosage_t:

  Numeric matrix (samples x markers) — transposed, ready for most GS
  packages.

- map:

  Tibble of variant positions (marker, chrom, pos).

- samples:

  Character vector of retained sample names.

- ploidy:

  Detected or specified ploidy level.

- filter_summary:

  Tibble summarizing how many variants/samples were removed at each
  step.

- vcf:

  The original parsed VCF data.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- vcf2dosage("my_genotypes.vcf.gz", ploidy = 2)
result

# Use with rrBLUP
Z <- export_rrblup(result$dosage)
# rrBLUP::mixed.solve(y = pheno, Z = Z)

# Tetraploid potato
result <- vcf2dosage("potato.vcf.gz", ploidy = 4, min_maf = 0.01)
} # }
```
