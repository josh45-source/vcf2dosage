# Convert Genotype Calls to Dosage Matrix

Converts a character matrix of VCF genotype calls (e.g. "0/0", "0/1",
"1/1" for diploids, or "0/0/1/1" for tetraploids) into a numeric dosage
matrix counting the number of alternate alleles.

## Usage

``` r
gt_to_dosage(gt_matrix, missing_value = NA_real_)
```

## Arguments

- gt_matrix:

  A character matrix of genotype calls (variants x samples), as returned
  by
  [`extract_gt()`](https://josh45-source.github.io/vcf2dosage/reference/extract_gt.md).

- missing_value:

  Numeric. Value to use for missing genotypes. Default `NA_real_`.

## Value

A numeric matrix (variants x samples) of alternate allele dosages. For a
diploid: 0, 1, or 2. For a tetraploid: 0, 1, 2, 3, or 4.

## Details

Handles both phased (`|`) and unphased (`/`) separators, missing
genotypes (`./.`, `.`), and arbitrary ploidy levels.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf <- read_vcf("genotypes.vcf.gz")
gt <- extract_gt(vcf)
dosage <- gt_to_dosage(gt)
} # }
```
