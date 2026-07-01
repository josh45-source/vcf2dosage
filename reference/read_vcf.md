# Read a VCF File

Reads a VCF file and returns a structured list containing the header
metadata, variant information (CHROM, POS, ID, REF, ALT, QUAL, FILTER,
INFO, FORMAT), and genotype matrix.

## Usage

``` r
read_vcf(file, n_max = Inf)
```

## Arguments

- file:

  Character. Path to a VCF file (`.vcf` or `.vcf.gz`).

- n_max:

  Integer or Inf. Maximum number of variants to read. Default `Inf`
  (read all).

## Value

An S3 object of class `"vcf_data"` with components:

- meta:

  Character vector of header lines (lines starting with `##`).

- fix:

  Tibble of fixed variant fields: CHROM, POS, ID, REF, ALT, QUAL,
  FILTER, INFO.

- gt:

  Character matrix of genotype fields (samples x variants after
  transposing), with the FORMAT column parsed.

- samples:

  Character vector of sample names.

- format:

  Character vector of FORMAT field entries per variant.

- n_variants:

  Integer. Number of variants.

- n_samples:

  Integer. Number of samples.

## Details

Uses
[`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
for fast reading of large files.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf <- read_vcf("my_genotypes.vcf.gz")
vcf
} # }
```
