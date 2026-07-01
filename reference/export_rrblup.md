# Export Dosage Matrix for rrBLUP

Exports the dosage matrix in the format expected by
`rrBLUP::mixed.solve()` and `rrBLUP::A.mat()`: a numeric matrix with
samples as rows and markers as columns, coded as -1, 0, 1 (centered on
heterozygote).

## Usage

``` r
export_rrblup(dosage, file = NULL)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples) from
  [`vcf2dosage()`](https://josh45-source.github.io/vcf2dosage/reference/vcf2dosage.md).

- file:

  Character or NULL. If provided, writes to a CSV file.

## Value

A numeric matrix (samples x markers) coded as -1, 0, 1. Returned
invisibly if `file` is provided.
