# Filter Variants by Minor Allele Frequency

Removes variants where the minor allele frequency is below a threshold.
Calculates MAF from the dosage matrix.

## Usage

``` r
filter_maf(dosage, maf = 0.05, ploidy = 2L)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- maf:

  Numeric. Minimum minor allele frequency. Default 0.05.

- ploidy:

  Integer. Ploidy level for MAF calculation. Default 2.

## Value

A list with `dosage` (filtered matrix), `kept` (logical vector), and
`maf_values` (numeric vector of MAF per variant).
