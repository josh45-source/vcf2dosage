# Summary of Variant Quality Metrics

Calculates per-variant quality metrics: MAF, missingness rate,
heterozygosity rate.

## Usage

``` r
variant_summary(dosage, variant_info = NULL, ploidy = 2L)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- variant_info:

  Tibble or NULL. Variant information (CHROM, POS, etc.) to join with
  metrics.

- ploidy:

  Integer. Ploidy level. Default 2.

## Value

A tibble with per-variant quality metrics.
