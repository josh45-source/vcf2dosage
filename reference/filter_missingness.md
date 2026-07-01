# Filter Variants by Missingness

Removes variants where the proportion of missing genotypes exceeds a
threshold.

## Usage

``` r
filter_missingness(dosage, max_missing = 0.2)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- max_missing:

  Numeric. Maximum proportion of missing genotypes per variant (0-1).
  Default 0.2 (20%).

## Value

A list with `dosage` (filtered matrix), `kept` (logical vector), and
`miss_rates` (numeric vector of missingness per variant).
