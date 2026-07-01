# Filter Samples by Missingness

Removes samples where the proportion of missing genotypes exceeds a
threshold.

## Usage

``` r
filter_samples(dosage, max_missing = 0.3)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- max_missing:

  Numeric. Maximum proportion of missing genotypes per sample (0-1).
  Default 0.3 (30%).

## Value

A list with `dosage` (filtered matrix), `kept` (logical vector), and
`miss_rates` (numeric vector of missingness per sample).
