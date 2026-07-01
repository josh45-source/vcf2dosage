# Export Dosage Matrix for BGLR

Exports the dosage matrix in the format expected by `BGLR::BGLR()`: a
numeric matrix with samples as rows and markers as columns. Missing
values are imputed with column means (BGLR requires complete data).

## Usage

``` r
export_bglr(dosage, impute = TRUE, file = NULL)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- impute:

  Logical. Whether to impute missing values with marker means. Default
  TRUE.

- file:

  Character or NULL. If provided, writes to a CSV file.

## Value

A numeric matrix (samples x markers).
