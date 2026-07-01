# Export Dosage Matrix for sommer

Exports the dosage matrix for use with `sommer::mmer()`: a numeric
matrix with samples as rows and markers as columns. Optionally computes
the additive relationship matrix (A = ZZ'/p).

## Usage

``` r
export_sommer(dosage, as_relationship = FALSE, file = NULL)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- as_relationship:

  Logical. If TRUE, returns the genomic relationship matrix instead of
  the raw dosage. Default FALSE.

- file:

  Character or NULL. If provided, writes to a CSV file.

## Value

A numeric matrix (samples x markers), or a relationship matrix if
`as_relationship = TRUE`.
