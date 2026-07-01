# Export Dosage Matrix for GWASpoly

Exports the dosage matrix in the format expected by
`GWASpoly::read.GWASpoly()`: a data frame with columns for marker info
(Marker, Chrom, Position) followed by one column per sample containing
integer dosage values.

## Usage

``` r
export_gwaspoly(dosage, variant_info, file = NULL)
```

## Arguments

- dosage:

  Numeric matrix. Dosage matrix (variants x samples).

- variant_info:

  Tibble. Variant information with at least CHROM and POS columns (from
  `vcf$fix`).

- file:

  Character or NULL. If provided, writes to a CSV file.

## Value

A data frame in GWASpoly format. Returned invisibly if `file` is
provided.
