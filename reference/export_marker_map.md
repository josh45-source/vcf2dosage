# Export Marker Map

Exports a marker map tibble with columns for marker name, chromosome,
and position — compatible with most GS and GWAS packages.

## Usage

``` r
export_marker_map(variant_info, file = NULL)
```

## Arguments

- variant_info:

  Tibble. Variant information from `vcf$fix`.

- file:

  Character or NULL. If provided, writes to a CSV file.

## Value

A tibble with columns: `marker`, `chrom`, `pos`.
