# Extract Genotype Calls from VCF Data

Extracts the GT (genotype) field from the full genotype strings in a VCF
object. Handles both simple GT-only formats and complex multi-field
formats (e.g. GT:DP:GQ).

## Usage

``` r
extract_gt(vcf)
```

## Arguments

- vcf:

  A `vcf_data` object from
  [`read_vcf()`](https://josh45-source.github.io/vcf2dosage/reference/read_vcf.md).

## Value

A character matrix (variants x samples) containing only the GT field
(e.g. "0/0", "0/1/1/1", "0\|0\|1\|1").
