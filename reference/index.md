# Package index

## Main Function

Convert VCF to dosage in one call.

- [`vcf2dosage()`](https://josh45-source.github.io/vcf2dosage/reference/vcf2dosage.md)
  : Convert a VCF File to a Dosage Matrix

## Reading & Parsing

Read VCF files and extract genotype data.

- [`read_vcf()`](https://josh45-source.github.io/vcf2dosage/reference/read_vcf.md)
  : Read a VCF File
- [`extract_gt()`](https://josh45-source.github.io/vcf2dosage/reference/extract_gt.md)
  : Extract Genotype Calls from VCF Data
- [`gt_to_dosage()`](https://josh45-source.github.io/vcf2dosage/reference/gt_to_dosage.md)
  : Convert Genotype Calls to Dosage Matrix
- [`detect_ploidy()`](https://josh45-source.github.io/vcf2dosage/reference/detect_ploidy.md)
  : Detect Ploidy from Genotype Calls

## Filtering

Filter variants and samples by quality metrics.

- [`filter_maf()`](https://josh45-source.github.io/vcf2dosage/reference/filter_maf.md)
  : Filter Variants by Minor Allele Frequency
- [`filter_missingness()`](https://josh45-source.github.io/vcf2dosage/reference/filter_missingness.md)
  : Filter Variants by Missingness
- [`filter_samples()`](https://josh45-source.github.io/vcf2dosage/reference/filter_samples.md)
  : Filter Samples by Missingness
- [`variant_summary()`](https://josh45-source.github.io/vcf2dosage/reference/variant_summary.md)
  : Summary of Variant Quality Metrics

## Export

Export dosage matrices for GS and GWAS packages.

- [`export_rrblup()`](https://josh45-source.github.io/vcf2dosage/reference/export_rrblup.md)
  : Export Dosage Matrix for rrBLUP
- [`export_gwaspoly()`](https://josh45-source.github.io/vcf2dosage/reference/export_gwaspoly.md)
  : Export Dosage Matrix for GWASpoly
- [`export_bglr()`](https://josh45-source.github.io/vcf2dosage/reference/export_bglr.md)
  : Export Dosage Matrix for BGLR
- [`export_sommer()`](https://josh45-source.github.io/vcf2dosage/reference/export_sommer.md)
  : Export Dosage Matrix for sommer
- [`export_marker_map()`](https://josh45-source.github.io/vcf2dosage/reference/export_marker_map.md)
  : Export Marker Map

## Utilities

- [`create_example_vcf()`](https://josh45-source.github.io/vcf2dosage/reference/create_example_vcf.md)
  : Create a Minimal Example VCF File

## Print Methods

S3 print methods for vcf2dosage’s data structures.

- [`print(`*`<vcf_data>`*`)`](https://josh45-source.github.io/vcf2dosage/reference/print.vcf_data.md)
  : Print a VCF Data Object
- [`print(`*`<vcf2dosage_result>`*`)`](https://josh45-source.github.io/vcf2dosage/reference/print.vcf2dosage_result.md)
  : Print a vcf2dosage Result
