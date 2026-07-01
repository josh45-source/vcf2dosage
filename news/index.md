# Changelog

## vcf2dosage (development version)

### vcf2dosage 0.1.0

#### New features

- [`read_vcf()`](https://josh45-source.github.io/vcf2dosage/reference/read_vcf.md)
  reads VCF files using
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  for speed.
- [`extract_gt()`](https://josh45-source.github.io/vcf2dosage/reference/extract_gt.md)
  parses genotype fields from complex FORMAT strings.
- [`gt_to_dosage()`](https://josh45-source.github.io/vcf2dosage/reference/gt_to_dosage.md)
  converts genotype calls to numeric dosage for any ploidy.
- [`detect_ploidy()`](https://josh45-source.github.io/vcf2dosage/reference/detect_ploidy.md)
  auto-detects ploidy from genotype calls.
- [`filter_maf()`](https://josh45-source.github.io/vcf2dosage/reference/filter_maf.md),
  [`filter_missingness()`](https://josh45-source.github.io/vcf2dosage/reference/filter_missingness.md),
  [`filter_samples()`](https://josh45-source.github.io/vcf2dosage/reference/filter_samples.md)
  for variant/sample QC.
- [`variant_summary()`](https://josh45-source.github.io/vcf2dosage/reference/variant_summary.md)
  computes per-variant quality metrics.
- [`vcf2dosage()`](https://josh45-source.github.io/vcf2dosage/reference/vcf2dosage.md)
  main wrapper: read, parse, filter, convert in one call.
- Export functions for rrBLUP, GWASpoly, BGLR, sommer, and generic
  marker maps.
- [`create_example_vcf()`](https://josh45-source.github.io/vcf2dosage/reference/create_example_vcf.md)
  generates test VCF files (diploid, tetraploid, mixed).
