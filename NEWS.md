# vcf2dosage (development version)

## vcf2dosage 0.1.0

### New features

* `read_vcf()` reads VCF files using `data.table::fread()` for speed.
* `extract_gt()` parses genotype fields from complex FORMAT strings.
* `gt_to_dosage()` converts genotype calls to numeric dosage for any ploidy.
* `detect_ploidy()` auto-detects ploidy from genotype calls.
* `filter_maf()`, `filter_missingness()`, `filter_samples()` for variant/sample QC.
* `variant_summary()` computes per-variant quality metrics.
* `vcf2dosage()` main wrapper: read, parse, filter, convert in one call.
* Export functions for rrBLUP, GWASpoly, BGLR, sommer, and generic marker maps.
* `create_example_vcf()` generates test VCF files (diploid, tetraploid, mixed).
