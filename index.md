# vcf2dosage

**vcf2dosage** converts VCF files to numeric dosage matrices for genomic
selection and GWAS. It handles diploid and polyploid species, filters by
MAF and missingness, and exports in formats ready for rrBLUP, GWASpoly,
BGLR, sommer, and AGHmatrix.

## The Problem

Getting from a VCF file to a usable dosage matrix in R is surprisingly
painful, especially for polyploids. You either cobble together `vcfR` +
manual parsing + custom scripts, or use species-specific tools that
don’t generalize. vcf2dosage does it in one function call.

## Quick Start

``` r

library(vcf2dosage)

# One function does everything
result <- vcf2dosage("genotypes.vcf.gz")
result
#> vcf2dosage Result
#> Variants: 12,847
#> Samples:  342
#> Ploidy:   2
#> Chroms:   10

# Use with rrBLUP
Z <- export_rrblup(result$dosage)
# rrBLUP::mixed.solve(y = pheno$yield, Z = Z)

# Tetraploid potato? Just specify ploidy
result <- vcf2dosage("potato.vcf.gz", ploidy = 4)

# Export for GWASpoly
geno_df <- export_gwaspoly(result$dosage, result$vcf$fix)
```

## Polyploid Support

Most VCF tools assume diploid genotypes (`0/0`, `0/1`, `1/1`).
vcf2dosage handles arbitrary ploidy:

``` r

# Tetraploid: 0/0/0/0 -> dosage 0, 0/0/1/1 -> dosage 2, 1/1/1/1 -> dosage 4
result <- vcf2dosage("tetraploid.vcf.gz", ploidy = 4)

# Hexaploid wheat
result <- vcf2dosage("wheat.vcf.gz", ploidy = 6)

# Auto-detect ploidy from genotype calls
result <- vcf2dosage("unknown_ploidy.vcf.gz")
result$ploidy
```

## Filtering

``` r

result <- vcf2dosage(
  "genotypes.vcf.gz",
  min_maf = 0.05,              # Remove rare variants
  max_missing_variant = 0.2,   # Max 20% missing per variant
  max_missing_sample = 0.5     # Max 50% missing per sample
)
```

## Export Formats

``` r

# rrBLUP: centered (-1, 0, 1), samples x markers
Z <- export_rrblup(result$dosage)

# GWASpoly: Marker/Chrom/Position + dosage columns
geno <- export_gwaspoly(result$dosage, result$vcf$fix)

# BGLR: samples x markers, missing imputed with marker means
X <- export_bglr(result$dosage)

# sommer: samples x markers, or genomic relationship matrix
G <- export_sommer(result$dosage, as_relationship = TRUE)

# Marker map for any package
map <- export_marker_map(result$vcf$fix)
```

## Part of the Plant Breeding Toolkit

vcf2dosage works alongside:

- [brapiR2](https://github.com/josh45-source/brapiR2) - Pull genotype
  data from BrAPI servers
- [phenoQC](https://github.com/josh45-source/phenoQC) - Quality control
  for phenotypic trial data

Together: **retrieve** (brapiR2) -\> **clean** (phenoQC) -\> **prepare
genotypes** (vcf2dosage) -\> **analyze**.

## Installation

``` r

remotes::install_github("josh45-source/vcf2dosage")
```

## License

MIT

## Support This Project

If vcf2dosage has been useful to you, please consider sponsoring its
development on Patreon — it helps keep the project maintained.

[![Support on
Patreon](https://img.shields.io/badge/Patreon-Support-f96854?logo=patreon&logoColor=white)](https://www.patreon.com/c/Joshfarm/membership)
