# Detect Ploidy from Genotype Calls

Examines a genotype matrix to determine the ploidy level of each sample
based on the number of alleles per genotype call.

## Usage

``` r
detect_ploidy(gt_matrix)
```

## Arguments

- gt_matrix:

  A character matrix of genotype calls (variants x samples).

## Value

A named integer vector of ploidy per sample. Returns the most common
ploidy observed across non-missing genotypes.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf <- read_vcf("tetraploid.vcf.gz")
gt <- extract_gt(vcf)
detect_ploidy(gt)
# sample1 sample2 sample3
#       4       4       4
} # }
```
