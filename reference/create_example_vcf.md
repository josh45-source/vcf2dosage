# Create a Minimal Example VCF File

Writes a small VCF file with known genotypes for testing. Includes
diploid, tetraploid, and mixed-ploidy examples.

## Usage

``` r
create_example_vcf(file, type = "diploid", n_variants = 10L, n_samples = 5L)
```

## Arguments

- file:

  Character. Path to write the VCF file.

- type:

  Character. Type of example: `"diploid"`, `"tetraploid"`, or `"mixed"`.
  Default `"diploid"`.

- n_variants:

  Integer. Number of variants. Default 10.

- n_samples:

  Integer. Number of samples. Default 5.

## Value

The file path (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
tmp <- tempfile(fileext = ".vcf")
create_example_vcf(tmp, type = "diploid")
vcf <- read_vcf(tmp)
} # }
```
