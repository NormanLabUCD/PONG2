# Convert from PLINK BED format

To convert a PLINK BED file to an object of
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md).

## Usage

``` r
hlaBED2Geno(bed.fn, fam.fn, bim.fn, rm.invalid.allele=FALSE,
  import.chr="xMHC", assembly=c("auto", "hg18", "hg19", "unknown"),
  verbose=TRUE)
```

## Arguments

- bed.fn:

  binary file, genotype information

- fam.fn:

  family, individual information, etc

- bim.fn:

  extended MAP file: two extra cols = allele names

- rm.invalid.allele:

  if TRUE, remove SNPs with invalid alleles

- import.chr:

  the chromosome, "1" .. "22", "X", "Y", "XY", "MT", "xMHC", or "",
  where "xMHC" implies the extended MHC on chromosome 6, and "" for all
  SNPs

- assembly:

  the human genome reference: "hg19" (default), "hg18", "auto" refers to
  "hg19"

- verbose:

  if TRUE, show information

## Value

Return an object of
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md).

## Author

Xiuwen Zheng

## See also

[`hlaGeno2PED`](https://normanlabucd.github.io/PONG2/reference/hlaGeno2PED.md),
[`hlaGDS2Geno`](https://normanlabucd.github.io/PONG2/reference/hlaGDS2Geno.md)

## Examples

``` r
# Import a PLINK BED file
bed.fn <- system.file("extdata", "HapMap_CEU.bed", package="HIBAG")
fam.fn <- system.file("extdata", "HapMap_CEU.fam", package="HIBAG")
bim.fn <- system.file("extdata", "HapMap_CEU.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")
#> Open "/usr/local/lib/R/site-library/HIBAG/extdata/HapMap_CEU.bed" in the individual-major mode.
#> Open "/usr/local/lib/R/site-library/HIBAG/extdata/HapMap_CEU.fam".
#> Open "/usr/local/lib/R/site-library/HIBAG/extdata/HapMap_CEU.bim".
#> Import 0 SNPs within the KIR gene cluster on chromosome 19.
#> Error in hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly = "hg19"): There is no SNP imported.

summary(hapmap.ceu)
#> Error: object 'hapmap.ceu' not found
```
