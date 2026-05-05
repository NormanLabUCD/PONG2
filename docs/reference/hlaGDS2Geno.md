# Convert from SNP GDS format

To convert a SNP GDS file to an object of
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md).

## Usage

``` r
hlaGDS2Geno(gds.fn, rm.invalid.allele=FALSE,
  import.chr="xMHC", assembly=c("auto", "hg18", "hg19", "unknown"),
  verbose=TRUE)
```

## Arguments

- gds.fn:

  the SNP GDS file used by the `SNPRelate` package

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
[`hlaBED2Geno`](https://normanlabucd.github.io/PONG2/reference/hlaBED2Geno.md)
