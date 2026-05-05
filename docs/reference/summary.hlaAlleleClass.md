# Summarize a “hlaAlleleClass” object

Show the information of a
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
object.

## Usage

``` r
# S3 method for class 'hlaAlleleClass'
summary(object, show=TRUE, ...)
```

## Arguments

- object:

  an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)

- show:

  if TRUE, show information

- ...:

  further arguments passed to or from other methods

## Value

Return a `data.frame` of count and frequency for each HLA allele.

## Author

Xiuwen Zheng

## See also

[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")
summary(hla)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 60
#> # of unique KIR3DL1/S1 alleles: 14
#> # of unique KIR3DL1/S1 genotypes: 29
```
