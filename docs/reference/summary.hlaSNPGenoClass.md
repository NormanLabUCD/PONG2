# Summarize the genotypic dataset

Summarize the genotypic dataset.

## Usage

``` r
# S3 method for class 'hlaSNPGenoClass'
summary(object, show=TRUE, ...)
```

## Arguments

- object:

  a genotype object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- show:

  if TRUE, print information

- ...:

  further arguments passed to or from other methods

## Value

None.

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md),
[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md),
[`hlaGenoSubset`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSubset.md),
[`hlaHaploSubset`](https://normanlabucd.github.io/PONG2/reference/hlaHaploSubset.md),
[`summary.hlaSNPHaploClass`](https://normanlabucd.github.io/PONG2/reference/summary.hlaSNPHaploClass.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

summary(HapMap_CEU_Geno)
#> SNP genotypes: 
#>  60 samples X 1564 SNPs
#>  SNPs range from 25769023bp to 33421576bp on hg19
#> Missing rate per SNP:
#>  min: 0, max: 0.0666667, mean: 0.0651215, median: 0.0666667, sd: 0.00968287
#> Missing rate per sample:
#>  min: 0, max: 0.969949, mean: 0.0651215, median: 0.000639386, sd: 0.243737
#> Minor allele frequency:
#>  min: 0, max: 0.5, mean: 0.227582, median: 0.205357, sd: 0.1389
#> Allelic information:
#> A/C A/G C/T G/T 
#> 136 655 632 141 
```
