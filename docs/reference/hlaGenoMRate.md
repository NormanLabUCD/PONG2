# Missing Rates Per SNP

To calculate the missing rates from genotypes or haplotypes per SNP.

## Usage

``` r
hlaGenoMRate(obj)
```

## Arguments

- obj:

  an object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)
  or
  [`hlaSNPHaploClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPHaploClass.md)

## Value

Return missing rates per SNP.

## Author

Xiuwen Zheng

## See also

[`hlaGenoAFreq`](https://normanlabucd.github.io/PONG2/reference/hlaGenoAFreq.md),
[`hlaGenoMFreq`](https://normanlabucd.github.io/PONG2/reference/hlaGenoMFreq.md),
`hlaGenoMRate`,
[`hlaGenoMRate_Samp`](https://normanlabucd.github.io/PONG2/reference/hlaGenoMRate_Samp.md)

## Examples

``` r
# load SNP genotypes
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

summary(hlaGenoMRate(HapMap_CEU_Geno))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.06667 0.06667 0.06512 0.06667 0.06667 
```
