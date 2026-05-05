# Get a subset of genotypes

To get a subset of genotypes from a
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)
object.

## Usage

``` r
hlaGenoSubset(genoobj, samp.sel=NULL, snp.sel=NULL)
```

## Arguments

- genoobj:

  a genotype object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- samp.sel:

  a logical vector, or an integer vector of indices

- snp.sel:

  a logical vector, or an integer vector of indices

## Details

`genoobj$genotype` is a numeric matrix, with an entry value 0 standing
for BB (ZERO A allele), 1 for AB (ONE A allele), 2 for AA (TWO A
alleles) and others for missing values (missing genotypes are usually
set to be NA).

## Value

Return a
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)
object, and it is a list:

- genotype:

  a genotype matrix, “# of SNPs” - by - “# of individuals”

- sample.id:

  a vector of sample IDs

- snp.id:

  a vector of SNP IDs

- snp.position:

  a vector of SNP positions in basepair

- snp.allele:

  a vector of characters with the format of “A allele/B allele”

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md),
[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md),
[`hlaHaploSubset`](https://normanlabucd.github.io/PONG2/reference/hlaHaploSubset.md),
[`hlaGenoCombine`](https://normanlabucd.github.io/PONG2/reference/hlaGenoCombine.md)

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

geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = (hlaGenoMFreq(HapMap_CEU_Geno)>0.10))
summary(geno)
#> SNP genotypes: 
#>  60 samples X 1197 SNPs
#>  SNPs range from 25769023bp to 33421576bp on hg19
#> Missing rate per SNP:
#>  min: 0, max: 0.0666667, mean: 0.0657059, median: 0.0666667, sd: 0.00757446
#> Missing rate per sample:
#>  min: 0, max: 0.978279, mean: 0.0657059, median: 0.000835422, sd: 0.245786
#> Minor allele frequency:
#>  min: 0.101695, max: 0.5, mean: 0.278734, median: 0.267857, sd: 0.117338
#> Allelic information:
#> A/C A/G C/T G/T 
#> 105 511 476 105 
```
