# Make a SNP genotype object

To create a
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)
object (SNP genotypic object).

## Usage

``` r
hlaMakeSNPGeno(genotype, sample.id, snp.id, snp.position,
  A.allele, B.allele, assembly=c("auto", "hg18", "hg19", "unknown"))
```

## Arguments

- genotype:

  a genotype matrix, “# of SNPs” - by - “# of individuals”

- sample.id:

  a vector of sample IDs

- snp.id:

  a vector of SNP IDs

- snp.position:

  a vector of SNP positions

- A.allele:

  a vector of A alleles in the SNP list

- B.allele:

  a vector of B alleles in the SNP list

- assembly:

  the human genome reference: "hg19" (default), "hg18", "auto" refers to
  "hg19"

## Details

`genotype` is a numeric matrix, with an entry value 0 standing for BB
(ZERO A allele), 1 for AB (ONE A allele), 2 for AA (TWO A alleles) and
others for missing values (missing genotypes are usually set to be NA).

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

- assembly:

  the human genome reference

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md),
[`hlaGenoSubset`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSubset.md),
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

allele <- strsplit(HapMap_CEU_Geno$snp.allele, "/")
A.allele <- sapply(allele, function(x) { x[1] })
B.allele <- sapply(allele, function(x) { x[2] })

geno <- hlaMakeSNPGeno(HapMap_CEU_Geno$genotype, HapMap_CEU_Geno$sample.id,
  HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position, A.allele, B.allele,
  assembly="hg19")

summary(geno)
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
