# Make a SNP haplotype object

To create a
[`hlaSNPHaploClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPHaploClass.md)
object (SNP haplotypic object).

## Usage

``` r
hlaMakeSNPHaplo(haplotype, sample.id, snp.id, snp.position,
  A.allele, B.allele, assembly=c("auto", "hg18", "hg19", "unknown"))
```

## Arguments

- haplotype:

  a haplotype matrix, “# of SNPs” - by - “2 x \# of individuals”

- sample.id:

  a vector of sample ids

- snp.id:

  a vector of SNP ids

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

`haplotype` is a numeric matrix, with an entry value 0 standing for B
(ZERO A allele), 1 for A (ONE A allele) and others for missing values
(missing genotypes are usually set to be NA).

## Value

Return a
[`hlaSNPHaploClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPHaploClass.md)
object, and it is a list:

- haplotype:

  a haplotype matrix, “# of SNPs” - by - “2 x \# of individuals”

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

[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md),
[`hlaGenoSubset`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSubset.md),
[`hlaHaploSubset`](https://normanlabucd.github.io/PONG2/reference/hlaHaploSubset.md)
