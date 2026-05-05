# The class of SNP haplotypes

The class of SNP haplotypes, and its instance is returned from
[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md).

## Value

There are five components:

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

  the human genome reference, such like "hg19"

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md)
