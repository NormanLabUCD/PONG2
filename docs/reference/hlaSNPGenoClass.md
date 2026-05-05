# The class of SNP genotypes

The class of SNP genotypes, and its instance is returned from
[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md).

## Value

There are five components:

- genotype:

  a genotype matrix, “# of SNPs” - by - “# of individuals”

- sample.id:

  a vector of sample IDs

- snp.id:

  a vector of SNP IDs

- snp.position:

  a vector of SNP positions in basepair

- snp.allele:

  a vector of characters with a format of “A allele/B allele”

- assembly:

  the human genome reference, such like "hg19"

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md)
