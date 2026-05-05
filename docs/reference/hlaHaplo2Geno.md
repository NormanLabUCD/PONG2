# Get a genotype object from a specified haplotype object

To get a genotype object from a specified haplotype object by combining
SNP alleles.

## Usage

``` r
hlaHaplo2Geno(hapobj)
```

## Arguments

- hapobj:

  a haplotype object of
  [`hlaSNPHaploClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPHaploClass.md)

## Details

`hapobj$haplotype` is a numeric matrix, with an entry value 0 standing
for B (ZERO A allele), 1 for A (ONE A allele) and others for missing
values (missing SNP alleles are usually set to be NA).

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

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md),
[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md),
[`hlaGenoSubset`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSubset.md),
[`hlaHaploSubset`](https://normanlabucd.github.io/PONG2/reference/hlaHaploSubset.md)
