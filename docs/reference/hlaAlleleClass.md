# Class of HLA Type

The definition of a class for HLA types, returned from
[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md).

## Value

There are following components:

- locus:

  HLA locus

- pos.start:

  the starting position in basepair

- pos.end:

  the end position in basepair

- value:

  a data frame

- assembly:

  the human genome reference, such like "hg19"

- postprob:

  a matrix of all posterior probabilities

The component `value` includes:

- sample.id:

  sample ID

- allele1:

  HLA allele

- allele2:

  HLA allele

- prob:

  the posterior probability

## Author

Xiuwen Zheng

## See also

[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md)
