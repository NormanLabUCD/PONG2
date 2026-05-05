# Summarize a “hlaAttrBagClass” object

Show the information of a
[`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
object.

## Usage

``` r
# S3 method for class 'hlaAttrBagClass'
summary(object, show=TRUE, ...)
```

## Arguments

- object:

  an object of
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)

- show:

  if TRUE, show information

- ...:

  further arguments passed to or from other methods

## Value

Return a `list`:

- num.classifier:

  the total number of classifiers

- num.snp:

  the total number of SNPs

- snp.id:

  SNP IDs

- snp.position:

  SNP position in basepair

- snp.hist:

  the number of classifier for each SNP, and it could be used for SNP
  importance

- info:

  a `data.frame` for the average number of SNPs (`num.snp`), haplotypes
  (`num.haplo`), out-of-bag accuracies (`accuracy`) among all
  classifiers: mean, standard deviation, min, max

## Author

Xiuwen Zheng

## See also

[`plot.hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/plot.hlaAttrBagClass.md),
[`print.hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/print.hlaAttrBagClass.md)
