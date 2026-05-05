# Summarize a “hlaAttrBagObj” object

Show the information of a
[`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)
object

## Usage

``` r
# S3 method for class 'hlaAttrBagObj'
summary(object, show=TRUE, ...)
```

## Arguments

- object:

  an object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

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

[`plot.hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/plot.hlaAttrBagObj.md),
[`print.hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/print.hlaAttrBagObj.md)
