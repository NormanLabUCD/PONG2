# Check the SNP predictors in a HIBAG model

Check the SNP predictors in a HIBAG model, by calculating the
overlapping between the model and SNP genotypes.

## Usage

``` r
hlaCheckSNPs(model, object,
  match.type=c("RefSNP+Position", "RefSNP", "Position"), verbose=TRUE)
```

## Arguments

- model:

  an object of
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md),
  or an object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- object:

  a genotype object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md),
  or a character vector like c("rs2523442", "rs9257863", ...)

- match.type:

  `"RefSNP+Position"` (by default) – using both of RefSNP IDs and
  positions; `"RefSNP"` – using RefSNP IDs only; `"Position"` – using
  positions only

- verbose:

  if TRUE, show information

## Value

return a `data.frame` for individual classifiers:

- NumOfValidSNP:

  the number of non-missing SNPs in an individual classifier

- NumOfSNP:

  the number of SNP predictors in an individual classifier

- fraction:

  NumOfValidSNP / NumOfSNP

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`predict.hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/predict.hlaAttrBagClass.md)
