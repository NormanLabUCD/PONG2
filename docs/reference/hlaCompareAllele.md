# Evaluate prediction accuracies

To evaluate the overall accuracy, sensitivity, specificity, positive
predictive value, negative predictive value.

## Usage

``` r
hlaCompareAllele(TrueHLA, PredHLA, allele.limit=NULL, call.threshold=NaN,
  max.resolution="", output.individual=FALSE, verbose=TRUE)
```

## Arguments

- TrueHLA:

  an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md),
  the true HLA types

- PredHLA:

  an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md),
  the predicted HLA types

- allele.limit:

  a list of HLA alleles, the validation samples are limited to those
  having HLA alleles in `allele.limit`, or `NULL` for no limit.
  `allele.limit` could be character-type,
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
  or
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- call.threshold:

  the call threshold for posterior probability, i.e., call or no call is
  determined by whether `prob >= call.threshold` or not

- max.resolution:

  "2-digit", "4-digit", "6-digit", "8-digit", "allele", "protein", "2",
  "4", "6", "8", "full" or "": "allele" = "2-digit", "protein" =
  "4-digit", "full" and "" indicating no limit on resolution

- output.individual:

  if TRUE, output accuracy for each individual

- verbose:

  if TRUE, show information

## Value

Return a `list(overall, confusion, detail)`, or
`list(overall, confusion, detail, individual)` if
`output.individual=TRUE`.

`overall` (data.frame):

- total.num.ind:

  the total number of individuals

- crt.num.ind:

  the number of individuals with correct HLA types

- crt.num.haplo:

  the number of chromosomes with correct HLA alleles

- acc.ind:

  the proportion of individuals with correctly predicted HLA types
  (i.e., both of alleles are correct, the accuracy of an individual is 0
  or 1.)

- acc.haplo:

  the proportion of chromosomes with correctly predicted HLA alleles
  (i.e., the accuracy of an individual is 0, 0.5 or 1, since an
  individual has two alleles.)

- call.threshold:

  call threshold, if it is `NaN`, no call threshold is executed

- n.call:

  the number of individuals with call

- call.rate:

  overall call rate

`confusion` (matrix): a confusion matrix.

`detail` (data.frame):

- allele:

  HLA alleles

- train.num:

  the number of training haplotypes

- train.freq:

  the training haplotype frequencies

- valid.num:

  the number of validation haplotypes

- valid.freq:

  the validation haplotype frequencies

- call.rate:

  the call rates for HLA alleles

- accuracy:

  allele accuracy

- sensitivity:

  sensitivity

- specificity:

  specificity

- ppv:

  positive predictive value

- npv:

  negative predictive value

- miscall:

  the most likely miss-called alleles

- miscall.prop:

  the proportions of the most likely miss-called allele in all
  miss-called alleles

`individual` (data.frame):

- sample.id:

  sample id

- true.hla:

  the true HLA type

- pred.hla:

  the prediction of HLA type

- accuracy:

  accuracy, 0, 0.5, or 1

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`predict.hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/predict.hlaAttrBagClass.md),
[`hlaReport`](https://normanlabucd.github.io/PONG2/reference/hlaReport.md)
