# HIBAG model prediction (in parallel)

To predict HLA type based on a HIBAG model (in parallel).

## Usage

``` r
# S3 method for class 'hlaAttrBagClass'
predict(object, snp, cl,
  type=c("response", "prob", "response+prob"), vote=c("prob", "majority"),
  allele.check=TRUE, match.type=c("RefSNP+Position", "RefSNP", "Position"),
  same.strand=FALSE, verbose=TRUE, ...)
```

## Arguments

- object:

  a model of
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)

- snp:

  a genotypic object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- cl:

  a cluster object, created by the package
  [parallel](https://rdrr.io/r/parallel/parallel-package.html) or
  [snow](http://CRAN.R-project.org/package=snow); if `NULL` is given, a
  uniprocessor implementation will be performed

- type:

  "response": return the best-guess type plus its posterior probability;
  "prob": return all posterior probabilities; "response+prob": return
  the best-guess and all posterior probabilities

- vote:

  `"prob"` (default behavior) – make a prediction based on the averaged
  posterior probabilities from all individual classifiers; `"majority"`
  – majority voting from all individual classifiers, where each
  classifier votes for an HLA type

- allele.check:

  if `TRUE`, check and then switch allele pairs if needed

- match.type:

  `"RefSNP+Position"` (by default) – using both of RefSNP IDs and
  positions; `"RefSNP"` – using RefSNP IDs only; `"Position"` – using
  positions only

- same.strand:

  `TRUE` assuming alleles are on the same strand (e.g., forward strand);
  otherwise, `FALSE` not assuming whether on the same strand or not

- verbose:

  if TRUE, show information

- ...:

  further arguments passed to or from other methods

## Value

Return a
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
object with posterior probabilities of predicted HLA types, or a matrix
of pairwise possible HLA types with all posterior probabilities. If
`type = "response+prob"`, return a
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
object with a matrix of `postprob` for the probabilities of all pairs of
alleles. If a probability matrix is returned, `colnames` is `sample.id`
and `rownames` is an unordered pair of HLA alleles.

## Details

If more than 50% of SNP predictors are missing, a warning will be given.

When `match.type="RefSNP+Position"`, the matching of SNPs requires both
RefSNP IDs and positions. A lower missing fraction maybe gained by
matching RefSNP IDs or positions only. Call
`predict(..., match.type="RefSNP")` or
`predict(..., match.type="Position")` for this purpose. It might be safe
to assume that the SNPs with the same positions on the same genome
reference (e.g., hg19) are the same variant albeit the different RefSNP
IDs. Any concern about SNP mismatching should be emailed to the
genotyping platform provider.

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md),
[`hlaCompareAllele`](https://normanlabucd.github.io/PONG2/reference/hlaCompareAllele.md),
[`hlaParallelAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaParallelAttrBagging.md)
