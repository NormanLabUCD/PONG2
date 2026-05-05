# Merge prediction results from multiple HIBAG models

Return an object of
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md),
which contains predicted HLA types.

## Usage

``` r
hlaPredMerge(..., weight=NULL, equivalence=NULL)
```

## Arguments

- ...:

  The object(s) of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md),
  having a field of 'postprob', and returned by
  `predict(..., type="response+prob", vote="majority")`

- weight:

  the weight used for each prediction; if `NULL`, equal weights

- equivalence:

  a `data.frame` with two columns, the first column for new equivalent
  alleles, and the second for the alleles possibly existed in the
  object(s) passed to this function

## Details

Calculate a new probability matrix for each pair of HLA alleles, by
averaging (posterior) probabilities from all models with specified
weights. If `equivalence` is specified, multiple alleles might be
collapsed into one class.

## Value

Return a
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
object.

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md),
[`predict.hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/predict.hlaAttrBagClass.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
#> [1] "training"   "validation"
# "training"   "validation"
summary(hlatab$training)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 34
#> # of unique KIR3DL1/S1 alleles: 14
#> # of unique KIR3DL1/S1 genotypes: 23
summary(hlatab$validation)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 26
#> # of unique KIR3DL1/S1 alleles: 12
#> # of unique KIR3DL1/S1 genotypes: 14

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"
length(snpid)  # 275
#> Error: object 'snpid' not found

# training and validation genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel=match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel=match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'snpid' not found
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))

# train HIBAG models
set.seed(100)

# please use "nclassifier=100" when you use HIBAG for real data
m1 <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=2,
  verbose.detail=TRUE)
#> Error: object 'train.geno' not found
m2 <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=2,
  verbose.detail=TRUE)
#> Error: object 'train.geno' not found


# validation
pd1 <- predict(m1, test.geno, type="response+prob", vote="majority")
#> Error: object 'm1' not found
pd2 <- predict(m2, test.geno, type="response+prob", vote="majority")
#> Error: object 'm2' not found

hlaCompareAllele(hlatab$validation, pd1)$overall
#> Error: object 'pd1' not found
# acc.haplo = 0.8269231

hlaCompareAllele(hlatab$validation, pd2)$overall
#> Error: object 'pd2' not found
# acc.haplo = 0.8076923


# merge predictions from multiple models, by voting from all classifiers
pd <- hlaPredMerge(pd1, pd2, weight=c(1,1))
#> Error: object 'pd1' not found

hlaCompareAllele(hlatab$validation, pd)$overall
#> Error: object 'pd' not found
# acc.haplo = 0.9230769
```
