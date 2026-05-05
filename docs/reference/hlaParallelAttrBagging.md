# Build a HIBAG model via parallel computation

To build a HIBAG model for predicting HLA types via parallel
computation.

## Usage

``` r
hlaParallelAttrBagging(cl, hla, snp, auto.save="",
  nclassifier=100, mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
  stop.cluster=FALSE, verbose=TRUE)
```

## Arguments

- cl:

  a cluster object, created by the package
  [parallel](https://rdrr.io/r/parallel/parallel-package.html) or
  [snow](http://CRAN.R-project.org/package=snow); if `NULL` is given, a
  uniprocessor implementation will be performed

- hla:

  training HLA types, an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)

- snp:

  training SNP genotypes, an object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- auto.save:

  specify a autosaved file, see details

- nclassifier:

  the total number of individual classifiers

- mtry:

  a character or a numeric value, the number of variables randomly
  sampled as candidates for each selection. See details

- prune:

  if TRUE, to perform a parsimonious forward variable selection,
  otherwise, exhaustive forward variable selection. See details

- rm.na:

  if TRUE, remove the samples with missing HLA types

- stop.cluster:

  `TRUE`: stop cluster nodes after computing

- verbose:

  if TRUE, show information

## Details

`mtry` (the number of variables randomly sampled as candidates for each
selection): `"sqrt"`, using the square root of the total number of
candidate SNPs; `"all"`, using all candidate SNPs; `"one"`, using one
SNP; `an integer`, specifying the number of candidate SNPs; `0 < r < 1`,
the number of candidate SNPs is "r \* the total number of SNPs".

`prune`: there is no significant difference on accuracy between
parsimonious and exhaustive forward variable selections. If
`prune = TRUE`, the searching algorithm performs a parsimonious forward
variable selection: if a new SNP predictor reduces the current
out-of-bag accuracy, then it is removed from the candidate SNP set for
future searching. Parsimonious selection helps to improve the
computational efficiency by reducing the searching times of
non-informative SNP markers.

If `auto.save=""`, the function returns a HIBAG model (an object of
[`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md));
otherwise, there is no return.

## Value

Return an object of
[`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
if `auto.save` is specified.

## References

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS; HIBAG –
HLA Genotype Imputation with Attribute Bagging; (Abstract 294,
Platform/Oral Talk); Present at the 62nd Annual Meeting of the American
Society of Human Genetics, November 9, 2012 in San Francisco,
California.

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS; HIBAG –
HLA Genotype Imputation with Attribute Bagging. Pharmacogenomics
Journal. doi: 10.1038/tpj.2013.18.
<http://dx.doi.org/10.1038/tpj.2013.18>

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaClose`](https://normanlabucd.github.io/PONG2/reference/hlaClose.md)

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
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel = match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'snpid' not found
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))


#############################################################################

library(parallel)

# use option cl.core to choose an appropriate cluster size.
cl <- makeCluster(getOption("cl.cores", 2))
set.seed(100)

# train a HIBAG model in parallel
# please use "nclassifier=100" when you use HIBAG for real data
hlaParallelAttrBagging(cl, hlatab$training, train.geno, nclassifier=4,
  auto.save="tmp_model.RData", stop.cluster=TRUE)
#> Error in hlaParallelAttrBagging(cl, hlatab$training, train.geno, nclassifier = 4,     auto.save = "tmp_model.RData", stop.cluster = TRUE): could not find function "hlaParallelAttrBagging"

mobj <- get(load("tmp_model.RData"))
#> Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
summary(mobj)
#> Gene: KIR2DL23
#> Training dataset: 589 samples X 2770 SNPs
#>  # of KIR3DL1/S1 alleles: 23
#>  # of individual classifiers: 100
#>  total # of SNPs used: 2092
#>  average # of SNPs in an individual classifier: 103.39, sd: 14.99, min: 63, max: 128
#>  average # of haplotypes in an individual classifier: 421.30, sd: 62.46, min: 288, max: 614
#>  average out-of-bag accuracy: 95.74%, sd: 0.89%, min: 93.09%, max: 98.39%
#> Genome assembly: hg38
model <- hlaModelFromObj(mobj)

# validation
pred <- predict(model, test.geno)
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"
summary(pred)
#> Gene: KIR2DS5
#> Range: [55281035bp, 55296300bp] on hg19
#> # of samples: 225
#> # of unique KIR3DL1/S1 alleles: 3
#> # of unique KIR3DL1/S1 genotypes: 4
#> Posterior probability:
#>    [0,0.25)  [0.25,0.5)  [0.5,0.75)    [0.75,1] 
#>    0 (0.0%)   11 (4.9%)  97 (43.1%) 117 (52.0%) 

# compare
hlaCompareAllele(hlatab$validation, pred, allele.limit=model)$overall
#> Calling 'hlaCompareAllele': there are 0 individuals in common.
#>   total.num.ind crt.num.ind crt.num.haplo acc.ind acc.haplo call.threshold
#> 1             0           0             0     NaN       NaN              0
#>   n.call call.rate
#> 1      0         1


# since 'stop.cluster=TRUE' used in 'hlaParallelAttrBagging'
# need a new cluster
cl <- makeCluster(getOption("cl.cores", 2))

pred <- predict(model, test.geno, cl=cl)
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"
summary(pred)
#> Gene: KIR2DS5
#> Range: [55281035bp, 55296300bp] on hg19
#> # of samples: 225
#> # of unique KIR3DL1/S1 alleles: 3
#> # of unique KIR3DL1/S1 genotypes: 4
#> Posterior probability:
#>    [0,0.25)  [0.25,0.5)  [0.5,0.75)    [0.75,1] 
#>    0 (0.0%)   11 (4.9%)  97 (43.1%) 117 (52.0%) 

# stop parallel nodes
stopCluster(cl)
```
