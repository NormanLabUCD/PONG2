# Build a HIBAG model

To build a HIBAG model for predicting HLA types.

## Usage

``` r
hlaAttrBagging(hla, snp, nclassifier=100, mtry=c("sqrt", "all", "one"),
  prune=TRUE, rm.na=TRUE, verbose=TRUE, verbose.detail=FALSE)
```

## Arguments

- hla:

  the training HLA types, an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)

- snp:

  the training SNP genotypes, an object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

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

- verbose:

  if TRUE, show information

- verbose.detail:

  if TRUE, show more information

## Details

`mtry` (the number of variables randomly sampled as candidates for each
selection): `"sqrt"`, using the square root of the total number of
candidate SNPs; `"all"`, using all candidate SNPs; `"one"`, using one
SNP; `an integer`, specifying the number of candidate SNPs; `0 < r < 1`,
the number of candidate SNPs is "r \* the total number of SNPs".

`prune`: there is no significant difference on accuracy between
parsimonious and exhaustive forward variable selections. If
`prune=TRUE`, the searching algorithm performs a parsimonious forward
variable selection: if a new SNP predictor reduces the current
out-of-bag accuracy, then it is removed from the candidate SNP set for
future searching. Parsimonious selection helps to improve the
computational efficiency by reducing the searching times on
non-informative SNP markers.

A parallel version of `hlaAttrBagging` is
[`hlaParallelAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaParallelAttrBagging.md).

## Value

Return an object of
[`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md):

- n.samp:

  the total number of training samples

- n.snp:

  the total number of candidate SNP predictors

- sample.id:

  the sample IDs

- snp.id:

  the SNP IDs

- snp.position:

  SNP position in basepair

- snp.allele:

  a vector of characters with the format of “A allele/B allele”

- snp.allele.freq:

  the allele frequencies

- hla.locus:

  the name of HLA locus

- hla.allele:

  the HLA alleles used in the model

- hla.freq:

  the HLA allele frequencies

- assembly:

  the human genome reference, such like "hg19"

- model:

  internal use

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

[`hlaClose`](https://normanlabucd.github.io/PONG2/reference/hlaClose.md),
[`hlaParallelAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaParallelAttrBagging.md),
[`summary.hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/summary.hlaAttrBagClass.md),
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

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4,
  verbose.detail=TRUE)
#> Error: object 'train.geno' not found
summary(model)
#> Gene: KIR2DL23
#> Training dataset: 589 samples X 2770 SNPs
#>  # of KIR3DL1/S1 alleles: 23
#>  # of individual classifiers: 100
#>  total # of SNPs used: 2092
#>  average # of SNPs in an individual classifier: 103.39, sd: 14.99, min: 63, max: 128
#>  average # of haplotypes in an individual classifier: 421.30, sd: 62.46, min: 288, max: 614
#>  average out-of-bag accuracy: 95.74%, sd: 0.89%, min: 93.09%, max: 98.39%
#> Genome assembly: hg38

# validation
pred <- predict(model, test.geno)
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"
# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
  call.threshold=0))
#> Calling 'hlaCompareAllele': there are 0 individuals in common.
#> $overall
#>   total.num.ind crt.num.ind crt.num.haplo acc.ind acc.haplo call.threshold
#> 1             0           0             0     NaN       NaN              0
#>   n.call call.rate
#> 1      0       NaN
#> 
#> $confusion
#>                     True
#> Predict              KIR2DL2*00101 KIR2DL2*00301 KIR2DL2*005 KIR2DL2*00602
#>   KIR2DL2*00101                  0             0           0             0
#>   KIR2DL2*00301                  0             0           0             0
#>   KIR2DL2*005                    0             0           0             0
#>   KIR2DL2*00602                  0             0           0             0
#>   KIR2DL2*new                    0             0           0             0
#>   KIR2DL2*null                   0             0           0             0
#>   KIR2DL23*0000                  0             0           0             0
#>   KIR2DL3*00101                  0             0           0             0
#>   KIR2DL3*00102                  0             0           0             0
#>   KIR2DL3*00110                  0             0           0             0
#>   KIR2DL3*00201                  0             0           0             0
#>   KIR2DL3*003                    0             0           0             0
#>   KIR2DL3*00501                  0             0           0             0
#>   KIR2DL3*006                    0             0           0             0
#>   KIR2DL3*009                    0             0           0             0
#>   KIR2DL3*011                    0             0           0             0
#>   KIR2DL3*01202                  0             0           0             0
#>   KIR2DL3*013                    0             0           0             0
#>   KIR2DL3*015                    0             0           0             0
#>   KIR2DL3*019                    0             0           0             0
#>   KIR2DL3*022                    0             0           0             0
#>   KIR2DL3*null                   0             0           0             0
#>   KIR2DL3*unresolved             0             0           0             0
#>   ...                            0             0           0             0
#>                     True
#> Predict              KIR2DL2*new KIR2DL2*null KIR2DL23*0000 KIR2DL3*00101
#>   KIR2DL2*00101                0            0             0             0
#>   KIR2DL2*00301                0            0             0             0
#>   KIR2DL2*005                  0            0             0             0
#>   KIR2DL2*00602                0            0             0             0
#>   KIR2DL2*new                  0            0             0             0
#>   KIR2DL2*null                 0            0             0             0
#>   KIR2DL23*0000                0            0             0             0
#>   KIR2DL3*00101                0            0             0             0
#>   KIR2DL3*00102                0            0             0             0
#>   KIR2DL3*00110                0            0             0             0
#>   KIR2DL3*00201                0            0             0             0
#>   KIR2DL3*003                  0            0             0             0
#>   KIR2DL3*00501                0            0             0             0
#>   KIR2DL3*006                  0            0             0             0
#>   KIR2DL3*009                  0            0             0             0
#>   KIR2DL3*011                  0            0             0             0
#>   KIR2DL3*01202                0            0             0             0
#>   KIR2DL3*013                  0            0             0             0
#>   KIR2DL3*015                  0            0             0             0
#>   KIR2DL3*019                  0            0             0             0
#>   KIR2DL3*022                  0            0             0             0
#>   KIR2DL3*null                 0            0             0             0
#>   KIR2DL3*unresolved           0            0             0             0
#>   ...                          0            0             0             0
#>                     True
#> Predict              KIR2DL3*00102 KIR2DL3*00110 KIR2DL3*00201 KIR2DL3*003
#>   KIR2DL2*00101                  0             0             0           0
#>   KIR2DL2*00301                  0             0             0           0
#>   KIR2DL2*005                    0             0             0           0
#>   KIR2DL2*00602                  0             0             0           0
#>   KIR2DL2*new                    0             0             0           0
#>   KIR2DL2*null                   0             0             0           0
#>   KIR2DL23*0000                  0             0             0           0
#>   KIR2DL3*00101                  0             0             0           0
#>   KIR2DL3*00102                  0             0             0           0
#>   KIR2DL3*00110                  0             0             0           0
#>   KIR2DL3*00201                  0             0             0           0
#>   KIR2DL3*003                    0             0             0           0
#>   KIR2DL3*00501                  0             0             0           0
#>   KIR2DL3*006                    0             0             0           0
#>   KIR2DL3*009                    0             0             0           0
#>   KIR2DL3*011                    0             0             0           0
#>   KIR2DL3*01202                  0             0             0           0
#>   KIR2DL3*013                    0             0             0           0
#>   KIR2DL3*015                    0             0             0           0
#>   KIR2DL3*019                    0             0             0           0
#>   KIR2DL3*022                    0             0             0           0
#>   KIR2DL3*null                   0             0             0           0
#>   KIR2DL3*unresolved             0             0             0           0
#>   ...                            0             0             0           0
#>                     True
#> Predict              KIR2DL3*00501 KIR2DL3*006 KIR2DL3*009 KIR2DL3*011
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*01202 KIR2DL3*013 KIR2DL3*015 KIR2DL3*019
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*022 KIR2DL3*null KIR2DL3*unresolved
#>   KIR2DL2*00101                0            0                  0
#>   KIR2DL2*00301                0            0                  0
#>   KIR2DL2*005                  0            0                  0
#>   KIR2DL2*00602                0            0                  0
#>   KIR2DL2*new                  0            0                  0
#>   KIR2DL2*null                 0            0                  0
#>   KIR2DL23*0000                0            0                  0
#>   KIR2DL3*00101                0            0                  0
#>   KIR2DL3*00102                0            0                  0
#>   KIR2DL3*00110                0            0                  0
#>   KIR2DL3*00201                0            0                  0
#>   KIR2DL3*003                  0            0                  0
#>   KIR2DL3*00501                0            0                  0
#>   KIR2DL3*006                  0            0                  0
#>   KIR2DL3*009                  0            0                  0
#>   KIR2DL3*011                  0            0                  0
#>   KIR2DL3*01202                0            0                  0
#>   KIR2DL3*013                  0            0                  0
#>   KIR2DL3*015                  0            0                  0
#>   KIR2DL3*019                  0            0                  0
#>   KIR2DL3*022                  0            0                  0
#>   KIR2DL3*null                 0            0                  0
#>   KIR2DL3*unresolved           0            0                  0
#>   ...                          0            0                  0
#> 
#> $detail
#>                allele train.num   train.freq valid.num valid.freq call.rate
#> 1       KIR2DL2*00101       176 0.1494057725         0        NaN         0
#> 2       KIR2DL2*00301       103 0.0874363328         0        NaN         0
#> 3         KIR2DL2*005         1 0.0008488964         0        NaN         0
#> 4       KIR2DL2*00602        16 0.0135823430         0        NaN         0
#> 5         KIR2DL2*new         1 0.0008488964         0        NaN         0
#> 6        KIR2DL2*null         2 0.0016977929         0        NaN         0
#> 7       KIR2DL23*0000        33 0.0280135823         0        NaN         0
#> 8       KIR2DL3*00101       560 0.4753820034         0        NaN         0
#> 9       KIR2DL3*00102         4 0.0033955857         0        NaN         0
#> 10      KIR2DL3*00110         4 0.0033955857         0        NaN         0
#> 11      KIR2DL3*00201       169 0.1434634975         0        NaN         0
#> 12        KIR2DL3*003         1 0.0008488964         0        NaN         0
#> 13      KIR2DL3*00501        55 0.0466893039         0        NaN         0
#> 14        KIR2DL3*006        12 0.0101867572         0        NaN         0
#> 15        KIR2DL3*009         7 0.0059422750         0        NaN         0
#> 16        KIR2DL3*011         3 0.0025466893         0        NaN         0
#> 17      KIR2DL3*01202         9 0.0076400679         0        NaN         0
#> 18        KIR2DL3*013         5 0.0042444822         0        NaN         0
#> 19        KIR2DL3*015         1 0.0008488964         0        NaN         0
#> 20        KIR2DL3*019         3 0.0025466893         0        NaN         0
#> 21        KIR2DL3*022         1 0.0008488964         0        NaN         0
#> 22       KIR2DL3*null         9 0.0076400679         0        NaN         0
#> 23 KIR2DL3*unresolved         3 0.0025466893         0        NaN         0
#>    accuracy sensitivity specificity ppv npv miscall miscall.prop
#> 1       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 2       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 3       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 4       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 5       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 6       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 7       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 8       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 9       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 10      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 11      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 12      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 13      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 14      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 15      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 16      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 17      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 18      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 19      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 20      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 21      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 22      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 23      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
  call.threshold=0.5))
#> Calling 'hlaCompareAllele': there are 0 individuals in common.
#> $overall
#>   total.num.ind crt.num.ind crt.num.haplo acc.ind acc.haplo call.threshold
#> 1             0           0             0     NaN       NaN            0.5
#>   n.call call.rate
#> 1      0       NaN
#> 
#> $confusion
#>                     True
#> Predict              KIR2DL2*00101 KIR2DL2*00301 KIR2DL2*005 KIR2DL2*00602
#>   KIR2DL2*00101                  0             0           0             0
#>   KIR2DL2*00301                  0             0           0             0
#>   KIR2DL2*005                    0             0           0             0
#>   KIR2DL2*00602                  0             0           0             0
#>   KIR2DL2*new                    0             0           0             0
#>   KIR2DL2*null                   0             0           0             0
#>   KIR2DL23*0000                  0             0           0             0
#>   KIR2DL3*00101                  0             0           0             0
#>   KIR2DL3*00102                  0             0           0             0
#>   KIR2DL3*00110                  0             0           0             0
#>   KIR2DL3*00201                  0             0           0             0
#>   KIR2DL3*003                    0             0           0             0
#>   KIR2DL3*00501                  0             0           0             0
#>   KIR2DL3*006                    0             0           0             0
#>   KIR2DL3*009                    0             0           0             0
#>   KIR2DL3*011                    0             0           0             0
#>   KIR2DL3*01202                  0             0           0             0
#>   KIR2DL3*013                    0             0           0             0
#>   KIR2DL3*015                    0             0           0             0
#>   KIR2DL3*019                    0             0           0             0
#>   KIR2DL3*022                    0             0           0             0
#>   KIR2DL3*null                   0             0           0             0
#>   KIR2DL3*unresolved             0             0           0             0
#>   ...                            0             0           0             0
#>                     True
#> Predict              KIR2DL2*new KIR2DL2*null KIR2DL23*0000 KIR2DL3*00101
#>   KIR2DL2*00101                0            0             0             0
#>   KIR2DL2*00301                0            0             0             0
#>   KIR2DL2*005                  0            0             0             0
#>   KIR2DL2*00602                0            0             0             0
#>   KIR2DL2*new                  0            0             0             0
#>   KIR2DL2*null                 0            0             0             0
#>   KIR2DL23*0000                0            0             0             0
#>   KIR2DL3*00101                0            0             0             0
#>   KIR2DL3*00102                0            0             0             0
#>   KIR2DL3*00110                0            0             0             0
#>   KIR2DL3*00201                0            0             0             0
#>   KIR2DL3*003                  0            0             0             0
#>   KIR2DL3*00501                0            0             0             0
#>   KIR2DL3*006                  0            0             0             0
#>   KIR2DL3*009                  0            0             0             0
#>   KIR2DL3*011                  0            0             0             0
#>   KIR2DL3*01202                0            0             0             0
#>   KIR2DL3*013                  0            0             0             0
#>   KIR2DL3*015                  0            0             0             0
#>   KIR2DL3*019                  0            0             0             0
#>   KIR2DL3*022                  0            0             0             0
#>   KIR2DL3*null                 0            0             0             0
#>   KIR2DL3*unresolved           0            0             0             0
#>   ...                          0            0             0             0
#>                     True
#> Predict              KIR2DL3*00102 KIR2DL3*00110 KIR2DL3*00201 KIR2DL3*003
#>   KIR2DL2*00101                  0             0             0           0
#>   KIR2DL2*00301                  0             0             0           0
#>   KIR2DL2*005                    0             0             0           0
#>   KIR2DL2*00602                  0             0             0           0
#>   KIR2DL2*new                    0             0             0           0
#>   KIR2DL2*null                   0             0             0           0
#>   KIR2DL23*0000                  0             0             0           0
#>   KIR2DL3*00101                  0             0             0           0
#>   KIR2DL3*00102                  0             0             0           0
#>   KIR2DL3*00110                  0             0             0           0
#>   KIR2DL3*00201                  0             0             0           0
#>   KIR2DL3*003                    0             0             0           0
#>   KIR2DL3*00501                  0             0             0           0
#>   KIR2DL3*006                    0             0             0           0
#>   KIR2DL3*009                    0             0             0           0
#>   KIR2DL3*011                    0             0             0           0
#>   KIR2DL3*01202                  0             0             0           0
#>   KIR2DL3*013                    0             0             0           0
#>   KIR2DL3*015                    0             0             0           0
#>   KIR2DL3*019                    0             0             0           0
#>   KIR2DL3*022                    0             0             0           0
#>   KIR2DL3*null                   0             0             0           0
#>   KIR2DL3*unresolved             0             0             0           0
#>   ...                            0             0             0           0
#>                     True
#> Predict              KIR2DL3*00501 KIR2DL3*006 KIR2DL3*009 KIR2DL3*011
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*01202 KIR2DL3*013 KIR2DL3*015 KIR2DL3*019
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*022 KIR2DL3*null KIR2DL3*unresolved
#>   KIR2DL2*00101                0            0                  0
#>   KIR2DL2*00301                0            0                  0
#>   KIR2DL2*005                  0            0                  0
#>   KIR2DL2*00602                0            0                  0
#>   KIR2DL2*new                  0            0                  0
#>   KIR2DL2*null                 0            0                  0
#>   KIR2DL23*0000                0            0                  0
#>   KIR2DL3*00101                0            0                  0
#>   KIR2DL3*00102                0            0                  0
#>   KIR2DL3*00110                0            0                  0
#>   KIR2DL3*00201                0            0                  0
#>   KIR2DL3*003                  0            0                  0
#>   KIR2DL3*00501                0            0                  0
#>   KIR2DL3*006                  0            0                  0
#>   KIR2DL3*009                  0            0                  0
#>   KIR2DL3*011                  0            0                  0
#>   KIR2DL3*01202                0            0                  0
#>   KIR2DL3*013                  0            0                  0
#>   KIR2DL3*015                  0            0                  0
#>   KIR2DL3*019                  0            0                  0
#>   KIR2DL3*022                  0            0                  0
#>   KIR2DL3*null                 0            0                  0
#>   KIR2DL3*unresolved           0            0                  0
#>   ...                          0            0                  0
#> 
#> $detail
#>                allele train.num   train.freq valid.num valid.freq call.rate
#> 1       KIR2DL2*00101       176 0.1494057725         0        NaN         0
#> 2       KIR2DL2*00301       103 0.0874363328         0        NaN         0
#> 3         KIR2DL2*005         1 0.0008488964         0        NaN         0
#> 4       KIR2DL2*00602        16 0.0135823430         0        NaN         0
#> 5         KIR2DL2*new         1 0.0008488964         0        NaN         0
#> 6        KIR2DL2*null         2 0.0016977929         0        NaN         0
#> 7       KIR2DL23*0000        33 0.0280135823         0        NaN         0
#> 8       KIR2DL3*00101       560 0.4753820034         0        NaN         0
#> 9       KIR2DL3*00102         4 0.0033955857         0        NaN         0
#> 10      KIR2DL3*00110         4 0.0033955857         0        NaN         0
#> 11      KIR2DL3*00201       169 0.1434634975         0        NaN         0
#> 12        KIR2DL3*003         1 0.0008488964         0        NaN         0
#> 13      KIR2DL3*00501        55 0.0466893039         0        NaN         0
#> 14        KIR2DL3*006        12 0.0101867572         0        NaN         0
#> 15        KIR2DL3*009         7 0.0059422750         0        NaN         0
#> 16        KIR2DL3*011         3 0.0025466893         0        NaN         0
#> 17      KIR2DL3*01202         9 0.0076400679         0        NaN         0
#> 18        KIR2DL3*013         5 0.0042444822         0        NaN         0
#> 19        KIR2DL3*015         1 0.0008488964         0        NaN         0
#> 20        KIR2DL3*019         3 0.0025466893         0        NaN         0
#> 21        KIR2DL3*022         1 0.0008488964         0        NaN         0
#> 22       KIR2DL3*null         9 0.0076400679         0        NaN         0
#> 23 KIR2DL3*unresolved         3 0.0025466893         0        NaN         0
#>    accuracy sensitivity specificity ppv npv miscall miscall.prop
#> 1       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 2       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 3       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 4       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 5       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 6       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 7       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 8       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 9       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 10      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 11      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 12      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 13      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 14      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 15      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 16      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 17      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 18      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 19      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 20      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 21      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 22      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 23      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 


# save the parameter file
mobj <- hlaModelToObj(model)
save(mobj, file="HIBAG_model.RData")
save(test.geno, file="testgeno.RData")
save(hlatab, file="HLASplit.RData")

# Clear Workspace
hlaClose(model)  # release all resources of model
rm(list = ls())


######################################################################

# NOW, load a HIBAG model from the parameter file
mobj <- get(load("HIBAG_model.RData"))
model <- hlaModelFromObj(mobj)

# validation
test.geno <- get(load("testgeno.RData"))
hlatab <- get(load("HLASplit.RData"))

pred <- predict(model, test.geno, type="response")
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
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
  call.threshold=0))
#> Calling 'hlaCompareAllele': there are 0 individuals in common.
#> $overall
#>   total.num.ind crt.num.ind crt.num.haplo acc.ind acc.haplo call.threshold
#> 1             0           0             0     NaN       NaN              0
#>   n.call call.rate
#> 1      0       NaN
#> 
#> $confusion
#>                     True
#> Predict              KIR2DL2*00101 KIR2DL2*00301 KIR2DL2*005 KIR2DL2*00602
#>   KIR2DL2*00101                  0             0           0             0
#>   KIR2DL2*00301                  0             0           0             0
#>   KIR2DL2*005                    0             0           0             0
#>   KIR2DL2*00602                  0             0           0             0
#>   KIR2DL2*new                    0             0           0             0
#>   KIR2DL2*null                   0             0           0             0
#>   KIR2DL23*0000                  0             0           0             0
#>   KIR2DL3*00101                  0             0           0             0
#>   KIR2DL3*00102                  0             0           0             0
#>   KIR2DL3*00110                  0             0           0             0
#>   KIR2DL3*00201                  0             0           0             0
#>   KIR2DL3*003                    0             0           0             0
#>   KIR2DL3*00501                  0             0           0             0
#>   KIR2DL3*006                    0             0           0             0
#>   KIR2DL3*009                    0             0           0             0
#>   KIR2DL3*011                    0             0           0             0
#>   KIR2DL3*01202                  0             0           0             0
#>   KIR2DL3*013                    0             0           0             0
#>   KIR2DL3*015                    0             0           0             0
#>   KIR2DL3*019                    0             0           0             0
#>   KIR2DL3*022                    0             0           0             0
#>   KIR2DL3*null                   0             0           0             0
#>   KIR2DL3*unresolved             0             0           0             0
#>   ...                            0             0           0             0
#>                     True
#> Predict              KIR2DL2*new KIR2DL2*null KIR2DL23*0000 KIR2DL3*00101
#>   KIR2DL2*00101                0            0             0             0
#>   KIR2DL2*00301                0            0             0             0
#>   KIR2DL2*005                  0            0             0             0
#>   KIR2DL2*00602                0            0             0             0
#>   KIR2DL2*new                  0            0             0             0
#>   KIR2DL2*null                 0            0             0             0
#>   KIR2DL23*0000                0            0             0             0
#>   KIR2DL3*00101                0            0             0             0
#>   KIR2DL3*00102                0            0             0             0
#>   KIR2DL3*00110                0            0             0             0
#>   KIR2DL3*00201                0            0             0             0
#>   KIR2DL3*003                  0            0             0             0
#>   KIR2DL3*00501                0            0             0             0
#>   KIR2DL3*006                  0            0             0             0
#>   KIR2DL3*009                  0            0             0             0
#>   KIR2DL3*011                  0            0             0             0
#>   KIR2DL3*01202                0            0             0             0
#>   KIR2DL3*013                  0            0             0             0
#>   KIR2DL3*015                  0            0             0             0
#>   KIR2DL3*019                  0            0             0             0
#>   KIR2DL3*022                  0            0             0             0
#>   KIR2DL3*null                 0            0             0             0
#>   KIR2DL3*unresolved           0            0             0             0
#>   ...                          0            0             0             0
#>                     True
#> Predict              KIR2DL3*00102 KIR2DL3*00110 KIR2DL3*00201 KIR2DL3*003
#>   KIR2DL2*00101                  0             0             0           0
#>   KIR2DL2*00301                  0             0             0           0
#>   KIR2DL2*005                    0             0             0           0
#>   KIR2DL2*00602                  0             0             0           0
#>   KIR2DL2*new                    0             0             0           0
#>   KIR2DL2*null                   0             0             0           0
#>   KIR2DL23*0000                  0             0             0           0
#>   KIR2DL3*00101                  0             0             0           0
#>   KIR2DL3*00102                  0             0             0           0
#>   KIR2DL3*00110                  0             0             0           0
#>   KIR2DL3*00201                  0             0             0           0
#>   KIR2DL3*003                    0             0             0           0
#>   KIR2DL3*00501                  0             0             0           0
#>   KIR2DL3*006                    0             0             0           0
#>   KIR2DL3*009                    0             0             0           0
#>   KIR2DL3*011                    0             0             0           0
#>   KIR2DL3*01202                  0             0             0           0
#>   KIR2DL3*013                    0             0             0           0
#>   KIR2DL3*015                    0             0             0           0
#>   KIR2DL3*019                    0             0             0           0
#>   KIR2DL3*022                    0             0             0           0
#>   KIR2DL3*null                   0             0             0           0
#>   KIR2DL3*unresolved             0             0             0           0
#>   ...                            0             0             0           0
#>                     True
#> Predict              KIR2DL3*00501 KIR2DL3*006 KIR2DL3*009 KIR2DL3*011
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*01202 KIR2DL3*013 KIR2DL3*015 KIR2DL3*019
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*022 KIR2DL3*null KIR2DL3*unresolved
#>   KIR2DL2*00101                0            0                  0
#>   KIR2DL2*00301                0            0                  0
#>   KIR2DL2*005                  0            0                  0
#>   KIR2DL2*00602                0            0                  0
#>   KIR2DL2*new                  0            0                  0
#>   KIR2DL2*null                 0            0                  0
#>   KIR2DL23*0000                0            0                  0
#>   KIR2DL3*00101                0            0                  0
#>   KIR2DL3*00102                0            0                  0
#>   KIR2DL3*00110                0            0                  0
#>   KIR2DL3*00201                0            0                  0
#>   KIR2DL3*003                  0            0                  0
#>   KIR2DL3*00501                0            0                  0
#>   KIR2DL3*006                  0            0                  0
#>   KIR2DL3*009                  0            0                  0
#>   KIR2DL3*011                  0            0                  0
#>   KIR2DL3*01202                0            0                  0
#>   KIR2DL3*013                  0            0                  0
#>   KIR2DL3*015                  0            0                  0
#>   KIR2DL3*019                  0            0                  0
#>   KIR2DL3*022                  0            0                  0
#>   KIR2DL3*null                 0            0                  0
#>   KIR2DL3*unresolved           0            0                  0
#>   ...                          0            0                  0
#> 
#> $detail
#>                allele train.num   train.freq valid.num valid.freq call.rate
#> 1       KIR2DL2*00101       176 0.1494057725         0        NaN         0
#> 2       KIR2DL2*00301       103 0.0874363328         0        NaN         0
#> 3         KIR2DL2*005         1 0.0008488964         0        NaN         0
#> 4       KIR2DL2*00602        16 0.0135823430         0        NaN         0
#> 5         KIR2DL2*new         1 0.0008488964         0        NaN         0
#> 6        KIR2DL2*null         2 0.0016977929         0        NaN         0
#> 7       KIR2DL23*0000        33 0.0280135823         0        NaN         0
#> 8       KIR2DL3*00101       560 0.4753820034         0        NaN         0
#> 9       KIR2DL3*00102         4 0.0033955857         0        NaN         0
#> 10      KIR2DL3*00110         4 0.0033955857         0        NaN         0
#> 11      KIR2DL3*00201       169 0.1434634975         0        NaN         0
#> 12        KIR2DL3*003         1 0.0008488964         0        NaN         0
#> 13      KIR2DL3*00501        55 0.0466893039         0        NaN         0
#> 14        KIR2DL3*006        12 0.0101867572         0        NaN         0
#> 15        KIR2DL3*009         7 0.0059422750         0        NaN         0
#> 16        KIR2DL3*011         3 0.0025466893         0        NaN         0
#> 17      KIR2DL3*01202         9 0.0076400679         0        NaN         0
#> 18        KIR2DL3*013         5 0.0042444822         0        NaN         0
#> 19        KIR2DL3*015         1 0.0008488964         0        NaN         0
#> 20        KIR2DL3*019         3 0.0025466893         0        NaN         0
#> 21        KIR2DL3*022         1 0.0008488964         0        NaN         0
#> 22       KIR2DL3*null         9 0.0076400679         0        NaN         0
#> 23 KIR2DL3*unresolved         3 0.0025466893         0        NaN         0
#>    accuracy sensitivity specificity ppv npv miscall miscall.prop
#> 1       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 2       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 3       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 4       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 5       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 6       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 7       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 8       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 9       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 10      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 11      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 12      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 13      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 14      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 15      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 16      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 17      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 18      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 19      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 20      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 21      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 22      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 23      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
  call.threshold=0.5))
#> Calling 'hlaCompareAllele': there are 0 individuals in common.
#> $overall
#>   total.num.ind crt.num.ind crt.num.haplo acc.ind acc.haplo call.threshold
#> 1             0           0             0     NaN       NaN            0.5
#>   n.call call.rate
#> 1      0       NaN
#> 
#> $confusion
#>                     True
#> Predict              KIR2DL2*00101 KIR2DL2*00301 KIR2DL2*005 KIR2DL2*00602
#>   KIR2DL2*00101                  0             0           0             0
#>   KIR2DL2*00301                  0             0           0             0
#>   KIR2DL2*005                    0             0           0             0
#>   KIR2DL2*00602                  0             0           0             0
#>   KIR2DL2*new                    0             0           0             0
#>   KIR2DL2*null                   0             0           0             0
#>   KIR2DL23*0000                  0             0           0             0
#>   KIR2DL3*00101                  0             0           0             0
#>   KIR2DL3*00102                  0             0           0             0
#>   KIR2DL3*00110                  0             0           0             0
#>   KIR2DL3*00201                  0             0           0             0
#>   KIR2DL3*003                    0             0           0             0
#>   KIR2DL3*00501                  0             0           0             0
#>   KIR2DL3*006                    0             0           0             0
#>   KIR2DL3*009                    0             0           0             0
#>   KIR2DL3*011                    0             0           0             0
#>   KIR2DL3*01202                  0             0           0             0
#>   KIR2DL3*013                    0             0           0             0
#>   KIR2DL3*015                    0             0           0             0
#>   KIR2DL3*019                    0             0           0             0
#>   KIR2DL3*022                    0             0           0             0
#>   KIR2DL3*null                   0             0           0             0
#>   KIR2DL3*unresolved             0             0           0             0
#>   ...                            0             0           0             0
#>                     True
#> Predict              KIR2DL2*new KIR2DL2*null KIR2DL23*0000 KIR2DL3*00101
#>   KIR2DL2*00101                0            0             0             0
#>   KIR2DL2*00301                0            0             0             0
#>   KIR2DL2*005                  0            0             0             0
#>   KIR2DL2*00602                0            0             0             0
#>   KIR2DL2*new                  0            0             0             0
#>   KIR2DL2*null                 0            0             0             0
#>   KIR2DL23*0000                0            0             0             0
#>   KIR2DL3*00101                0            0             0             0
#>   KIR2DL3*00102                0            0             0             0
#>   KIR2DL3*00110                0            0             0             0
#>   KIR2DL3*00201                0            0             0             0
#>   KIR2DL3*003                  0            0             0             0
#>   KIR2DL3*00501                0            0             0             0
#>   KIR2DL3*006                  0            0             0             0
#>   KIR2DL3*009                  0            0             0             0
#>   KIR2DL3*011                  0            0             0             0
#>   KIR2DL3*01202                0            0             0             0
#>   KIR2DL3*013                  0            0             0             0
#>   KIR2DL3*015                  0            0             0             0
#>   KIR2DL3*019                  0            0             0             0
#>   KIR2DL3*022                  0            0             0             0
#>   KIR2DL3*null                 0            0             0             0
#>   KIR2DL3*unresolved           0            0             0             0
#>   ...                          0            0             0             0
#>                     True
#> Predict              KIR2DL3*00102 KIR2DL3*00110 KIR2DL3*00201 KIR2DL3*003
#>   KIR2DL2*00101                  0             0             0           0
#>   KIR2DL2*00301                  0             0             0           0
#>   KIR2DL2*005                    0             0             0           0
#>   KIR2DL2*00602                  0             0             0           0
#>   KIR2DL2*new                    0             0             0           0
#>   KIR2DL2*null                   0             0             0           0
#>   KIR2DL23*0000                  0             0             0           0
#>   KIR2DL3*00101                  0             0             0           0
#>   KIR2DL3*00102                  0             0             0           0
#>   KIR2DL3*00110                  0             0             0           0
#>   KIR2DL3*00201                  0             0             0           0
#>   KIR2DL3*003                    0             0             0           0
#>   KIR2DL3*00501                  0             0             0           0
#>   KIR2DL3*006                    0             0             0           0
#>   KIR2DL3*009                    0             0             0           0
#>   KIR2DL3*011                    0             0             0           0
#>   KIR2DL3*01202                  0             0             0           0
#>   KIR2DL3*013                    0             0             0           0
#>   KIR2DL3*015                    0             0             0           0
#>   KIR2DL3*019                    0             0             0           0
#>   KIR2DL3*022                    0             0             0           0
#>   KIR2DL3*null                   0             0             0           0
#>   KIR2DL3*unresolved             0             0             0           0
#>   ...                            0             0             0           0
#>                     True
#> Predict              KIR2DL3*00501 KIR2DL3*006 KIR2DL3*009 KIR2DL3*011
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*01202 KIR2DL3*013 KIR2DL3*015 KIR2DL3*019
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*022 KIR2DL3*null KIR2DL3*unresolved
#>   KIR2DL2*00101                0            0                  0
#>   KIR2DL2*00301                0            0                  0
#>   KIR2DL2*005                  0            0                  0
#>   KIR2DL2*00602                0            0                  0
#>   KIR2DL2*new                  0            0                  0
#>   KIR2DL2*null                 0            0                  0
#>   KIR2DL23*0000                0            0                  0
#>   KIR2DL3*00101                0            0                  0
#>   KIR2DL3*00102                0            0                  0
#>   KIR2DL3*00110                0            0                  0
#>   KIR2DL3*00201                0            0                  0
#>   KIR2DL3*003                  0            0                  0
#>   KIR2DL3*00501                0            0                  0
#>   KIR2DL3*006                  0            0                  0
#>   KIR2DL3*009                  0            0                  0
#>   KIR2DL3*011                  0            0                  0
#>   KIR2DL3*01202                0            0                  0
#>   KIR2DL3*013                  0            0                  0
#>   KIR2DL3*015                  0            0                  0
#>   KIR2DL3*019                  0            0                  0
#>   KIR2DL3*022                  0            0                  0
#>   KIR2DL3*null                 0            0                  0
#>   KIR2DL3*unresolved           0            0                  0
#>   ...                          0            0                  0
#> 
#> $detail
#>                allele train.num   train.freq valid.num valid.freq call.rate
#> 1       KIR2DL2*00101       176 0.1494057725         0        NaN         0
#> 2       KIR2DL2*00301       103 0.0874363328         0        NaN         0
#> 3         KIR2DL2*005         1 0.0008488964         0        NaN         0
#> 4       KIR2DL2*00602        16 0.0135823430         0        NaN         0
#> 5         KIR2DL2*new         1 0.0008488964         0        NaN         0
#> 6        KIR2DL2*null         2 0.0016977929         0        NaN         0
#> 7       KIR2DL23*0000        33 0.0280135823         0        NaN         0
#> 8       KIR2DL3*00101       560 0.4753820034         0        NaN         0
#> 9       KIR2DL3*00102         4 0.0033955857         0        NaN         0
#> 10      KIR2DL3*00110         4 0.0033955857         0        NaN         0
#> 11      KIR2DL3*00201       169 0.1434634975         0        NaN         0
#> 12        KIR2DL3*003         1 0.0008488964         0        NaN         0
#> 13      KIR2DL3*00501        55 0.0466893039         0        NaN         0
#> 14        KIR2DL3*006        12 0.0101867572         0        NaN         0
#> 15        KIR2DL3*009         7 0.0059422750         0        NaN         0
#> 16        KIR2DL3*011         3 0.0025466893         0        NaN         0
#> 17      KIR2DL3*01202         9 0.0076400679         0        NaN         0
#> 18        KIR2DL3*013         5 0.0042444822         0        NaN         0
#> 19        KIR2DL3*015         1 0.0008488964         0        NaN         0
#> 20        KIR2DL3*019         3 0.0025466893         0        NaN         0
#> 21        KIR2DL3*022         1 0.0008488964         0        NaN         0
#> 22       KIR2DL3*null         9 0.0076400679         0        NaN         0
#> 23 KIR2DL3*unresolved         3 0.0025466893         0        NaN         0
#>    accuracy sensitivity specificity ppv npv miscall miscall.prop
#> 1       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 2       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 3       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 4       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 5       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 6       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 7       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 8       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 9       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 10      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 11      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 12      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 13      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 14      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 15      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 16      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 17      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 18      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 19      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 20      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 21      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 22      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 23      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 
```
