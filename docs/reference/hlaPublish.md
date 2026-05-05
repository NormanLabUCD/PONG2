# Finalize a HIBAG model

Finalize a HIBAG model by removing unused SNP predictors and adding
appendix information (platform, training set, authors, warning, etc)

## Usage

``` r
hlaPublish(mobj, platform=NULL, information=NULL, warning=NULL,
  rm.unused.snp=TRUE, anonymize=TRUE, verbose=TRUE)
```

## Arguments

- mobj:

  an object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)
  or
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)

- platform:

  the text of platform information

- information:

  the other information, like authors

- warning:

  any warning message

- rm.unused.snp:

  if `TRUE`, remove unused SNPs from the model

- anonymize:

  if `TRUE`, remove sample IDs

- verbose:

  if TRUE, show information

## Value

Returns a new object of
[`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md).

## Author

Xiuwen Zheng

## See also

[`hlaModelFromObj`](https://normanlabucd.github.io/PONG2/reference/hlaModelfromObj.md),
[`hlaModelToObj`](https://normanlabucd.github.io/PONG2/reference/hlaModelToObj.md)

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

# training genotypes
region <- 250   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel = match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'snpid' not found


#
# train a HIBAG model
#
set.seed(1000)

# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
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
length(model$snp.id)
#> [1] 2770

mobj <- hlaPublish(model,
  platform = "Illumina 1M Duo",
  information = "Training set -- HapMap Phase II")
#> Remove 678 unused SNPs.
model2 <- hlaModelFromObj(mobj)
length(mobj$snp.id)
#> [1] 2092
mobj$appendix
#> $platform
#> [1] "Illumina 1M Duo"
#> 
#> $information
#> [1] "Training set -- HapMap Phase II"
#> 
#> $warning
#> NULL
#> 
summary(mobj)
#> Gene: KIR2DL23
#> Training dataset: 589 samples X 2092 SNPs
#>  # of KIR3DL1/S1 alleles: 23
#>  # of individual classifiers: 100
#>  total # of SNPs used: 2092
#>  average # of SNPs in an individual classifier: 103.39, sd: 14.99, min: 63, max: 128
#>  average # of haplotypes in an individual classifier: 421.30, sd: 62.46, min: 288, max: 614
#>  average out-of-bag accuracy: 95.74%, sd: 0.89%, min: 93.09%, max: 98.39%
#> Genome assembly: hg38
#> Platform: Illumina 1M Duo 
#> Information: Training set -- HapMap Phase II 

p1 <- predict(model, train.geno)
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"
p2 <- predict(model2, train.geno)
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"

# check
cbind(p1$value, p2$value)
#> Error: object 'p1' not found
```
