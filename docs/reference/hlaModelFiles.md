# Load a model object from files

To load HIBAG models from a list of files, and merge all together.

## Usage

``` r
hlaModelFiles(fn.list, action.missingfile=c("ignore", "stop"), verbose=TRUE)
```

## Arguments

- fn.list:

  a vector of file names

- action.missingfile:

  "ignore", ignore the missing files, by default; "stop", stop if
  missing

- verbose:

  if TRUE, show information

## Value

Return
[`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md).

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaModelToObj`](https://normanlabucd.github.io/PONG2/reference/hlaModelToObj.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")

# training genotypes
region <- 100   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel = match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'snpid' not found

#
# train HIBAG models
#
set.seed(1000)

# please use "nclassifier=100" when you use HIBAG for real data
model1 <- hlaAttrBagging(hla, train.geno, nclassifier=1, verbose.detail=TRUE)
#> Error: object 'train.geno' not found
mobj1 <- hlaModelToObj(model1)
#> Error in hlaModelToObj(model1): The handle of PONG2 model has been closed.
save(mobj1, file="tm1.RData")

model2 <- hlaAttrBagging(hla, train.geno, nclassifier=1, verbose.detail=TRUE)
#> Error: object 'train.geno' not found
mobj2 <- hlaModelToObj(model2)
#> Error in hlaModelToObj(model2): The handle of PONG2 model has been closed.
save(mobj2, file="tm2.RData")

model3 <- hlaAttrBagging(hla, train.geno, nclassifier=1, verbose.detail=TRUE)
#> Error: object 'train.geno' not found
mobj3 <- hlaModelToObj(model3)
#> Error in hlaModelToObj(model3): The handle of PONG2 model has been closed.
save(mobj3, file="tm3.RData")

# load all of mobj1, mobj2 and mobj3
mobj <- hlaModelFiles(c("tm1.RData", "tm2.RData", "tm3.RData"))
#> Error in hlaCombineModelObj(rv, tmp): identical(obj1$snp.id, obj2$snp.id) is not TRUE
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
```
