# Plot a "hlaAttrBagObj" object

To show a scatterplot of the numbers of individual classifiers and SNP
positions.

## Usage

``` r
# S3 method for class 'hlaAttrBagObj'
plot(x, xlab=NULL, ylab=NULL,
  locus.color="red", locus.lty=2, locus.cex=1.25,
  assembly=c("auto", "hg18", "hg19", "unknown"), ...)
```

## Arguments

- x:

  an object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- xlab:

  the label of X-axis

- ylab:

  the label of Y-axis

- locus.color:

  the color of text and line for HLA locus

- locus.lty:

  the type of line for HLA locus

- locus.cex:

  the font size of HLA locus

- assembly:

  the human genome reference: "hg19" (default), "hg18", "auto" refers to
  "hg19"

- ...:

  further arguments passed to or from other methods

## Author

Xiuwen Zheng

## See also

[`print.hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/print.hlaAttrBagObj.md),
[`summary.hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/summary.hlaAttrBagObj.md)

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
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))
#> Error: object 'snpid' not found

# train a HIBAG model
set.seed(1000)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, train.geno, nclassifier=2, verbose.detail=TRUE)
#> Error: object 'train.geno' not found
plot(model)
```
