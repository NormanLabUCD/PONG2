# Combine two HIBAG models together

Merge two objects of
[`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)
together, which is useful for building an ensemble model in parallel.

## Usage

``` r
hlaCombineModelObj(obj1, obj2)
```

## Arguments

- obj1:

  an object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- obj2:

  an object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

## Value

Return an object of
[`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md).

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaModelFiles`](https://normanlabucd.github.io/PONG2/reference/hlaModelFiles.md)

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

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"
length(snpid)  # 275
#> Error: object 'snpid' not found

# training genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))
#> Error: object 'snpid' not found

# train a HIBAG model
set.seed(100)
m1 <- hlaAttrBagging(hla, train.geno, nclassifier=1)
#> Error: object 'train.geno' not found
m2 <- hlaAttrBagging(hla, train.geno, nclassifier=1)
#> Error: object 'train.geno' not found

m1.obj <- hlaModelToObj(m1)
#> Error: object 'm1' not found
m2.obj <- hlaModelToObj(m2)
#> Error: object 'm2' not found

m.obj <- hlaCombineModelObj(m1.obj, m2.obj)
#> Error: object 'm1.obj' not found
summary(m.obj)
#> Error: object 'm.obj' not found
```
