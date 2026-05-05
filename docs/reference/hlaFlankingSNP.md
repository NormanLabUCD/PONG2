# SNP IDs in Flanking Region

To select SNPs in the flanking region of a specified HLA locus.

## Usage

``` r
hlaFlankingSNP(snp.id, position, hla.id, flank.bp=500*1000,
  assembly=c("auto", "hg18", "hg19", "unknown"))
```

## Arguments

- snp.id:

  a vector of SNP IDs

- position:

  a vector of positions

- hla.id:

  the name of HLA locus

- flank.bp:

  the size of flanking region on each side in basepair

- assembly:

  the human genome reference: "hg19" (default), "hg18", "auto" refers to
  "hg19"

## Details

`hla.id` is "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1", "DPB1" or
"any".

## Value

Return allele frequecies.

## Author

Xiuwen Zheng

## See also

[`hlaGenoSubset`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSubset.md),
[`hlaLociInfo`](https://normanlabucd.github.io/PONG2/reference/hlaLociInfo.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "KIR3DLS1"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")
#> Error in `[.data.frame`(HLA_Type_Table, , paste(hla.id, ".1", sep = "")): undefined columns selected

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"

train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel  = match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel = match(hla$training$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'hla' not found
summary(train.geno)
#> Error: object 'train.geno' not found
```
