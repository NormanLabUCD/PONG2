# Convert to PLINK PED format

Convert an object of
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)
to a file of PLINK PED format.

## Usage

``` r
hlaGeno2PED(geno, out.fn)
```

## Arguments

- geno:

  a genotype object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- out.fn:

  the file name of output ped file

## Details

Two files ".map" and ".ped" are created.

## Value

None.

## Author

Xiuwen Zheng

## See also

[`hlaBED2Geno`](https://normanlabucd.github.io/PONG2/reference/hlaBED2Geno.md)

## Examples

``` r
# load SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  max.resolution=4, locus=hla.id, assembly="hg19")

# training genotypes
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"

train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id))
#> Error: object 'snpid' not found

hlaGeno2PED(train.geno, "test")
#> Error: object 'train.geno' not found
```
