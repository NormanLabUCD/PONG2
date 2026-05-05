# Combine two datasets of HLA types

Get a subset of HLA types from an object of
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md).

## Usage

``` r
hlaCombineAllele(H1, H2)
```

## Arguments

- H1:

  the first
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
  object

- H2:

  the second
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
  object

## Value

Return
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md).

## Author

Xiuwen Zheng

## See also

[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md),
[`hlaAlleleSubset`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleSubset.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

head(HLA_Type_Table)
#>   sample.id   A.1   A.2   B.1   B.2   C.1   C.2 DQA1.1 DQA1.2 DQB1.1 DQB1.2
#> 1   NA11882 01:01 29:02 15:01 44:03 06:02 16:01  01:02  03:01  03:02  06:02
#> 2   NA11881 03:01 26:01 07:02 07:02 07:02 07:02  01:02  01:02  06:02  06:02
#> 3   NA11993 26:01 29:02 44:03  <NA> 16:01 16:01  01:01  01:02  05:01  06:02
#> 4   NA11992 01:01 02:01 08:01 35:01 04:01 07:01  01:01  05:01  02:01  05:01
#> 5   NA11995 01:01 01:01 08:01 57:01 06:02 07:01  01:02  01:03  06:02  06:03
#> 6   NA11994 01:01 11:01 07:02 51:01 07:02 15:02  03:01  03:01  03:02  03:02
#>   DRB1.1 DRB1.2
#> 1  04:01  15:01
#> 2  15:01  15:01
#> 3  01:01  15:01
#> 4  01:01  03:01
#> 5  13:01  15:01
#> 6  04:02  04:04
dim(HLA_Type_Table)  # 60 13
#> [1] 60 13

# make a "hlaAlleleClass" object
hla.id <- "C"
hla <- hlaAllele(HLA_Type_Table$sample.id, HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  HLA_Type_Table[, paste(hla.id, ".2", sep="")], locus=hla.id, assembly="hg19")
summary(hla)
#> Gene: HLA - C
#> Range: [NAbp, NAbp]
#> # of samples: 60
#> # of unique KIR3DL1/S1 alleles: 17
#> # of unique KIR3DL1/S1 genotypes: 35

subhla1 <- hlaAlleleSubset(hla,   1:100)
summary(subhla1)
#> Gene: HLA - C
#> Range: [NAbp, NAbp]
#> # of samples: 100
#> # of unique KIR3DL1/S1 alleles: 17
#> # of unique KIR3DL1/S1 genotypes: 35
subhla2 <- hlaAlleleSubset(hla, 201:300)
summary(subhla2)
#> Gene: HLA - C
#> Range: [NAbp, NAbp]
#> # of samples: 100
#> # of unique KIR3DL1/S1 alleles: 0
#> # of unique KIR3DL1/S1 genotypes: 0

H <- hlaCombineAllele(subhla1, subhla2)
#> Error in hlaCombineAllele(subhla1, subhla2): H1$pos.start == H2$pos.start is not TRUE
summary(H)
#> Error: object 'H' not found
```
