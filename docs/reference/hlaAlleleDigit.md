# Trim HLA alleles

Trim HLA alleles to specified width.

## Usage

``` r
hlaAlleleDigit(obj, max.resolution="4-digit", rm.suffix=FALSE)
```

## Arguments

- obj:

  should be a
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
  object or characters

- max.resolution:

  "2-digit", "4-digit", "6-digit", "8-digit", "allele", "protein", "2",
  "4", "6", "8", "full" or "": "allele" = "2-digit", "protein" =
  "4-digit", "full" and "" indicating no limit on resolution

- rm.suffix:

  whether remove the suffix, e.g., for "01:22N", "N" is a suffix

## Details

Either `HLAtypes` or `H1` `H2` should be specified. The format of
`HLAtypes` is "allele 1 / allele 2", e.g., "0512/0102". If
`max.resolution` is specified, the HLA alleles will be trimmed with the
maximum resolution.

## Value

Return a
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
object if `obj` is
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)-type,
or characters if `obj` is character-type.

## Author

Xiuwen Zheng

## See also

[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md)

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
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus = hla.id, assembly="hg19")
summary(hla)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 60
#> # of unique KIR3DL1/S1 alleles: 14
#> # of unique KIR3DL1/S1 genotypes: 29

hla2 <- hlaAlleleDigit(hla, "2-digit")
summary(hla2)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 60
#> # of unique KIR3DL1/S1 alleles: 12
#> # of unique KIR3DL1/S1 genotypes: 28
```
