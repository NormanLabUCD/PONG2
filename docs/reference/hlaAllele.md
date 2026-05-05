# A list of HLA types

Return an object of
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md),
which contains HLA types.

## Usage

``` r
hlaAllele(sample.id, H1, H2, max.resolution="", locus="any",
  assembly=c("auto", "auto-silent", "hg18", "hg19", "unknown"),
  locus.pos.start=NA, locus.pos.end=NA, prob=NULL, na.rm=TRUE)
```

## Arguments

- sample.id:

  sample IDs

- H1:

  a vector of HLA alleles

- H2:

  a vector of HLA alleles

- max.resolution:

  "2-digit", "4-digit", "6-digit", "8-digit", "allele", "protein", "2",
  "4", "6", "8", "full" or "": "allele" = "2-digit", "protein" =
  "4-digit", "full" and "" indicating no limit on resolution

- locus:

  the name of HLA locus: "A", "B", "C", "DRB1", "DRB5", "DQA1", "DQB1",
  "DPB1", or "any", where "any" indicates any other multiallelic locus

- assembly:

  the human genome reference: "hg19" (default), "hg18", "auto" refers to
  "hg19"; "auto-silent" refers to "hg19" without any warning

- locus.pos.start:

  the starting position in basepair

- locus.pos.end:

  the end position in basepair

- prob:

  the probabilities assigned to the samples

- na.rm:

  if TRUE, remove the samples without valid HLA types

## Details

The format of `H1` and `H2` is "allele group : different protein :
synonymous mutations in exons : synonymous mutations in introns"L, where
the suffix L is express level (N, null; L, low; S, secreted; A,
aberrant; Q: questionable). For example, "44:02:01:02L". If
`max.resolution` is specified, the HLA alleles will be trimmed with a
possible maximum resolution.

## Value

Return a
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
object, and it is a list:

- locus:

  HLA locus

- pos.start:

  the starting position in basepair

- pos.end:

  the end position in basepair

- value:

  a data frame

- assembly:

  the human genome reference, such like "hg19"

The component `value` includes:

- sample.id:

  sample ID

- allele1:

  HLA allele

- allele2:

  HLA allele

- prob:

  the posterior probability

## Author

Xiuwen Zheng

## See also

[`hlaAlleleDigit`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleDigit.md),
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
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")
summary(hla)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 60
#> # of unique KIR3DL1/S1 alleles: 14
#> # of unique KIR3DL1/S1 genotypes: 29
```
