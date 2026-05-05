# Get unique HLA alleles

Get unique HLA alleles, which are in ascending order.

## Usage

``` r
hlaUniqueAllele(hla)
```

## Arguments

- hla:

  character-type HLA alleles, or a
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)
  object

## Details

Each HLA allele name has a unique number corresponding to up to four
sets of digits separated by colons. The name designation depends on the
sequence of the allele and that of its nearest relative. The digits
before the first colon describe the type, which often corresponds to the
serological antigen carried by an allotype. The next set of digits are
used to list the subtypes, numbers being assigned in the order in which
DNA sequences have been determined. Alleles whose numbers differ in the
two sets of digits must differ in one or more nucleotide substitutions
that change the amino acid sequence of the encoded protein. Alleles that
differ only by synonymous nucleotide substitutions (also called silent
or non-coding substitutions) within the coding sequence are
distinguished by the use of the third set of digits. Alleles that only
differ by sequence polymorphisms in the introns or in the 5' or 3'
untranslated regions that flank the exons and introns are distinguished
by the use of the fourth set of digits.

In addition to the unique allele number there are additional optional
suffixes that may be added to an allele to indicate its expression
status. Alleles that have been shown not to be expressed, 'Null' alleles
have been given the suffix 'N'. Those alleles which have been shown to
be alternatively expressed may have the suffix 'L', 'S', 'C', 'A' or
'Q'.

[http://hla.alleles.org/nomenclature/index.html](http://hla.alleles.org/nomenclature/index.md)

## Value

Return a vector of HLA alleles

## Author

Xiuwen Zheng

## See also

[`hlaAllele`](https://normanlabucd.github.io/PONG2/reference/hlaAllele.md),
[`hlaAlleleDigit`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleDigit.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")

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
hlaUniqueAllele(hla)
#>  [1] "01:01" "02:01" "02:06" "03:01" "11:01" "23:01" "24:02" "24:03" "25:01"
#> [10] "26:01" "29:02" "31:01" "32:01" "68:01"

hlaUniqueAllele(c("01", "01:03", "01:01", "03:05", "03:01g", "03:01", "104:01"))
#> [1] "01"     "01:01"  "01:03"  "03:01"  "03:01g" "03:05"  "104:01"
```
