# HLA Locus Information

To get the starting and ending positions in basepair of HLA loci.

## Usage

``` r
hlaLociInfo(assembly=c("auto", "auto-silent", "hg18", "hg19", "unknown"))
```

## Arguments

- assembly:

  the human genome reference: "hg19" (default), "hg18", "auto" refers to
  "hg19"; "auto-silent" refers to "hg19" without any warning

## Value

Return a list:

- loci:

  the names of HLA classic I, II and III genes, such like A, B, C, DRB1,
  DQA1, DQB1 and DPB1

- pos.HLA.start:

  the starting position in basepair

- pos.HLA.end:

  the ending position in basepair

- length.HLA:

  the length of HLA genes in basepair

## Author

Xiuwen Zheng

## Examples

``` r
hlaLociInfo()
#> using the default genome assembly (assembly="hg38")
#>                name chrom    start      end suggest.pos
#> KIR2DL1     KIR2DL1    19 54769793 54784322          NA
#> KIR2DL2     KIR2DL2    19 54724442 55342622          NA
#> KIR2DL23   KIR2DL23    19 54724442 55342622          NA
#> KIR2DL3     KIR2DL3    19 54724442 54736632          NA
#> KIR2DL4     KIR2DL4    19 54803535 54814517          NA
#> KIR2DL5A   KIR2DL5A    19    60196    69656          -1
#> KIR2DL5B   KIR2DL5B    19   166421   175937          -1
#> KIR2DP1     KIR2DP1    19 54755023 54767371          NA
#> KIR2DS1     KIR2DS1    19 54724442 54867207          NA
#> KIR2DS2     KIR2DS2    19 54724442 54867207          NA
#> KIR2DS3     KIR2DS3    19   150177   164532          -1
#> KIR2DS4     KIR2DS4    19 54832676 54848569          NA
#> KIR2DS5     KIR2DS5    19 54724442 54867207          NA
#> KIR3DL1     KIR3DL1    19 54816438 54830778          NA
#> KIR3DL2     KIR3DL2    19 54850443 54867207          NA
#> KIR3DL3     KIR3DL3    19 54724447 54736536          NA
#> KIR3DP1     KIR3DP1    19 54786355 54790418          NA
#> KIR3DL1S1 KIR3DL1S1    19 54724442 54867207          NA
#> KIR3DS1     KIR3DS1    19 54724442 54867207          NA
```
