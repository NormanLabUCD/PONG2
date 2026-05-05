# Get SNP IDs and positions

Get the information of SNP ID with or without position.

## Usage

``` r
hlaSNPID(obj, type=c("RefSNP+Position", "RefSNP", "Position"))
```

## Arguments

- obj:

  a genotypic object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md),
  a haplotypic object of
  [`hlaSNPHaploClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPHaploClass.md),
  a model object of
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
  or a model object of
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- type:

  `"RefSNP+Position"` (by default), `"RefSNP"` or `"Position"`

## Value

If `type = "RefSNP+Position"`, return
`paste(obj$snp.id, obj$snp.position, sep="-")`; if `type = "RefSNP"`,
return `obj$snp.id`; otherwise, return `obj$snp.position`.

## Author

Xiuwen Zheng

## See also

[`hlaGenoSwitchStrand`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSwitchStrand.md),
[`hlaGenoCombine`](https://normanlabucd.github.io/PONG2/reference/hlaGenoCombine.md)

## Examples

``` r
# load SNP genotypes
data(HapMap_CEU_Geno, package="HIBAG")

x <- hlaSNPID(HapMap_CEU_Geno)
head(x)
#> [1] 25769023 25770238 25776948 25779735 25783314 25800196

x <- hlaSNPID(HapMap_CEU_Geno, "RefSNP")
head(x)
#> [1] "rs1892250"  "rs2328893"  "rs11754288" "rs3734526"  "rs3923"    
#> [6] "rs3823151" 

x <- hlaSNPID(HapMap_CEU_Geno, "Position")
head(x)
#> [1] 25769023 25770238 25776948 25779735 25783314 25800196
```
