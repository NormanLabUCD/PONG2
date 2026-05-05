# Get sample IDs from HLA types with a filter

Get sample IDs from HLA types limited to a set of HLA alleles.

## Usage

``` r
hlaSampleAllele(TrueHLA, allele.limit=NULL, max.resolution="")
```

## Arguments

- TrueHLA:

  an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)

- allele.limit:

  a list of HLA alleles, the validation samples are limited to those
  having HLA alleles in `allele.limit`, or `NULL` for no limit.
  `allele.limit` could be character-type,
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
  or
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- max.resolution:

  "2-digit", "4-digit", "6-digit", "8-digit", "allele", "protein", "2",
  "4", "6", "8", "full" or "": "allele" = "2-digit", "protein" =
  "4-digit", "full" and "" mean no limit on resolution

## Value

Return a list of sample IDs.

## Author

Xiuwen Zheng

## See also

[`hlaCompareAllele`](https://normanlabucd.github.io/PONG2/reference/hlaCompareAllele.md)
