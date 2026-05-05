# The class of HIBAG model

The class of a HIBAG model, and its instance is returned from
[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md).

## Value

Return a list of:

- n.samp:

  the total number of training samples

- n.snp:

  the total number of candidate SNP predictors

- sample.id:

  the sample IDs

- snp.id:

  the SNP IDs

- snp.position:

  SNP position in basepair

- snp.allele:

  a vector of characters with the format of “A allele/B allele”

- snp.allele.freq:

  the allele frequencies

- hla.locus:

  the name of HLA locus

- hla.allele:

  the HLA alleles used in the model

- hla.freq:

  the HLA allele frequencies

- assembly:

  the human genome reference, such like "hg19"

- model:

  internal use

- appendix:

  an optional list: `platform` – supported platform(s); `information` –
  other information, like training sets, authors; `warning` – any
  warning message

## References

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS; HIBAG –
HLA Genotype Imputation with Attribute Bagging; (Abstract 294,
Platform/Oral Talk); Present at the 62nd Annual Meeting of the American
Society of Human Genetics, November 9, 2012 in San Francisco,
California.

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS; HIBAG –
HLA Genotype Imputation with Attribute Bagging. To appear in the
Pharmacogenomics Journal.

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaParallelAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaParallelAttrBagging.md),
[`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)
