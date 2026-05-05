# The class of HIBAG object

The class of a HIBAG model.

## Value

A list of:

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

- classifiers:

  a list of all classifiers (described as follows)

- appendix:

  `platform` – supported platform(s); `information` – other information,
  like training sets, authors; `warning` – any warning message

`classifiers` has the following components:

- samp.num:

  the number of copies of samples in a bootstrap sample

- haplos:

  a data.frame of haplotype frequencies

- :

  `freq` – haplotype frequency

- :

  `hla` – a HLA allele

- :

  `haplo` – a SNP haplotype, with an entry value 0 standing for B (ZERO
  A allele), 1 for A (ONE A allele)

- snpidx:

  the SNP indices used in this classifier

- outofbag.acc:

  the out-of-bag accuracy of this classifier

## Author

Xiuwen Zheng

## See also

[`hlaAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagging.md),
[`hlaParallelAttrBagging`](https://normanlabucd.github.io/PONG2/reference/hlaParallelAttrBagging.md),
[`hlaModelToObj`](https://normanlabucd.github.io/PONG2/reference/hlaModelToObj.md),
[`hlaModelFiles`](https://normanlabucd.github.io/PONG2/reference/hlaModelFiles.md),
[`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
