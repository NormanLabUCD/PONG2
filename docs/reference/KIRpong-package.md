# HLA Genotype Imputation with Attribute Bagging

To impute HLA types from unphased SNP data using an attribute bagging
method.

## Details

|          |               |
|----------|---------------|
| Package: | PONG          |
| Type:    | Package       |
| Version: | 1.2.4         |
| License: | GPL version 3 |

PONG is a state of the art software package for imputing HLA types using
SNP data, and it uses the R statistical programming language. PONG is
highly accurate, computationally tractable, and can be used by
researchers with published parameter estimates instead of requiring
access to large training sample datasets. It combines the concepts of
attribute bagging, an ensemble classifier method, with haplotype
inference for SNPs and HLA types. Attribute bagging is a technique which
improves the accuracy and stability of classifier ensembles using
bootstrap aggregating and random variable selection.

**Features:**  
1) PONG can be used by researchers with published parameter estimates
(<http://www.biostat.washington.edu/~bsweir/PONG/>) instead of requiring
access to large training sample datasets.  
2) A typical PONG parameter file contains only haplotype frequencies at
different SNP subsets rather than individual training genotypes.  
3) SNPs within the xMHC region (chromosome 6) are used for imputation.  
4) PONG employs unphased genotypes of unrelated individuals as a
training set.  
5) PONG supports parallel computing with R.  

## References

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS; HIBAG –
HLA Genotype Imputation with Attribute Bagging; (Abstract 294,
Platform/Oral Talk); Present at the 62nd Annual Meeting of the American
Society of Human Genetics, November 9, 2012 in San Francisco,
California.

Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS; HIBAG –
HLA Genotype Imputation with Attribute Bagging. Pharmacogenomics
Journal. doi: 10.1038/tpj.2013.18.
<http://dx.doi.org/10.1038/tpj.2013.18>

## Author

Xiuwen Zheng <zhengx@u.washington.edu>, Bruce S. Weir
<bsweir@u.washington.edu>

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="PONG")
#> Error in find.package(package, lib.loc, verbose = verbose): there is no package called ‘PONG’
data(HapMap_CEU_Geno, package="PONG")
#> Error in find.package(package, lib.loc, verbose = verbose): there is no package called ‘PONG’

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

summary(HapMap_CEU_Geno)
#> SNP genotypes: 
#>  60 samples X 1564 SNPs
#>  SNPs range from 25769023bp to 33421576bp on hg19
#> Missing rate per SNP:
#>  min: 0, max: 0.0666667, mean: 0.0651215, median: 0.0666667, sd: 0.00968287
#> Missing rate per sample:
#>  min: 0, max: 0.969949, mean: 0.0651215, median: 0.000639386, sd: 0.243737
#> Minor allele frequency:
#>  min: 0, max: 0.5, mean: 0.227582, median: 0.205357, sd: 0.1389
#> Allelic information:
#> A/C A/G C/T G/T 
#> 136 655 632 141 


######################################################################
# Temporal Fixed
HLA_Type_Table$KIR3DLS1.1 <-HLA_Type_Table$A.1
HLA_Type_Table$KIR3DLS1.2 <-HLA_Type_Table$A.2
# make a "hlaAlleleClass" object
hla.id <- "KIR3DLS1"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
#> [1] "training"   "validation"
# "training"   "validation"
summary(hlatab$training)
#> Gene: KIR3DLS1
#> Range: [NAbp, NAbp]
#> # of samples: 34
#> # of unique KIR3DL1/S1 alleles: 14
#> # of unique KIR3DL1/S1 genotypes: 23
summary(hlatab$validation)
#> Gene: KIR3DLS1
#> Range: [NAbp, NAbp]
#> # of samples: 26
#> # of unique KIR3DL1/S1 alleles: 12
#> # of unique KIR3DL1/S1 genotypes: 14


```
