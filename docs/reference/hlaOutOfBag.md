# Out-of-bag estimation of overall accuracy, per-allele sensitivity, etc

Out-of-bag estimation of overall accuracy, per-allele sensitivity,
specificity, positive predictive value, negative predictive value and
call rate.

## Usage

``` r
hlaOutOfBag(model, hla, snp, call.threshold=NaN, verbose=TRUE)
```

## Arguments

- model:

  an object of
  [`hlaAttrBagClass`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagClass.md)
  or
  [`hlaAttrBagObj`](https://normanlabucd.github.io/PONG2/reference/hlaAttrBagObj.md)

- hla:

  the training HLA types, an object of
  [`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md)

- snp:

  the training SNP genotypes, an object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- call.threshold:

  the specified call threshold; if `NaN`, no threshold is used

- verbose:

  if TRUE, show information

## Value

Return
[`hlaAlleleClass`](https://normanlabucd.github.io/PONG2/reference/hlaAlleleClass.md).

## Author

Xiuwen Zheng

## See also

[`hlaCompareAllele`](https://normanlabucd.github.io/PONG2/reference/hlaCompareAllele.md),
[`hlaReport`](https://normanlabucd.github.io/PONG2/reference/hlaReport.md)

## Examples

``` r
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="HIBAG")
data(HapMap_CEU_Geno, package="HIBAG")

# make a "hlaAlleleClass" object
hla.id <- "A"
hla <- hlaAllele(HLA_Type_Table$sample.id,
  H1 = HLA_Type_Table[, paste(hla.id, ".1", sep="")],
  H2 = HLA_Type_Table[, paste(hla.id, ".2", sep="")],
  locus=hla.id, assembly="hg19")

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"
length(snpid)  # 275
#> Error: object 'snpid' not found

# training and validation genotypes
geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel = match(hla$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'snpid' not found

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hla, geno, nclassifier=4)
#> 13565 monomorphic SNPs have been removed.
#> Build a PONG2 model with 4 individual classifiers:
#> # of SNPs randomly sampled as candidates for each selection: 38
#> # of SNPs: 1378, # of samples: 4
#> # of unique KIR3DL1/S1 alleles: 4
#> Tue May  5 19:51:22 2026,   1 individual classifier, out-of-bag acc: 100.00%, # of SNPs: 2, # of haplo: 4
#> Tue May  5 19:51:22 2026,   2 individual classifier, out-of-bag acc: 100.00%, # of SNPs: 2, # of haplo: 4
#> Tue May  5 19:51:22 2026,   3 individual classifier, out-of-bag acc: 50.00%, # of SNPs: 1, # of haplo: 2
#> Tue May  5 19:51:22 2026,   4 individual classifier, out-of-bag acc: 50.00%, # of SNPs: 1, # of haplo: 3
summary(model)
#> Gene: HLA - A
#> Training dataset: 4 samples X 1378 SNPs
#>  # of KIR3DL1/S1 alleles: 4
#>  # of individual classifiers: 4
#>  total # of SNPs used: 6
#>  average # of SNPs in an individual classifier: 1.50, sd: 0.58, min: 1, max: 2
#>  average # of haplotypes in an individual classifier: 3.25, sd: 0.96, min: 2, max: 4
#>  average out-of-bag accuracy: 75.00%, sd: 28.87%, min: 50.00%, max: 100.00%
#> Genome assembly: hg19

# out-of-bag estimation
(comp <- hlaOutOfBag(model, hla, geno, call.threshold=NaN, verbose=TRUE))
#> Gene: HLA - A
#> Training dataset: 4 samples X 1378 SNPs
#>  # of KIR3DL1/S1 alleles: 4
#>  # of individual classifiers: 4
#>  total # of SNPs used: 6
#>  average # of SNPs in an individual classifier: 1.50, sd: 0.58, min: 1, max: 2
#>  average # of haplotypes in an individual classifier: 3.25, sd: 0.96, min: 2, max: 4
#>  average out-of-bag accuracy: 75.00%, sd: 28.87%, min: 50.00%, max: 100.00%
#> Genome assembly: hg19
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"

# report
hlaReport(comp, type="txt")
#> Allele   Num.    Freq.   Num.    Freq.   CR  ACC SEN SPE PPV NPV Miscall
#>  Train   Train   Valid.  Valid.  (%) (%) (%) (%) (%) (%) (%)
#> ----
#> Overall accuracy: 84.8%, Call rate: 95.1%
#> KIR2DS5*0000 568 0.8304  372 0.8267  96.2    84.8    98.0    17.1    85.8    82.9    KIR2DS5*00201 (86)
#> KIR2DS5*00201    106 0.1550  70  0.1556  91.4    85.7    14.1    98.4    60.0    87.4    KIR2DS5*0000 (100)
#> KIR2DS5*00701    5   0.0073  3   0.0067  100.0   99.8    100.0   99.8    75.0    100.0   --
#> KIR2DS5*009  5   0.0073  5   0.0111  60.0    99.3    0.0 100.0   --  99.3    KIR2DS5*0000 (100)

hlaReport(comp, type="tex")
#> \title{Imputation Evaluation}
#> 
#> \documentclass[12pt]{article}
#> 
#> \usepackage{fullpage}
#> \usepackage{longtable}
#> 
#> \begin{document}
#> 
#> \maketitle
#> 
#> \setlength{\LTcapwidth}{6.4in}
#> 
#> % -------- BEGIN TABLE --------
#> \begin{longtable}{rrrrr | rrrrrrl}
#> \caption{The sensitivity (SEN), specificity (SPE), positive predictive value (PPV), negative predictive value (NPV) and call rate (CR)}
#> \label{tab:accuracy} \\
#> Allele & Num. & Freq. & Num. & Freq. & CR & ACC & SEN & SPE & PPV & NPV & Miscall \\
#>  & Train & Train & Valid. & Valid. & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) \\
#> \hline\hline
#> \endfirsthead
#> \multicolumn{12}{c}{{\normalsize \tablename\ \thetable{} -- Continued from previous page}} \\
#> Allele & Num. & Freq. & Num. & Freq. & CR & ACC & SEN & SPE & PPV & NPV & Miscall \\
#>  & Train & Train & Valid. & Valid. & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) & (\%) \\
#> \hline\hline
#> \endhead
#> \hline
#> \multicolumn{12}{r}{Continued on next page ...} \\
#> \hline
#> \endfoot
#> \hline\hline
#> \endlastfoot
#> \multicolumn{12}{l}{\it Overall accuracy: 84.8\%, Call rate: 95.1\%} \\
#> KIR2DS5*0000 & 568 & 0.8304 & 372 & 0.8267 & 96.2 & 84.8 & 98.0 & 17.1 & 85.8 & 82.9 & KIR2DS5*00201 (86) \\
#> KIR2DS5*00201 & 106 & 0.1550 & 70 & 0.1556 & 91.4 & 85.7 & 14.1 & 98.4 & 60.0 & 87.4 & KIR2DS5*0000 (100) \\
#> KIR2DS5*00701 & 5 & 0.0073 & 3 & 0.0067 & 100.0 & 99.8 & 100.0 & 99.8 & 75.0 & 100.0 & -- \\
#> KIR2DS5*009 & 5 & 0.0073 & 5 & 0.0111 & 60.0 & 99.3 & 0.0 & 100.0 & -- & 99.3 & KIR2DS5*0000 (100) \\
#> \end{longtable}
#> % -------- END TABLE --------
#> 
#> \end{document}

hlaReport(comp, type="html")
#> <!DOCTYPE html>
#> <html>
#> <head>
#>   <title>Imputation Evaluation</title>
#> </head>
#> <body>
#> <h1>Imputation Evaluation</h1>
#> <p></p>
#> <h3><b>Table 1:</b> The sensitivity (SEN), specificity (SPE), positive predictive value (PPV), negative predictive value (NPV) and call rate (CR).</h3>
#> <table id="TB-Acc" class="tabular" border="1"  CELLSPACING="1">
#> <tr>
#> <th>Allele </th> <th>Num. Train</th> <th>Freq. Train</th> <th>Num. Valid.</th> <th>Freq. Valid.</th> <th>CR (%)</th> <th>ACC (%)</th> <th>SEN (%)</th> <th>SPE (%)</th> <th>PPV (%)</th> <th>NPV (%)</th> <th>Miscall (%)</th>
#> </tr>
#> <tr>
#> <td colspan="12">
#> <i> Overall accuracy: 84.8%, Call rate: 95.1% </i>
#> </td>
#> </tr>
#> <tr>
#> <td>KIR2DS5*0000</td> <td>568</td> <td>0.8304</td> <td>372</td> <td>0.8267</td> <td>96.2</td> <td>84.8</td> <td>98.0</td> <td>17.1</td> <td>85.8</td> <td>82.9</td> <td>KIR2DS5*00201 (86)</td>
#> </tr>
#> <tr>
#> <td>KIR2DS5*00201</td> <td>106</td> <td>0.1550</td> <td>70</td> <td>0.1556</td> <td>91.4</td> <td>85.7</td> <td>14.1</td> <td>98.4</td> <td>60.0</td> <td>87.4</td> <td>KIR2DS5*0000 (100)</td>
#> </tr>
#> <tr>
#> <td>KIR2DS5*00701</td> <td>5</td> <td>0.0073</td> <td>3</td> <td>0.0067</td> <td>100.0</td> <td>99.8</td> <td>100.0</td> <td>99.8</td> <td>75.0</td> <td>100.0</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DS5*009</td> <td>5</td> <td>0.0073</td> <td>5</td> <td>0.0111</td> <td>60.0</td> <td>99.3</td> <td>0.0</td> <td>100.0</td> <td>--</td> <td>99.3</td> <td>KIR2DS5*0000 (100)</td>
#> </tr>
#> </table>
#> 
#> </body>
#> </html>
```
