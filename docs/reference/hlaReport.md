# Format a report

Create a report for evaluating prediction accuracies.

## Usage

``` r
hlaReport(object, export.fn="", type=c("txt", "tex", "html"),
  header=TRUE)
```

## Arguments

- object:

  an object returned by
  [`hlaCompareAllele`](https://normanlabucd.github.io/PONG2/reference/hlaCompareAllele.md)

- export.fn:

  a file name for output, or "" for `stdout`

- type:

  `"txt"` – tab-delimited text format; `"tex"` – tex format using the
  'longtable' package; `"html"` – html file

- header:

  if `TRUE`, output the header of text file associated corresponding
  format

## Value

None.

## Author

Xiuwen Zheng

## See also

[`hlaCompareAllele`](https://normanlabucd.github.io/PONG2/reference/hlaCompareAllele.md)

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

# divide HLA types randomly
set.seed(100)
hlatab <- hlaSplitAllele(hla, train.prop=0.5)
names(hlatab)
#> [1] "training"   "validation"
# "training"   "validation"
summary(hlatab$training)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 34
#> # of unique KIR3DL1/S1 alleles: 14
#> # of unique KIR3DL1/S1 genotypes: 23
summary(hlatab$validation)
#> Gene: HLA - A
#> Range: [NAbp, NAbp]
#> # of samples: 26
#> # of unique KIR3DL1/S1 alleles: 12
#> # of unique KIR3DL1/S1 genotypes: 14

# SNP predictors within the flanking region on each side
region <- 500   # kb
snpid <- hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,
  hla.id, region*1000, assembly="hg19")
#> Error in hlaFlankingSNP(HapMap_CEU_Geno$snp.id, HapMap_CEU_Geno$snp.position,     hla.id, region * 1000, assembly = "hg19"): could not find function "hlaFlankingSNP"
length(snpid)  # 275
#> Error: object 'snpid' not found

# training and validation genotypes
train.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  snp.sel = match(snpid, HapMap_CEU_Geno$snp.id),
  samp.sel = match(hlatab$training$value$sample.id, HapMap_CEU_Geno$sample.id))
#> Error: object 'snpid' not found
test.geno <- hlaGenoSubset(HapMap_CEU_Geno,
  samp.sel=match(hlatab$validation$value$sample.id, HapMap_CEU_Geno$sample.id))

# train a HIBAG model
set.seed(100)
# please use "nclassifier=100" when you use HIBAG for real data
model <- hlaAttrBagging(hlatab$training, train.geno, nclassifier=4,
  verbose.detail=TRUE)
#> Error: object 'train.geno' not found
summary(model)
#> Gene: KIR2DL23
#> Training dataset: 589 samples X 2770 SNPs
#>  # of KIR3DL1/S1 alleles: 23
#>  # of individual classifiers: 100
#>  total # of SNPs used: 2092
#>  average # of SNPs in an individual classifier: 103.39, sd: 14.99, min: 63, max: 128
#>  average # of haplotypes in an individual classifier: 421.30, sd: 62.46, min: 288, max: 614
#>  average out-of-bag accuracy: 95.74%, sd: 0.89%, min: 93.09%, max: 98.39%
#> Genome assembly: hg38

# validation
pred <- predict(model, test.geno)
#> Error in hlaPredict(object, snp, cl, type, vote, allele.check, match.type,     same.strand, verbose, verbose.match): could not find function "hlaPredict"
# compare
(comp <- hlaCompareAllele(hlatab$validation, pred, allele.limit=model,
  call.threshold=0))
#> Calling 'hlaCompareAllele': there are 0 individuals in common.
#> $overall
#>   total.num.ind crt.num.ind crt.num.haplo acc.ind acc.haplo call.threshold
#> 1             0           0             0     NaN       NaN              0
#>   n.call call.rate
#> 1      0       NaN
#> 
#> $confusion
#>                     True
#> Predict              KIR2DL2*00101 KIR2DL2*00301 KIR2DL2*005 KIR2DL2*00602
#>   KIR2DL2*00101                  0             0           0             0
#>   KIR2DL2*00301                  0             0           0             0
#>   KIR2DL2*005                    0             0           0             0
#>   KIR2DL2*00602                  0             0           0             0
#>   KIR2DL2*new                    0             0           0             0
#>   KIR2DL2*null                   0             0           0             0
#>   KIR2DL23*0000                  0             0           0             0
#>   KIR2DL3*00101                  0             0           0             0
#>   KIR2DL3*00102                  0             0           0             0
#>   KIR2DL3*00110                  0             0           0             0
#>   KIR2DL3*00201                  0             0           0             0
#>   KIR2DL3*003                    0             0           0             0
#>   KIR2DL3*00501                  0             0           0             0
#>   KIR2DL3*006                    0             0           0             0
#>   KIR2DL3*009                    0             0           0             0
#>   KIR2DL3*011                    0             0           0             0
#>   KIR2DL3*01202                  0             0           0             0
#>   KIR2DL3*013                    0             0           0             0
#>   KIR2DL3*015                    0             0           0             0
#>   KIR2DL3*019                    0             0           0             0
#>   KIR2DL3*022                    0             0           0             0
#>   KIR2DL3*null                   0             0           0             0
#>   KIR2DL3*unresolved             0             0           0             0
#>   ...                            0             0           0             0
#>                     True
#> Predict              KIR2DL2*new KIR2DL2*null KIR2DL23*0000 KIR2DL3*00101
#>   KIR2DL2*00101                0            0             0             0
#>   KIR2DL2*00301                0            0             0             0
#>   KIR2DL2*005                  0            0             0             0
#>   KIR2DL2*00602                0            0             0             0
#>   KIR2DL2*new                  0            0             0             0
#>   KIR2DL2*null                 0            0             0             0
#>   KIR2DL23*0000                0            0             0             0
#>   KIR2DL3*00101                0            0             0             0
#>   KIR2DL3*00102                0            0             0             0
#>   KIR2DL3*00110                0            0             0             0
#>   KIR2DL3*00201                0            0             0             0
#>   KIR2DL3*003                  0            0             0             0
#>   KIR2DL3*00501                0            0             0             0
#>   KIR2DL3*006                  0            0             0             0
#>   KIR2DL3*009                  0            0             0             0
#>   KIR2DL3*011                  0            0             0             0
#>   KIR2DL3*01202                0            0             0             0
#>   KIR2DL3*013                  0            0             0             0
#>   KIR2DL3*015                  0            0             0             0
#>   KIR2DL3*019                  0            0             0             0
#>   KIR2DL3*022                  0            0             0             0
#>   KIR2DL3*null                 0            0             0             0
#>   KIR2DL3*unresolved           0            0             0             0
#>   ...                          0            0             0             0
#>                     True
#> Predict              KIR2DL3*00102 KIR2DL3*00110 KIR2DL3*00201 KIR2DL3*003
#>   KIR2DL2*00101                  0             0             0           0
#>   KIR2DL2*00301                  0             0             0           0
#>   KIR2DL2*005                    0             0             0           0
#>   KIR2DL2*00602                  0             0             0           0
#>   KIR2DL2*new                    0             0             0           0
#>   KIR2DL2*null                   0             0             0           0
#>   KIR2DL23*0000                  0             0             0           0
#>   KIR2DL3*00101                  0             0             0           0
#>   KIR2DL3*00102                  0             0             0           0
#>   KIR2DL3*00110                  0             0             0           0
#>   KIR2DL3*00201                  0             0             0           0
#>   KIR2DL3*003                    0             0             0           0
#>   KIR2DL3*00501                  0             0             0           0
#>   KIR2DL3*006                    0             0             0           0
#>   KIR2DL3*009                    0             0             0           0
#>   KIR2DL3*011                    0             0             0           0
#>   KIR2DL3*01202                  0             0             0           0
#>   KIR2DL3*013                    0             0             0           0
#>   KIR2DL3*015                    0             0             0           0
#>   KIR2DL3*019                    0             0             0           0
#>   KIR2DL3*022                    0             0             0           0
#>   KIR2DL3*null                   0             0             0           0
#>   KIR2DL3*unresolved             0             0             0           0
#>   ...                            0             0             0           0
#>                     True
#> Predict              KIR2DL3*00501 KIR2DL3*006 KIR2DL3*009 KIR2DL3*011
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*01202 KIR2DL3*013 KIR2DL3*015 KIR2DL3*019
#>   KIR2DL2*00101                  0           0           0           0
#>   KIR2DL2*00301                  0           0           0           0
#>   KIR2DL2*005                    0           0           0           0
#>   KIR2DL2*00602                  0           0           0           0
#>   KIR2DL2*new                    0           0           0           0
#>   KIR2DL2*null                   0           0           0           0
#>   KIR2DL23*0000                  0           0           0           0
#>   KIR2DL3*00101                  0           0           0           0
#>   KIR2DL3*00102                  0           0           0           0
#>   KIR2DL3*00110                  0           0           0           0
#>   KIR2DL3*00201                  0           0           0           0
#>   KIR2DL3*003                    0           0           0           0
#>   KIR2DL3*00501                  0           0           0           0
#>   KIR2DL3*006                    0           0           0           0
#>   KIR2DL3*009                    0           0           0           0
#>   KIR2DL3*011                    0           0           0           0
#>   KIR2DL3*01202                  0           0           0           0
#>   KIR2DL3*013                    0           0           0           0
#>   KIR2DL3*015                    0           0           0           0
#>   KIR2DL3*019                    0           0           0           0
#>   KIR2DL3*022                    0           0           0           0
#>   KIR2DL3*null                   0           0           0           0
#>   KIR2DL3*unresolved             0           0           0           0
#>   ...                            0           0           0           0
#>                     True
#> Predict              KIR2DL3*022 KIR2DL3*null KIR2DL3*unresolved
#>   KIR2DL2*00101                0            0                  0
#>   KIR2DL2*00301                0            0                  0
#>   KIR2DL2*005                  0            0                  0
#>   KIR2DL2*00602                0            0                  0
#>   KIR2DL2*new                  0            0                  0
#>   KIR2DL2*null                 0            0                  0
#>   KIR2DL23*0000                0            0                  0
#>   KIR2DL3*00101                0            0                  0
#>   KIR2DL3*00102                0            0                  0
#>   KIR2DL3*00110                0            0                  0
#>   KIR2DL3*00201                0            0                  0
#>   KIR2DL3*003                  0            0                  0
#>   KIR2DL3*00501                0            0                  0
#>   KIR2DL3*006                  0            0                  0
#>   KIR2DL3*009                  0            0                  0
#>   KIR2DL3*011                  0            0                  0
#>   KIR2DL3*01202                0            0                  0
#>   KIR2DL3*013                  0            0                  0
#>   KIR2DL3*015                  0            0                  0
#>   KIR2DL3*019                  0            0                  0
#>   KIR2DL3*022                  0            0                  0
#>   KIR2DL3*null                 0            0                  0
#>   KIR2DL3*unresolved           0            0                  0
#>   ...                          0            0                  0
#> 
#> $detail
#>                allele train.num   train.freq valid.num valid.freq call.rate
#> 1       KIR2DL2*00101       176 0.1494057725         0        NaN         0
#> 2       KIR2DL2*00301       103 0.0874363328         0        NaN         0
#> 3         KIR2DL2*005         1 0.0008488964         0        NaN         0
#> 4       KIR2DL2*00602        16 0.0135823430         0        NaN         0
#> 5         KIR2DL2*new         1 0.0008488964         0        NaN         0
#> 6        KIR2DL2*null         2 0.0016977929         0        NaN         0
#> 7       KIR2DL23*0000        33 0.0280135823         0        NaN         0
#> 8       KIR2DL3*00101       560 0.4753820034         0        NaN         0
#> 9       KIR2DL3*00102         4 0.0033955857         0        NaN         0
#> 10      KIR2DL3*00110         4 0.0033955857         0        NaN         0
#> 11      KIR2DL3*00201       169 0.1434634975         0        NaN         0
#> 12        KIR2DL3*003         1 0.0008488964         0        NaN         0
#> 13      KIR2DL3*00501        55 0.0466893039         0        NaN         0
#> 14        KIR2DL3*006        12 0.0101867572         0        NaN         0
#> 15        KIR2DL3*009         7 0.0059422750         0        NaN         0
#> 16        KIR2DL3*011         3 0.0025466893         0        NaN         0
#> 17      KIR2DL3*01202         9 0.0076400679         0        NaN         0
#> 18        KIR2DL3*013         5 0.0042444822         0        NaN         0
#> 19        KIR2DL3*015         1 0.0008488964         0        NaN         0
#> 20        KIR2DL3*019         3 0.0025466893         0        NaN         0
#> 21        KIR2DL3*022         1 0.0008488964         0        NaN         0
#> 22       KIR2DL3*null         9 0.0076400679         0        NaN         0
#> 23 KIR2DL3*unresolved         3 0.0025466893         0        NaN         0
#>    accuracy sensitivity specificity ppv npv miscall miscall.prop
#> 1       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 2       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 3       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 4       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 5       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 6       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 7       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 8       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 9       NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 10      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 11      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 12      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 13      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 14      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 15      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 16      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 17      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 18      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 19      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 20      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 21      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 22      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 23      NaN         NaN         NaN NaN NaN    <NA>          NaN
#> 


# report
hlaReport(comp, type="txt")
#> Allele   Num.    Freq.   Num.    Freq.   CR  ACC SEN SPE PPV NPV Miscall
#>  Train   Train   Valid.  Valid.  (%) (%) (%) (%) (%) (%) (%)
#> ----
#> Overall accuracy: NaN%, Call rate: NaN%
#> KIR2DL2*00101    176 0.1494  0   0   --  --  --  --  --  --  --
#> KIR2DL2*00301    103 0.0874  0   0   --  --  --  --  --  --  --
#> KIR2DL2*005  1   0.0008  0   0   --  --  --  --  --  --  --
#> KIR2DL2*00602    16  0.0136  0   0   --  --  --  --  --  --  --
#> KIR2DL2*new  1   0.0008  0   0   --  --  --  --  --  --  --
#> KIR2DL2*null 2   0.0017  0   0   --  --  --  --  --  --  --
#> KIR2DL23*0000    33  0.0280  0   0   --  --  --  --  --  --  --
#> KIR2DL3*00101    560 0.4754  0   0   --  --  --  --  --  --  --
#> KIR2DL3*00102    4   0.0034  0   0   --  --  --  --  --  --  --
#> KIR2DL3*00110    4   0.0034  0   0   --  --  --  --  --  --  --
#> KIR2DL3*00201    169 0.1435  0   0   --  --  --  --  --  --  --
#> KIR2DL3*003  1   0.0008  0   0   --  --  --  --  --  --  --
#> KIR2DL3*00501    55  0.0467  0   0   --  --  --  --  --  --  --
#> KIR2DL3*006  12  0.0102  0   0   --  --  --  --  --  --  --
#> KIR2DL3*009  7   0.0059  0   0   --  --  --  --  --  --  --
#> KIR2DL3*011  3   0.0025  0   0   --  --  --  --  --  --  --
#> KIR2DL3*01202    9   0.0076  0   0   --  --  --  --  --  --  --
#> KIR2DL3*013  5   0.0042  0   0   --  --  --  --  --  --  --
#> KIR2DL3*015  1   0.0008  0   0   --  --  --  --  --  --  --
#> KIR2DL3*019  3   0.0025  0   0   --  --  --  --  --  --  --
#> KIR2DL3*022  1   0.0008  0   0   --  --  --  --  --  --  --
#> KIR2DL3*null 9   0.0076  0   0   --  --  --  --  --  --  --
#> KIR2DL3*unresolved   3   0.0025  0   0   --  --  --  --  --  --  --

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
#> \multicolumn{12}{l}{\it Overall accuracy: NaN\%, Call rate: NaN\%} \\
#> KIR2DL2*00101 & 176 & 0.1494 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL2*00301 & 103 & 0.0874 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL2*005 & 1 & 0.0008 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL2*00602 & 16 & 0.0136 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL2*new & 1 & 0.0008 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL2*null & 2 & 0.0017 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL23*0000 & 33 & 0.0280 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*00101 & 560 & 0.4754 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*00102 & 4 & 0.0034 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*00110 & 4 & 0.0034 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*00201 & 169 & 0.1435 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*003 & 1 & 0.0008 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*00501 & 55 & 0.0467 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*006 & 12 & 0.0102 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*009 & 7 & 0.0059 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*011 & 3 & 0.0025 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*01202 & 9 & 0.0076 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*013 & 5 & 0.0042 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*015 & 1 & 0.0008 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*019 & 3 & 0.0025 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*022 & 1 & 0.0008 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*null & 9 & 0.0076 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
#> KIR2DL3*unresolved & 3 & 0.0025 & 0 & 0 & -- & -- & -- & -- & -- & -- & -- \\
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
#> <i> Overall accuracy: NaN%, Call rate: NaN% </i>
#> </td>
#> </tr>
#> <tr>
#> <td>KIR2DL2*00101</td> <td>176</td> <td>0.1494</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL2*00301</td> <td>103</td> <td>0.0874</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL2*005</td> <td>1</td> <td>0.0008</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL2*00602</td> <td>16</td> <td>0.0136</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL2*new</td> <td>1</td> <td>0.0008</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL2*null</td> <td>2</td> <td>0.0017</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL23*0000</td> <td>33</td> <td>0.0280</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*00101</td> <td>560</td> <td>0.4754</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*00102</td> <td>4</td> <td>0.0034</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*00110</td> <td>4</td> <td>0.0034</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*00201</td> <td>169</td> <td>0.1435</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*003</td> <td>1</td> <td>0.0008</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*00501</td> <td>55</td> <td>0.0467</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*006</td> <td>12</td> <td>0.0102</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*009</td> <td>7</td> <td>0.0059</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*011</td> <td>3</td> <td>0.0025</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*01202</td> <td>9</td> <td>0.0076</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*013</td> <td>5</td> <td>0.0042</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*015</td> <td>1</td> <td>0.0008</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*019</td> <td>3</td> <td>0.0025</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*022</td> <td>1</td> <td>0.0008</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*null</td> <td>9</td> <td>0.0076</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> <tr>
#> <td>KIR2DL3*unresolved</td> <td>3</td> <td>0.0025</td> <td>0</td> <td>0</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td> <td>--</td>
#> </tr>
#> </table>
#> 
#> </body>
#> </html>
```
