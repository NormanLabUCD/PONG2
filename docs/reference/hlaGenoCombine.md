# Combine two genotypic data sets into one

To combine two genotypic data sets into one dataset.

## Usage

``` r
hlaGenoCombine(geno1, geno2,
  match.type=c("RefSNP+Position", "RefSNP", "Position"),
  allele.check=TRUE, same.strand=FALSE, verbose=TRUE)
```

## Arguments

- geno1:

  the first genotype object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- geno2:

  the second genotype object of
  [`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md)

- match.type:

  `"RefSNP+Position"` (by default) – using both of RefSNP IDs and
  positions; `"RefSNP"` – using RefSNP IDs only; `"Position"` – using
  positions only

- allele.check:

  if `TRUE`, call
  [`hlaGenoSwitchStrand`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSwitchStrand.md)
  to check and then switch allele pairs if needed

- same.strand:

  `TRUE` assuming alleles are on the same strand (e.g., forward strand);
  otherwise, `FALSE` not assuming whether on the same strand or not

- verbose:

  show information, if TRUE

## Details

The function merges two SNP dataset `geno1` and `geno2`, and returns a
SNP dataset consisting of the SNP intersect between `geno1` and `geno2`,
and having the same SNP information (allele and position) as `geno1`.

## Value

An object of
[`hlaSNPGenoClass`](https://normanlabucd.github.io/PONG2/reference/hlaSNPGenoClass.md).

## Author

Xiuwen Zheng

## See also

[`hlaMakeSNPGeno`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPGeno.md),
[`hlaMakeSNPHaplo`](https://normanlabucd.github.io/PONG2/reference/hlaMakeSNPHaplo.md),
[`hlaHaploSubset`](https://normanlabucd.github.io/PONG2/reference/hlaHaploSubset.md),
[`hlaGenoSubset`](https://normanlabucd.github.io/PONG2/reference/hlaGenoSubset.md)

## Examples

``` r
# load SNP genotypes
data(HapMap_CEU_Geno, package="HIBAG")

# import a PLINK BED file
bed.fn <- system.file("extdata", "HapMap_CEU.bed", package="HIBAG")
fam.fn <- system.file("extdata", "HapMap_CEU.fam", package="HIBAG")
bim.fn <- system.file("extdata", "HapMap_CEU.bim", package="HIBAG")
hapmap.ceu <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly="hg19")
#> Open "/usr/local/lib/R/site-library/HIBAG/extdata/HapMap_CEU.bed" in the individual-major mode.
#> Open "/usr/local/lib/R/site-library/HIBAG/extdata/HapMap_CEU.fam".
#> Open "/usr/local/lib/R/site-library/HIBAG/extdata/HapMap_CEU.bim".
#> Import 0 SNPs within the KIR gene cluster on chromosome 19.
#> Error in hlaBED2Geno(bed.fn, fam.fn, bim.fn, assembly = "hg19"): There is no SNP imported.

# combine two datasets together
geno <- hlaGenoCombine(HapMap_CEU_Geno, hapmap.ceu)
#> Error: object 'hapmap.ceu' not found
summary(geno)
#> SNP genotypes: 
#>  100 samples X 14943 SNPs
#>  SNPs range from 55000044bp to 55399932bp on hg19
#> Missing rate per SNP:
#>  min: 0, max: 0, mean: 0, median: 0, sd: 0
#> Missing rate per sample:
#>  min: 0, max: 0, mean: 0, median: 0, sd: 0
#> Minor allele frequency:
#>  min: 0, max: 0.5, mean: 0.0302637, median: 0, sd: 0.0826476
#> Allelic information:
#>                                      <CN0>/A 
#>                                            2 
#>                                      <CN0>/C 
#>                                            2 
#>                                      <CN0>/G 
#>                                            2 
#>                                      <CN2>/C 
#>                                            1 
#>                               <INS:ME:ALU>/C 
#>                                            1 
#>                             <INS:ME:LINE1>/C 
#>                                            1 
#>                                A/AAAAACAAAAC 
#>                                            1 
#>                      A/AAAAGAAAGAAAGAAAGAAAG 
#>                                            1 
#>                                       A/AAAC 
#>                                            1 
#>                                        A/AAC 
#>                                            1 
#>                                        A/AAG 
#>                                            2 
#>                                      A/AAGAG 
#>                                            1 
#>                                      A/AAGAT 
#>                                            2 
#>                                       A/AAGG 
#>                                            1 
#>                                      A/AATAG 
#>                                            1 
#>                                         A/AC 
#>                                            8 
#>                                       A/ACCT 
#>                                            1 
#>                            A/ACCTCAGATGGAGAG 
#>                                            1 
#>                                        A/ACT 
#>                                            1 
#>                                      A/ACTGT 
#>                                            1 
#>                                       A/ACTT 
#>                                            1 
#>                                         A/AG 
#>                                            6 
#>                                    A/AGCACAG 
#>                                            1 
#>                                       A/AGGG 
#>                                            1 
#>                                  A/AGTGTGTTG 
#>                                            1 
#>                                         A/AT 
#>                                           21 
#>                                     A/ATAAAT 
#>                                            1 
#>                                      A/ATAAG 
#>                                            1 
#>                                      A/ATACT 
#>                                            1 
#>                                       A/ATAT 
#>                                            1 
#>                                     A/ATATAT 
#>                                            1 
#>                                        A/ATC 
#>                                            1 
#>                                       A/ATCT 
#>                                            1 
#>                                        A/ATG 
#>                                            1 
#>                                      A/ATTAG 
#>                                            1 
#>                                       A/ATTC 
#>                                            1 
#>                                     A/ATTTGT 
#>                                            1 
#>                                          A/C 
#>                                          741 
#>                                          A/G 
#>                                         3331 
#>                                          A/T 
#>                                          475 
#>                                   AAAAAAAG/A 
#>                                            1 
#>                                      AAAAC/A 
#>                                            1 
#>                                      AAAAT/A 
#>                                            1 
#>                  AAATAAATAAATAAATAAATAAATT/A 
#>                                            1 
#>                                      AAATT/A 
#>                                            1 
#>                                         AC/A 
#>                                            6 
#>                                      ACAGG/A 
#>                                            1 
#>                                 ACATTTTCTT/A 
#>                                            1 
#>                                       ACCC/A 
#>                                            1 
#>                                        ACT/A 
#>                                            2 
#>                                         AG/A 
#>                                            4 
#>                                       AGAT/A 
#>                                            1 
#>                        AGCTTCAAGGCCTTAATGC/A 
#>                                            1 
#>                                  AGGACAAAG/A 
#>                                            1 
#>                                         AT/A 
#>                                           14 
#>                                        ATC/A 
#>                                            1 
#>                       ATGGGCCTGAAGTGGAGATC/A 
#>                                            1 
#>                             ATTCTTTTTTTTTT/A 
#>                                            1 
#>                                      ATTTT/A 
#>                                            1 
#>                                          C/A 
#>                                          492 
#>                                         C/CA 
#>                                            8 
#>                                        C/CAA 
#>                                            1 
#>                                    C/CAAAAAA 
#>                                            1 
#>                                      C/CAAAT 
#>                                            1 
#>                       C/CACATTTAAAATAAACAATA 
#>                                            1 
#>                                      C/CACTT 
#>                                            1 
#>                                        C/CAG 
#>                                            1 
#>                                      C/CAGAG 
#>                                            1 
#>                                        C/CAT 
#>                                            3 
#>                                       C/CATG 
#>                                            1 
#>                                     C/CATTAG 
#>                                            1 
#>                                      C/CATTT 
#>                                            1 
#>                                      C/CCACT 
#>                                            1 
#>                                    C/CCATCTT 
#>                                            1 
#>                                       C/CCCA 
#>                                            1 
#>                                      C/CCTAT 
#>                                            1 
#>                                       C/CCTT 
#>                                            3 
#>                                         C/CG 
#>                                            2 
#>                                        C/CGA 
#>                                            1 
#>                                        C/CGT 
#>                                            1 
#>                                         C/CT 
#>                                           17 
#>                                        C/CTA 
#>                                            1 
#>                                       C/CTAA 
#>                                            1 
#>                                      C/CTCTT 
#>                                            1 
#>                                        C/CTG 
#>                                            4 
#>                                       C/CTGA 
#>                                            1 
#>                                    C/CTGTGTG 
#>                                            1 
#>                                        C/CTT 
#>                                            1 
#>                                       C/CTTA 
#>                                            2 
#>                                      C/CTTAT 
#>                                            1 
#>                                       C/CTTT 
#>                                            1 
#>                                      C/CTTTT 
#>                                            1 
#>                                     C/CTTTTT 
#>                                            1 
#>                                          C/G 
#>                                          811 
#>                                          C/T 
#>                                         1650 
#>                                         CA/C 
#>                                           15 
#>                                        CAA/C 
#>                                            1 
#>                                     CAAAAA/C 
#>                                            1 
#>                                     CAAACA/C 
#>                                            1 
#>                                      CAATT/C 
#>                                            1 
#>                                CACATGTGTGT/C 
#>                                            1 
#>                                      CACCA/C 
#>                                            1 
#>                                        CAG/C 
#>                                            1 
#>                                        CAT/C 
#>                                            1 
#>                                       CATG/C 
#>                                            1 
#>                       CCTGCCCTGCGGTGTCCATG/C 
#>                                            1 
#>                                         CG/C 
#>                                            1 
#>                                    CGGTGTG/C 
#>                                            1 
#>                                        CGT/C 
#>                                            2 
#>                                       CGTA/C 
#>                                            1 
#>                                         CT/C 
#>                                            9 
#>                                        CTA/C 
#>                                            1 
#>                                      CTCAT/C 
#>                                            1 
#>                               CTCCAGGCAAAT/C 
#>                                            1 
#>                                        CTT/C 
#>                                            2 
#>                                          G/A 
#>                                         1638 
#>                                          G/C 
#>                                          717 
#>                                         G/GA 
#>                                            7 
#>                                        G/GAA 
#>                                            1 
#>                                      G/GAAGA 
#>                                            1 
#>                                      G/GAGAC 
#>                                            1 
#>                                        G/GAT 
#>                                            1 
#>                                       G/GATA 
#>                                            1 
#>                                         G/GC 
#>                                            8 
#>                                       G/GCCC 
#>                                            1 
#>                                       G/GCTC 
#>                                            1 
#>                                        G/GGA 
#>                                            1 
#> G/GGACCTCAGGCTCCTATGGTCTCCCCCTGTATGTTGGTATCT 
#>                                            1 
#>                    G/GGAGCCACTCCCCAGGGAAACAC 
#>                                            1 
#>                                    G/GGAGTTT 
#>                                            1 
#>                           G/GGGCGGGGGGGCGGGT 
#>                                            1 
#>                                        G/GGT 
#>                                            1 
#>                                         G/GT 
#>                                           12 
#>                                        G/GTA 
#>                                            2 
#>                                    G/GTCAGGC 
#>                                            1 
#>    G/GTGGATCACCTGAGGCCAGGAGTTCAAGATTAGTCTGGC 
#>                                            1 
#>                                      G/GTGTT 
#>                                            2 
#>                                        G/GTT 
#>                                            2 
#>                                    G/GTTACTA 
#>                                            1 
#>                                   G/GTTACTTT 
#>                                            1 
#>                                      G/GTTCT 
#>                                            1 
#>                                    G/GTTTATT 
#>                                            1 
#>                                     G/GTTTTC 
#>                                            1 
#>                              G/GTTTTGTTTCGTT 
#>                                            1 
#>                                          G/T 
#>                                          498 
#>                                         GA/G 
#>                                           13 
#>                                    GAAAAGA/G 
#>                                            1 
#>                                        GAC/G 
#>                                            1 
#>                                         GC/G 
#>                                            2 
#>                                        GCC/G 
#>                                            1 
#>                                      GCTAT/G 
#>                                            1 
#>                                 GCTCCCCCTT/G 
#>                                            1 
#>                                         GT/G 
#>                                           11 
#>                                        GTA/G 
#>                                            2 
#>                                      GTCAA/G 
#>                                            1 
#>                                       GTCC/G 
#>                                            1 
#>                                          T/A 
#>                                          442 
#>                                          T/C 
#>                                         2994 
#>                                          T/G 
#>                                          723 
#>                                         T/TA 
#>                                            6 
#>                                        T/TAA 
#>                                            1 
#>                     T/TAAAAAAAAAAAAAAAAAAAAA 
#>                                            1 
#>                                      T/TAATA 
#>                                            2 
#>                                      T/TAGAA 
#>                                            1 
#>                                      T/TAGAC 
#>                                            1 
#>                                      T/TAGAG 
#>                                            2 
#>                                     T/TAGCTA 
#>                                            1 
#>                                       T/TATA 
#>                                            1 
#>                                       T/TATC 
#>                                            1 
#>                                         T/TC 
#>                                           11 
#>                                      T/TCTTC 
#>                                            1 
#>                                         T/TG 
#>                                           16 
#>                                        T/TGA 
#>                                            6 
#>                                      T/TGAGA 
#>                                            1 
#>                                      T/TGATA 
#>                                            1 
#>                                        T/TGC 
#>                                            1 
#>                                     T/TGCTGG 
#>                                            1 
#>                                        T/TGG 
#>                                            1 
#>                                      T/TGTGA 
#>                                            2 
#>                                        T/TTA 
#>                                            1 
#>                                       T/TTAA 
#>                                            1 
#>                                    T/TTAAAAA 
#>                                            1 
#>                                      T/TTATC 
#>                                            1 
#>                                        T/TTC 
#>                                            3 
#>                                      T/TTCTC 
#>                                            1 
#>                                      T/TTCTG 
#>                                            1 
#>                                        T/TTG 
#>                                            2 
#>                                      T/TTGTA 
#>                                            3 
#>                                      T/TTGTC 
#>                                            1 
#>                          T/TTGTGTGTGAGTGTGCA 
#>                                            1 
#>                                       T/TTTA 
#>                                            1 
#>                               T/TTTTCTTCTTTC 
#>                                            1 
#>                                      T/TTTTG 
#>                                            2 
#>                                         TA/T 
#>                                           11 
#>                                        TAG/T 
#>                                            1 
#>                                       TATA/T 
#>                                            1 
#>                                     TATAAA/T 
#>                                            1 
#>                                         TC/T 
#>                                            1 
#>                                        TCA/T 
#>                                            2 
#>                                   TCAGGTCC/T 
#>                                            1 
#>                    TCCCGGAGCTCCTATGACATGTA/T 
#>                                            1 
#>                                         TG/T 
#>                                            9 
#>                                        TGA/T 
#>                                            1 
#>                                   TGTGGGAG/T 
#>                                            2 
#>                                        TTA/T 
#>                                            1 
#>                                      TTAGG/T 
#>                                            1 
#>                                     TTATAA/T 
#>                                            1 
#>                                TTCCCTGAGTC/T 
#>                                            1 
#>                                        TTG/T 
#>                                            3 
#>                                      TTTCC/T 
#>                                            1 
#>                                     TTTTTG/T 
#>                                            1 
```
