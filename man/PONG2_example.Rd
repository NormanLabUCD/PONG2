\name{PONG2_example}
\alias{PONG2_example}
\alias{example_snp}
\alias{example_kir}
\alias{example_mobj}
\title{PONG2 Example Dataset}
\description{
  A small example dataset for demonstrating PONG2 functions.
  Contains 50 samples and 200 SNPs in the KIR region (chr19),
  along with a pre-trained KIR3DL1 model with 10 classifiers.
}
\usage{
data(PONG2_example)
}
\format{
  Three objects are loaded:
  \describe{
    \item{example_snp}{A \code{hlaSNPGenoClass} object with 50 samples
      and 200 SNPs in the KIR region (chr19, hg19 assembly)}
    \item{example_kir}{A \code{hlaAlleleClass} object with KIR3DL1
      allele calls for 50 samples}
    \item{example_mobj}{A \code{hlaAttrBagObj} object — a pre-trained
      KIR3DL1 model with 10 ensemble classifiers}
  }
}
\examples{
data(PONG2_example)
# SNP data
cat("Samples:", ncol(example_snp$genotype), "\n")
cat("SNPs:   ", nrow(example_snp$genotype), "\n")
# KIR allele table
cat("Locus:  ", example_kir$locus, "\n")
# Model
cat("Classifiers:", length(example_mobj$classifiers), "\n")
}
\keyword{datasets}
