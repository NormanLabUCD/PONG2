\name{PONG-package}
\alias{PONG-package}
\alias{PONG}
\docType{package}
\title{HLA Genotype Imputation with Attribute Bagging}
\description{
	To impute HLA types from unphased SNP data using an attribute bagging
method.
}
\details{
	\tabular{ll}{
		Package: \tab PONG\cr
		Type:    \tab Package\cr
		Version: \tab 1.2.4\cr
		License: \tab GPL version 3\cr
	}

	PONG is a state of the art software package for imputing HLA types using
SNP data, and it uses the R statistical programming language. PONG is highly
accurate, computationally tractable, and can be used by researchers with
published parameter estimates instead of requiring access to large training
sample datasets. It combines the concepts of attribute bagging, an ensemble
classifier method, with haplotype inference for SNPs and HLA types. Attribute
bagging is a technique which improves the accuracy and stability of classifier
ensembles using bootstrap aggregating and random variable selection.

\bold{Features:} \cr
1) PONG can be used by researchers with published parameter estimates
	(\url{http://www.biostat.washington.edu/~bsweir/PONG/}) instead of
	requiring access to large training sample datasets. \cr
2) A typical PONG parameter file contains only haplotype frequencies at
	different SNP subsets rather than individual training genotypes. \cr
3) SNPs within the xMHC region (chromosome 6) are used for imputation. \cr
4) PONG employs unphased genotypes of unrelated individuals as a training
	set. \cr
5) PONG supports parallel computing with R. \cr
}
\references{
	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging;
	(Abstract 294, Platform/Oral Talk);
	Present at the 62nd Annual Meeting of the American Society of Human Genetics,
	November 9, 2012 in San Francisco, California.

	Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS;
	HIBAG -- HLA Genotype Imputation with Attribute Bagging.
	Pharmacogenomics Journal. doi: 10.1038/tpj.2013.18.
	\url{http://dx.doi.org/10.1038/tpj.2013.18}
}
\author{
	Xiuwen Zheng \email{zhengx@u.washington.edu},
	Bruce S. Weir \email{bsweir@u.washington.edu}
}

\examples{
# load HLA types and SNP genotypes
data(HLA_Type_Table, package="PONG")
data(HapMap_CEU_Geno, package="PONG")

head(HLA_Type_Table)
dim(HLA_Type_Table)  # 60 13

summary(HapMap_CEU_Geno)


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
# "training"   "validation"
summary(hlatab$training)
summary(hlatab$validation)



}

\keyword{HLA}
\keyword{SNP}
\keyword{genetics}
