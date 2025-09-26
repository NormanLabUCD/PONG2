PONG2 🧬
KIR Genotype Imputation and Analysis Toolkit

Accurate KIR genotype imputation from phased sequencing data

🚀 Quick Start

             
📋 Requirements
Phased imputation data (chromosome 19 KIR region)
PLINK-formatted input files/vcf file
R 4.0+ with essential packages
plink 2.0
minimac4 4.16

🎯 Supported KIR Loci
 KIR2DL1, KIR2DL3, KIR2DL4, KIR3DL1, KIR3DL2, KIR2DS1, KIR2DS2, KIR2DS4, KIR2DS5

# Install from source
R CMD INSTALL PONG2_1.0.0.tar.gz
R CMD INSTALL --library=path PONG2_1.0.0.tar.gz

# Or from R
install.packages("PONG2_1.0.0.tar.gz", repos = NULL, type = "source")


💡 Usage Examples
bash
# Standard analysis
pong2 impute -i sample_data -o results -l KIR3DL1 -a hg19




