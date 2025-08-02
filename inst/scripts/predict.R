#!/usr/bin/env Rscript
library(PONG)

source("preprocess.R")
args <- commandArgs(trailingOnly = TRUE)
action <- args[1]
geno <- args[2]
locus <- args[3]
out <- args[4]

user_dir <- path.expand("~")
dst_dir <- ""
pong_output="Pong_result"
# Check if the output directory exists and Create the output directory
if (!dir.exists(out)) {
# Create the full path to the folder
dst_dir <- file.path(user_dir, pong_output)
if (!dir.exists(dst_dir)) { dir.create(dst_dir) }
}else{
    dst_dir <- file.path(out, pong_output)
    if (!dir.exists(dst_dir)) {
    dir.create(dst_dir, recursive = TRUE)
    }
}

print("Reformating plink files to the right format")
geno_path <- data_format(geno, dst_dir)

#predict.hlaAttrBagClass()
bed.fn <- (paste0(geno_path,'.bed'))
fam.fn <- (paste0(geno_path,'.fam'))
bim.fn <- (paste0(geno_path,'.bim'))

genotype <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, import.chr='19', assembly="hg19")
print("SNP loaded..")
model_path <- system.file("models", paste0("KIR3DL1_55314840_55378697_model",".RData"), package = "PONG")


mobj <- get(load(model_path))
model <- hlaModelFromObj(mobj)

region <- 5000    # kb
geno <- hlaGenoSubsetFlank(genotype, locus, region*5000, assembly="hg19")

pred.guess <- kirPredict(model, geno, type="response+prob")

save(pred.guess, file=paste0(dst_dir, locus, "_file.RData"))
write.table(pred.guess$value, file=paste0(dst_dir, locus, ".txt"),  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
print(paste0("results saved in ",dst_dir,""))
print("Prediction done..")


#print(length(args))
# Check if arguments are provided
if (length(args) < 2) {
  stop("At least two arguments are required: --predict, --geno")
}

# Assign arguments to variables

#Downloads/Indigo/427ME_Viken061219Updated_ForKIRimpute_hg19.bim
#/Users/suraj/Downloads/indigo/427ME_Viken061219Updated_ForKIRimpute_hg19.bim