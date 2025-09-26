#!/usr/bin/env Rscript
library(PONG2)
args <- commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]
locus = args[3]
assembly = args[4]
filter = args[5]
threads = args[6]

# Verify arguments (critical!)
stopifnot(
  !is.na(input),
  !is.na(output),
  !is.na(locus),
  !is.na(threads)
)

modelObject <- function(locus, filter=0.005, assembly = c("hg38", "hg19")){

  assembly <- match.arg(assembly)
  locus <- toupper(locus)
  valid_filters <- c(0.01, 0.005)
  if (!filter %in% valid_filters) {
    stop("filter must be one of: ", paste(valid_filters, collapse = ", "))
  }

  rds_path <- system.file("data", "Rdata.rds", package = "PONG2")
  object <- readRDS(rds_path)
  getObject <- get(object$models)

  if(assembly=="hg19"){
    if (filter==0.01){
    mobj = getObject$hg19$allele_fileter_001[[locus]]
    }else{
      mobj = getObject$hg19$allele_fileter_0005[[locus]]
    }
  }

  if(assembly=="hg38"){
    if (filter==0.01){
      mobj = getObject$hg38$allele_fileter_001[[locus]]
    }else{
      mobj = getObject$hg38$allele_fileter_0005[[locus]]
    }
  }

  return(mobj)
  }

# input="/home/suraju/out2/tmp/chr19_freeze3_RISEMS_051623_deID.nodup"
# locus="KIR3DL2"
# filter="0.005"
# assembly="hg38"

mobj = modelObject(locus, filter, assembly)
model <- hlaModelFromObj(mobj)

bed.fn <- paste0(input,'.bed')
fam.fn <- paste0(input,'.fam')
bim.fn <- paste0(input,'.bim')

region <- 5000    # kb
genotype <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, import.chr='19', assembly=assembly)
geno <- hlaGenoSubsetFlank(genotype, locus, region*5000, assembly=assembly)
pred.guess <- kirPredict(model, geno, type="response+prob")

dir.create(output, "/KIR/", recursive = TRUE)
save(pred.guess, file=paste0(output, "/KIR/", locus, ".RData"))
write.table(pred.guess$value, file=paste0(output, "/KIR/", locus, ".csv"), row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
print(paste0("results saved in ",output, "/KIR/"))
print("Prediction done..")

