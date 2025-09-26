#!/usr/bin/env Rscript
library(PONG2)
args <- commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]
assembly = args[3]
locus = args[4]
filter = args[5]
PONG2_root =args[6]

# Verify arguments (critical!)
stopifnot(
  !is.na(input),
  !is.na(output),
  !is.na(locus),
  !is.na(assembly),
  !is.na(PONG2_root)
)

modelObject <- function(locus, filter=0.005, assembly = c("hg38", "hg19")){

  assembly <- match.arg(assembly)
  locus <- toupper(locus)
  valid_filters <- c(0.01, 0.005)
  if (!filter %in% valid_filters) {
    stop("filter must be one of: ", paste(valid_filters, collapse = ", "))
  }

  rds_path <- system.file("data", "Rdata.rds", package = "PONG2")
  #rds_path <- paste0(PONG2_root,"/data/", "Rdata.rds")
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

mobj = modelObject(locus, filter, assembly)
model <- hlaModelFromObj(mobj)

min = min(model$snp.position)
max = max(model$snp.position)

bim <- read.table(paste0(input, ".bim"), header = FALSE, sep = "",
                  col.names = c("CHR","SNP","CM","BP","A1","A2"),
                  stringsAsFactors = FALSE)
subset_bim <- subset(bim, BP >= min & BP <= max)

s1=model$snp.position
s2=subset_bim$BP
s <- unique(intersect(model$snp.position, subset_bim$BP))

cat(format(length(s)/length(s1), nsmall = 2), "\n") # print to stdout

