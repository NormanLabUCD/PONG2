library(PONG2)
library(readr)
library(parallel)
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
chr19       <- args[1]
kirfile     <- args[2]
locus       <- args[3]
assembly    <- args[4]
out         <- args[5]
nclassifier <- args[6]
threads     <- args[7]
kirSplit    <- args[8]
kirmaf      <- args[9]
mac         <- args[10]
locusregion <- args[11]

# Handle Defaults and Numeric Conversions
if(nclassifier == "Null"){
  nclassifier <- 100
} else {
  nclassifier <- as.numeric(nclassifier)
}

if(kirSplit == "Null"){
  kirSplit <- 0.8
} else {
  kirSplit <- as.numeric(kirSplit)
}

if(kirmaf == "Null"){
  kirmaf <- 0.00
} else {
  kirmaf <- as.numeric(kirmaf)
}

if(mac == "Null"){
  mac <- 3
} else {
  mac <- as.numeric(mac)
}

# Cluster setup using the threads argument
if(threads == "Null"){
  num_threads <- 4
} else {
  num_threads <- as.numeric(threads)
}
cluster <- makeCluster(num_threads)

# Print Summary
cat("\n--- Parameter Settings ---\n")
cat("Chromosome 19: ", chr19, "\n")
cat("KIR File:      ", kirfile, "\n")
cat("Locus:         ", locus, "\n")
cat("Assembly:      ", assembly, "\n")
cat("Output:        ", out, "\n")
cat("N Classifier:  ", nclassifier, "\n")
cat("Threads:       ", num_threads, "\n")
cat("KIR Split:     ", kirSplit, "\n")
cat("KIR MAF:       ", kirmaf, "\n")
cat("Genotype MAC:  ", mac, "\n")
cat("--------------------------\n\n")


#pong2 train --bfile ~/projects/PONG2/tests/example/chr19 --kfile ~/projects/PONG2/tests/example/kir_file.csv --output ~/projects/ --locus KIR2DL1 --assembly hg19


# Best performing positions (Default)
if(assembly == "hg19" && locusregion == "Null"){
  kir_position = list(`KIR3DL3` = c(55235681, 55248171), `KIR3DL2` = c(55361663, 55378697
  ), `KIR3DL1` = c(55327689, 55378569), `KIR2DL5A` = c(55271881, 55371344
  ), `KIR2DL5B` = c(55265986, 55337447), `KIR2DL4` = c(55314840, 55326052
  ), `KIR2DL3` = c(55249743, 55264553), `KIR2DL2` = c(55249711, 55264528
  ), `KIR2DL1` = c(55281035, 55295784), `KIR2DS1` = c(55281035, 55295755
  ), `KIR2DS2` = c(55249711, 55264292), `KIR2DS3` = c(55281003, 55296109
  ), `KIR2DS4` = c(55290062, 55360046), `KIR2DS5` = c(55281035, 55296300
  ), `KIR3DS1` = c(55327478, 55342622), `KIR2DP1` = c(55266208, 55279344
  ), `KIR3DP1` = c(55297540, 55303593), `KIR3DL1S1` = c(55327478, 55378569),
  `KIR2DL23` = c(55249743, 55264553))
  region = kir_position[[locus]]

}else if(assembly == "hg38" && locusregion == "Null"){
  kir_position = list(`KIR3DL3` = c(54724442, 54736632), `KIR3DL2` = c(54850443, 54867207
  ), `KIR3DL1` = c(54816468, 54830778), `KIR2DL5A` = c(55271881, 55371344
  ), `KIR2DL5B` = c(55265986, 55337447), `KIR2DL4` = c(54803610, 54814517
  ), `KIR2DL3` = c(54738513,54753052), `KIR2DL2` = c(55249711, 55264528
  ), `KIR2DL1` = c(54769793, 54784322), `KIR2DS1` = c(55281035, 55295755
  ), `KIR2DS2` = c(55249711, 55264292), `KIR2DS3` = c(55281003, 55296109
  ), `KIR2DS4` = c(54832676, 54848569), `KIR2DS5` = c(55281035, 55296300
  ), `KIR3DS1` = c(55327478, 55342622), `KIR2DP1` = c(55266208, 55279344
  ), `KIR3DP1` = c(54786362,54790417))
  region = kir_position[[locus]]
}else{
  region = as.numeric(unlist(regmatches(locusregion, gregexpr("\\d+", locusregion)))[-1])
}

frequent_filtered <- function(kir_data, locus, filtered_freq){
  allele1=paste0(locus, "_h1")
  allele2=paste0(locus, "_h2")
  kir_data = kir_data[, c("Sample", allele1, allele2)]

  for (col in names(kir_data)) {
    kir_data[[col]] <- gsub(paste0(locus, "\\*null"), paste0(locus, "\\*0000"), kir_data[[col]])
  }
  kir_data <-  kir_data[!apply(kir_data, 1, function(row) any(row == paste0(locus, "*new"))), ]
  kir_data <-  kir_data[!apply(kir_data, 1, function(row) any(row == paste0(locus, "*unresolved"))), ]

  allele_freqs <- kir_data %>%
    select(all_of(c(allele1, allele2))) %>%
    pivot_longer(everything()) %>%
    count(value) %>%
    mutate(freq = n / sum(n)) %>%
    filter(freq > filtered_freq)

  filtered_data <- kir_data[
    kir_data[[allele1]] %in% allele_freqs$value & kir_data[[allele2]] %in% allele_freqs$value,
  ]
  return(filtered_data)
}


filterSNP <- paste0(out, "/kir")
system(paste("plink2 --bfile", chr19," --mac", mac, "--silent --make-bed --out", filterSNP))
bed.fn <- paste0(filterSNP, '.bed')
fam.fn <- paste0(filterSNP, '.fam')
bim.fn <- paste0(filterSNP, '.bim')
geno19 <- hlaBED2Geno(bed.fn, fam.fn, bim.fn, import.chr='19', assembly=assembly, verbose=TRUE)

allele1 = paste0(locus, "_h1")
allele2 = paste0(locus, "_h2")

KIR_df <- read_csv(kirfile, show_col_types = FALSE)
type_filtered <- frequent_filtered(KIR_df, locus, filtered_freq=kirmaf)

KIR_type <- hlaAllele(
  type_filtered$Sample,
  H1=type_filtered[[allele1]],
  H2=type_filtered[[allele2]],
  locus=locus,
  assembly=assembly
)


# Load Best performing KIR regions
#load_postion <- dget("kir_positions.txt")

#Split into training and validation sets
kirtab <- hlaSplitAllele(KIR_type, train.prop=kirSplit)
#print(paste("KIR locus:", locus,", min:", min(region)))
#cat(paste("KIR locus:", locus,", max:", max(region)))

# Select best KIR region
kir_geno_pos <- geno19$snp.position[geno19$snp.position >= min(region) & geno19$snp.position <= max(region)]
kir_geno <- hlaGenoSubset(geno19, snp.sel= unique(match(kir_geno_pos, geno19$snp.position)))
snpId <- kir_geno$snp.id

# Divide into train and test sets
train.geno <- hlaGenoSubset(geno19, snp.sel = match(snpId, geno19$snp.id), samp.sel = na.omit(match(kirtab$training$value$sample.id, geno19$sample.id)))
test.geno <- hlaGenoSubset(geno19, samp.sel=na.omit(match(kirtab$validation$value$sample.id, geno19$sample.id)))

cat("\n--- Training Settings ---\n")
#cat("Training samples:  ", length(train_samples), "\n")
cat(paste0("SNPs in", locus," Region:    "), length(kir_geno_pos), "\n")
#cat("SNP IDs:           ", length(snpId), "\n")
cat("Train Samples:", length(train.geno$sample.id), "\n")
cat("Test Samples:   ", length(test.geno$sample.id), "\n")
cat("--------------------------\n\n")

# train a KIR model
set.seed(100)
model <- kirParallelAttrBagging(cl=cluster, kirtab$training, train.geno, nclassifier=nclassifier)


#save results
mobj <- hlaModelToObj(model)
save(mobj, file=paste0(out, "/", locus,"_model.RData"))
save(test.geno, file=paste0(out, "/", locus, "_test.RData"))
save(kirtab, file=paste0(out, "/", locus, "_split.RData"))

stopCluster(cluster)
