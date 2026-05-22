#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(PONG2))

args      <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: evaluate.R <model_dir> <locus> <threshold>")
}

model_dir <- args[1]
locus     <- args[2]
threshold <- as.numeric(args[3])

# Load saved objects
path      <- file.path(model_dir, locus)
mobj      <- get(load(paste0(path, "_model.RData")))
test.geno <- get(load(paste0(path, "_test.RData")))
kirtab    <- get(load(paste0(path, "_split.RData")))
model     <- hlaModelFromObj(mobj)

# Predict on test set
cat("\n--- Model Evaluation ---\n")
pred <- kirPredict(model, test.geno, type = "response+prob", verbose = FALSE)
comp <- hlaCompareAllele(kirtab$validation, pred,
                         allele.limit = model, call.threshold = threshold)

# Overall accuracy
cat(sprintf("Locus:              %s\n",     locus))
cat(sprintf("Test samples:       %d\n",     comp$overall$n.call))
cat(sprintf("Call threshold:     %.2f\n",   threshold))
cat(sprintf("Haplotype accuracy: %.1f%%\n", comp$overall$acc.haplo * 100))
cat(sprintf("Genotype accuracy:  %.1f%%\n", comp$overall$acc.geno  * 100))
cat(sprintf("Call rate:          %.1f%%\n", comp$overall$call.rate  * 100))
cat("------------------------\n")

# Per-allele accuracy
if (!is.null(comp$detail)) {
  cat("\nPer-allele accuracy:\n")
  detail <- comp$detail[order(comp$detail$acc.haplo), ]
  print(detail[, c("allele", "valid.num", "acc.haplo", "sensitivity", "specificity")])
}

# Save summary
out_file <- file.path(model_dir, paste0(locus, "_eval_summary.csv"))
write.csv(comp$detail, file = out_file, row.names = FALSE)
cat(paste0("\nEvaluation saved to: ", out_file, "\n"))
