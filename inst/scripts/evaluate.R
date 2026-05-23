#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(PONG2))

args      <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: evaluate.R <model_dir> <locus> <threshold>")
}

model_dir <- args[1]
locus     <- args[2]
threshold <- as.numeric(args[3])

# ── Check required files exist ───────────────────────────────────────────────
path       <- file.path(model_dir, locus)
model_file <- paste0(path, "_model.RData")
test_file  <- paste0(path, "_test.RData")
split_file <- paste0(path, "_split.RData")

if (!file.exists(model_file))
  stop("Model file not found: ", model_file)

if (!file.exists(test_file))
  stop("Test genotype file not found: ", test_file,
       "\nHint: re-run training with --split < 1 to generate a test set")

if (!file.exists(split_file))
  stop("Split file not found: ", split_file,
       "\nHint: re-run training with --split < 1 to generate a split object")

# ── Load saved objects ────────────────────────────────────────────────────────
cat("Loading model files...\n")
mobj      <- get(load(model_file))
test.geno <- get(load(test_file))
kirtab    <- get(load(split_file))
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
  out_file <- file.path(model_dir, paste0(locus, "_eval_summary.txt"))
  detail <- hlaReport(comp, out_file, type="txt")
}
cat(paste0("\nEvaluation saved to: ", out_file, "\n"))
