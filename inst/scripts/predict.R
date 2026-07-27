#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(PONG2))
# =============================================================================
# ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: predict <input> <output> <locus> <assembly> <filter> <threads> [model_path]")
}
input      <- args[1]
output     <- args[2]
locus      <- args[3]
assembly   <- args[4]
filter     <- as.numeric(args[5])
threads    <- as.numeric(args[6])
model_path <- if (length(args) >= 7 && args[7] != "Null") args[7] else NULL
# Chunk size is fixed internally — not user configurable
# Chunking applied automatically when n_samples > CHUNK_SIZE
CHUNK_SIZE <- 2000

# =============================================================================
# VALIDATE ARGUMENTS
# =============================================================================
stopifnot(
  "input is missing"    = !is.na(input)    && nzchar(input),
  "output is missing"   = !is.na(output)   && nzchar(output),
  "locus is missing"    = !is.na(locus)    && nzchar(locus),
  "assembly is missing" = !is.na(assembly) && nzchar(assembly),
  "filter is invalid"   = !is.na(filter),
  "threads is invalid"  = !is.na(threads)
)

# =============================================================================
# MODEL LOADER — built-in
# =============================================================================
modelObject <- function(locus, filter = 0.005, assembly = c("hg38", "hg19")) {
  assembly <- match.arg(assembly)
  locus    <- toupper(locus)
  valid_filters <- c(0.01, 0.005)
  if (!filter %in% valid_filters) {
    stop("filter must be one of: ", paste(valid_filters, collapse = ", "),
         " — received: ", filter)
  }
  rds_path   <- .get_model_path()
  getObject  <- readRDS(rds_path)
  filter_key <- ifelse(filter == 0.01, "allele_fileter_001", "allele_fileter_0005")
  mobj       <- getObject[[assembly]][[filter_key]][[locus]]
  if (is.null(mobj)) {
    stop("No model found for locus: ", locus,
         " | assembly: ", assembly,
         " | filter: ", filter)
  }
  return(mobj)
}

# =============================================================================
# PRINT SETTINGS
# =============================================================================
cat("\n--- Prediction Settings ---\n")
cat("Input:      ", input,    "\n")
cat("Output:     ", output,   "\n")
cat("Locus:      ", locus,    "\n")
cat("Assembly:   ", assembly, "\n")
cat("Filter:     ", filter,   "\n")
cat("Threads:    ", threads,  "\n")
cat("Model:      ", if (is.null(model_path))
  "Built-in (PONG2 package)" else model_path, "\n")
cat("---------------------------\n\n")

# =============================================================================
# LOAD MODEL — custom or built-in
# =============================================================================
if (!is.null(model_path)) {
  cat("Loading custom model:", model_path, "\n")
  if (!file.exists(model_path)) {
    stop("Custom model file not found: ", model_path)
  }
  local_env <- new.env()
  load(model_path, envir = local_env)
  if (!"mobj" %in% ls(local_env)) {
    stop(
      "Custom model file must contain an object named 'mobj'.\n",
      "Save your model with:\n",
      "  mobj <- hlaModelToObj(model)\n",
      "  save(mobj, file = '", model_path, "')"
    )
  }
  mobj <- local_env$mobj
  cat("Custom model loaded successfully\n\n")
} else {
  cat("Loading built-in PONG2 model...\n\n")
  mobj <- modelObject(locus, filter, assembly)
}
model <- hlaModelFromObj(mobj)

# =============================================================================
# SET THREADING
# =============================================================================
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = threads)
  cat(paste0("Threading: RcppParallel (", threads, " threads)\n"))
} else {
  cat("Threading: RcppParallel not available — using default threading\n")
  cat("Tip: install.packages('RcppParallel') for explicit thread control\n")
}

# =============================================================================
# LOAD GENOTYPE DATA & PREDICTION — unified chunking
# Both hlaBED2Geno and kirPredict run per chunk so memory stays flat
# throughout — essential for WGS data (83k+ SNPs) and large biobank cohorts
# =============================================================================
bed.fn <- paste0(input, ".bed")
fam.fn <- paste0(input, ".fam")
bim.fn <- paste0(input, ".bim")
region <- 5000

# Read FAM to get full sample list without loading genotypes
fam_data    <- read.table(fam.fn, header = FALSE)
all_samples <- fam_data[[2]]
n_samples   <- length(all_samples)
n_chunks    <- ceiling(n_samples / CHUNK_SIZE)

cat(paste0("Total samples: ", n_samples, "\n"))

if (n_chunks == 1) {
  cat("Running prediction (no chunking required)...\n\n")
} else {
  cat(paste0("Large dataset detected (n=", n_samples, ") — ",
             "auto-chunking into ", n_chunks,
             " chunks of up to ", CHUNK_SIZE, " samples\n\n"))
}

chunk_results <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  idx_start <- (i - 1) * CHUNK_SIZE + 1
  idx_end   <- min(i * CHUNK_SIZE, n_samples)
  chunk_ids <- all_samples[idx_start:idx_end]
  
  if (n_chunks > 1) {
    cat(sprintf("[Chunk %d/%d] Samples %d-%d (%d samples)...\n",
                i, n_chunks, idx_start, idx_end, length(chunk_ids)))
  }
  
  # Write keep file (FID IID) for plink2 --keep
  tmp_keep   <- tempfile()
  tmp_prefix <- tempfile()
  
  chunk_results[[i]] <- tryCatch({
    # Write keep file — plink2 needs FID + IID
    keep_df <- fam_data[fam_data[[2]] %in% chunk_ids, c(1, 2)]
    write.table(keep_df, tmp_keep,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Use plink2 to extract chunk into temporary BED/FAM/BIM
    plink2_cmd <- sprintf(
      "plink2 --bfile %s --keep %s --make-bed --out %s --silent",
      shQuote(input), shQuote(tmp_keep), shQuote(tmp_prefix)
    )
    ret <- system(plink2_cmd)
    if (ret != 0) stop("plink2 extraction failed for chunk ", i)
    
    # Now load the properly subsetted BED
    genotype   <- hlaBED2Geno(
      paste0(tmp_prefix, ".bed"),
      paste0(tmp_prefix, ".fam"),
      paste0(tmp_prefix, ".bim"),
      import.chr = "19",
      assembly   = assembly
    )
    chunk_geno <- hlaGenoSubsetFlank(genotype, locus,
                                     region * 5000,
                                     assembly = assembly)
    pred <- kirPredict(model, chunk_geno, type = "response+prob")
    
    if (n_chunks > 1)
      cat(sprintf("[Chunk %d/%d] Done\n", i, n_chunks))
    
    pred
  }, error = function(e) {
    warning(sprintf("Chunk %d failed: %s", i, e$message))
    NULL
  }, finally = {
    # Clean up all temp files
    unlink(tmp_keep)
    unlink(paste0(tmp_prefix, c(".bed", ".bim", ".fam", ".log")))
  })
}

# =============================================================================
# COMBINE RESULTS
# =============================================================================
failed <- which(sapply(chunk_results, is.null))
if (length(failed) > 0) {
  warning("Failed chunks: ", paste(failed, collapse = ", "),
          " — excluded from final results")
  chunk_results <- chunk_results[!sapply(chunk_results, is.null)]
}
if (length(chunk_results) == 0) {
  stop("All chunks failed — no results to save")
}

if (n_chunks == 1) {
  pred.guess <- chunk_results[[1]]
} else {
  cat("\nCombining chunk results...\n")
  pred.guess <- list(
    value = do.call(rbind, lapply(chunk_results, function(x) x$value)),
    prob  = do.call(rbind, lapply(chunk_results, function(x) x$prob))
  )
  cat(paste0("Combined: ", nrow(pred.guess$value),
             " / ", n_samples, " samples\n"))
}

# =============================================================================
# SAVE RESULTS
# =============================================================================
out_dir <- file.path(output, "KIR")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save(pred.guess,
     file = file.path(out_dir, paste0(locus, ".RData")))

write.table(pred.guess$value,
            file      = file.path(out_dir, paste0(locus, ".csv")),
            row.names = FALSE,
            col.names = TRUE,
            sep       = ",",
            quote     = FALSE)

cat("\nImputation complete. Results saved to:", out_dir, "\n")

# =============================================================================
# CLEAN UP
# =============================================================================
hlaClose(model)