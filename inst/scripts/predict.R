#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(PONG2))

# =============================================================================
# ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: predict.R <input> <output> <locus> <assembly> <filter> <threads> [model_path]")
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

  rds_path  <- .get_model_path()
  getObject <- readRDS(rds_path)

  filter_key <- ifelse(filter == 0.01, "allele_fileter_001", "allele_fileter_0005")
  mobj <- getObject[[assembly]][[filter_key]][[locus]]

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
cat("Input:      ", input,                               "\n")
cat("Output:     ", output,                              "\n")
cat("Locus:      ", locus,                               "\n")
cat("Assembly:   ", assembly,                            "\n")
cat("Filter:     ", filter,                              "\n")
cat("Threads:    ", threads,                             "\n")
cat("Model:      ", if (is.null(model_path))
  "Built-in (PONG2 package)"
  else
    model_path,                        "\n")
cat("---------------------------\n\n")

# =============================================================================
# LOAD MODEL — custom or built-in
# =============================================================================
if (!is.null(model_path)) {
  # ── User-supplied model ───────────────────────────────────────────────────
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
  # ── Built-in pre-trained model ────────────────────────────────────────────
  cat("Loading built-in PONG2 model...\n\n")
  mobj <- modelObject(locus, filter, assembly)
}

model <- hlaModelFromObj(mobj)

# =============================================================================
# LOAD GENOTYPE DATA
# =============================================================================
bed.fn <- paste0(input, ".bed")
fam.fn <- paste0(input, ".fam")
bim.fn <- paste0(input, ".bim")

region   <- 5000
genotype <- hlaBED2Geno(bed.fn, fam.fn, bim.fn,
                        import.chr = "19",
                        assembly   = assembly)

geno <- hlaGenoSubsetFlank(genotype, locus,
                           region * 5000,
                           assembly = assembly)

# =============================================================================
# SET THREADING
# RcppParallel controls HIBAG internal C++ threads — no makeCluster needed
# RcppParallel uses shared memory threads — handle stays valid
# =============================================================================
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = threads)
  cat(paste0("Threading: RcppParallel (", threads, " threads)\n"))
} else {
  cat("Threading: RcppParallel not available — using default threading\n")
  cat("Tip: install.packages('RcppParallel') for explicit thread control\n")
}

# =============================================================================
# PREDICTION — auto-chunked when n_samples > 2000
# =============================================================================
n_samples <- length(geno$sample.id)
cat(paste0("Total samples: ", n_samples, "\n"))

if (n_samples <= CHUNK_SIZE) {
  # ── No chunking needed ────────────────────────────────────────────────────
  cat("Running prediction (no chunking required)...\n\n")
  pred.guess <- kirPredict(model, geno, type = "response+prob")

} else {
  # ── Auto chunking ─────────────────────────────────────────────────────────
  n_chunks <- ceiling(n_samples / CHUNK_SIZE)
  cat(paste0("Large sample detected (n=", n_samples, ") — ",
             "auto-chunking into ", n_chunks,
             " chunks of ", CHUNK_SIZE, "\n\n"))

  chunk_results <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {

    idx_start <- (i - 1) * CHUNK_SIZE + 1
    idx_end   <- min(i * CHUNK_SIZE, n_samples)
    chunk_ids <- geno$sample.id[idx_start:idx_end]

    cat(sprintf("[Chunk %d/%d] Samples %d-%d (%d samples)...\n",
                i, n_chunks, idx_start, idx_end, length(chunk_ids)))

    chunk_geno <- hlaGenoSubset(
      geno,
      samp.sel = match(chunk_ids, geno$sample.id)
    )

    chunk_results[[i]] <- tryCatch({
      kirPredict(model, chunk_geno, type = "response+prob")
    }, error = function(e) {
      warning(sprintf("Chunk %d failed: %s", i, e$message))
      NULL
    })

    cat(sprintf("[Chunk %d/%d] Done\n", i, n_chunks))
  }

  # ── Remove failed chunks ──────────────────────────────────────────────────
  failed <- which(sapply(chunk_results, is.null))
  if (length(failed) > 0) {
    warning("Failed chunks: ", paste(failed, collapse = ", "),
            " — excluded from final results")
    chunk_results <- chunk_results[!sapply(chunk_results, is.null)]
  }

  if (length(chunk_results) == 0) {
    stop("All chunks failed — no results to save")
  }

  # ── Combine chunk results ─────────────────────────────────────────────────
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

