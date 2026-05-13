#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(PONG2))
suppressPackageStartupMessages(library(parallel))

# =============================================================================
# ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: predict.R <input> <output> <locus> <assembly> <filter> <threads>")
}

input    <- args[1]
output   <- args[2]
locus    <- args[3]
assembly <- args[4]
filter   <- as.numeric(args[5])
threads  <- as.numeric(args[6])

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
# MODEL LOADER
# =============================================================================
modelObject <- function(locus, filter = 0.005, assembly = c("hg38", "hg19")) {
  assembly <- match.arg(assembly)
  locus    <- toupper(locus)

  valid_filters <- c(0.01, 0.005)
  if (!filter %in% valid_filters) {
    stop("filter must be one of: ", paste(valid_filters, collapse = ", "),
         " — received: ", filter)
  }

  rds_path  <- system.file("data", "Rdata.rds", package = "PONG2")
  object    <- readRDS(rds_path)
  getObject <- get(object$models)

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
# LOAD GENOTYPE DATA
# =============================================================================
cat("\n--- Prediction Settings ---\n")
cat("Input:    ", input,    "\n")
cat("Output:   ", output,   "\n")
cat("Locus:    ", locus,    "\n")
cat("Assembly: ", assembly, "\n")
cat("Filter:   ", filter,   "\n")
cat("Threads:  ", threads,  "\n")
cat("---------------------------\n\n")

mobj  <- modelObject(locus, filter, assembly)
model <- hlaModelFromObj(mobj)

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
# PREDICTION — parallel + auto-chunked when n_samples > 2000
# =============================================================================
n_samples <- length(geno$sample.id)
cat(paste0("Total samples: ", n_samples, "\n"))

# Create cluster once — shared across all chunks
# on.exit guarantees cluster stops even if an error occurs
cl <- makeCluster(threads)
on.exit(stopCluster(cl))

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

#cat(paste0("\nResults saved in: ", out_dir, "\n"))
#cat("Prediction complete.\n")
