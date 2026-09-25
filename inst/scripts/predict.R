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
# VALIDATE INPUT FILES
# Fail here with a clear message rather than deep inside read.table or plink2.
# =============================================================================
for (.ext in c(".bed", ".bim", ".fam")) {
  .f <- paste0(input, .ext)
  if (!file.exists(.f)) stop("PLINK input file not found: ", .f)
}
rm(.ext, .f)

# =============================================================================
# MODEL LOADER — built-in
# =============================================================================
modelObject <- function(locus, filter = 0.005, assembly = c("hg38", "hg19")) {
  assembly <- match.arg(assembly)
  locus    <- toupper(locus)
  # Must accept the same set as snpmissingness.R, otherwise a run passes the
  # SNP-overlap check and then dies here on an unsupported --filter.
  valid_filters <- c(0, 0.01, 0.005)
  if (!filter %in% valid_filters) {
    stop("filter must be one of: ", paste(valid_filters, collapse = ", "),
         " — received: ", filter)
  }
  rds_path   <- .get_model_path()
  getObject  <- readRDS(rds_path)
  filter_key <- switch(as.character(filter),
                       "0"     = "allele_fileter_00",
                       "0.01"  = "allele_fileter_001",
                       "0.005" = "allele_fileter_0005")
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

# =============================================================================
# REBUILD VARIANT IDs FROM POSITION
# HIBAG pairs model to genotypes with match(model$snp.id, geno$snp.id) — plain
# string equality — so both sides must spell a site identically. Nothing in the
# incoming ID strings can be relied on for that: the shipped models mix bare
# chr:pos (19:55314880, hg19) with allele-bearing IDs (19:54724545:G:A, hg38)
# and a handful of source identifiers that encode no coordinate at all
# (HGSV_234687). User data arrives with rsIDs, chr-prefixed IDs, or whatever an
# imputation server produced.
#
# So IDs are not parsed, they are reconstructed as <chr>:<position> from each
# side's own snp.position, which is authoritative and present on both. This
# makes matching positional by construction and independent of naming.
#
# The rewrite is in place — same length, same order — because a HIBAG model's
# classifiers index into snp.id. Entries must never be added or removed here.
# =============================================================================
KIR_CHR <- "19"   # PONG2 is KIR-only; hlaBED2Geno below also imports chr 19

rebuild_ids <- function(pos) paste0(KIR_CHR, ":", pos)

id_shape <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("empty")
  s <- head(x, 1000)
  if      (all(grepl("^[^:]+:[0-9]+$", s)))             "chr:pos"
  else if (all(grepl("^[^:]+:[0-9]+:[^:]+:[^:]+$", s))) "chr:pos:ref:alt"
  else if (all(grepl("^rs[0-9]+$", s)))                 "rsID"
  else                                                  "mixed/other"
}

if (is.null(mobj$snp.position))
  stop("Model object has no snp.position; cannot rebuild variant IDs.")
if (length(mobj$snp.position) != length(mobj$snp.id))
  stop("Model snp.position (", length(mobj$snp.position), ") and snp.id (",
       length(mobj$snp.id), ") are not aligned; refusing to rewrite IDs.")

.shape_before <- id_shape(mobj$snp.id)
.example_from <- mobj$snp.id[1]
mobj$snp.id   <- rebuild_ids(mobj$snp.position)
.n_model_dup  <- sum(duplicated(mobj$snp.id))

cat("Model SNP IDs rebuilt from position: ", length(mobj$snp.id),
    " SNPs (was '", .shape_before, "', e.g. ", .example_from, " -> ",
    mobj$snp.id[1], ")\n", sep = "")

if (.n_model_dup > 0) {
  cat("*** WARNING:", .n_model_dup, "model SNPs share a position with another.\n")
  cat("    Their rebuilt IDs collide; HIBAG keeps the first and treats the rest\n")
  cat("    as missing. Expect that many SNPs to be unusable.\n\n")
}
rm(.shape_before, .example_from, .n_model_dup)

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

# =============================================================================
# PRE-FLIGHT: MODEL / DATA OVERLAP
# The model's IDs were rebuilt from position above, and the genotype's are
# rebuilt the same way inside the chunk loop, so overlap is checked here on the
# same basis: positions on chromosome 19. Doing it before the loop means a
# hopeless run fails immediately instead of once per chunk, and it covers direct
# calls to predict.R and custom --model files.
# =============================================================================
model_ids <- mobj$snp.id

# Read only the position column of the .bim and build IDs the same way
bim_pos <- read.table(bim.fn, header = FALSE, sep = "",
                      colClasses = c("NULL", "NULL", "NULL",
                                     "integer", "NULL", "NULL"),
                      stringsAsFactors = FALSE)[[1]]
bim_ids <- rebuild_ids(bim_pos)

n_id_hit <- sum(model_ids %in% bim_ids)

cat("--- Model / data overlap ---\n")
cat("Model SNPs:      ", length(model_ids), "  e.g. ", head(model_ids, 1), "\n", sep = "")
cat("Data variants:   ", length(bim_ids),   "  e.g. ", head(bim_ids, 1),   "\n", sep = "")
cat(sprintf("Model SNPs present: %d / %d (%.2f%%)\n",
            n_id_hit, length(model_ids), 100 * n_id_hit / length(model_ids)))
cat("----------------------------\n\n")

if (n_id_hit == 0) {
  stop("No model SNP positions are present in the data.\n",
       "  Model positions span ", min(mobj$snp.position), "-",
       max(mobj$snp.position), "\n",
       "  Data  positions span ", min(bim_pos), "-", max(bim_pos), "\n",
       "  Non-overlapping ranges like these usually mean the data and the model\n",
       "  are on different genome builds. Check --assembly (", assembly, ").")
}

# Duplicate positions are dangerous rather than merely untidy: match() returns
# the FIRST hit, so a model SNP would bind to an arbitrary variant there.
n_dup <- sum(duplicated(bim_ids))
if (n_dup > 0) {
  cat("*** WARNING:", n_dup, "variants share a position with another in", bim.fn, "\n")
  cat("    Each model SNP binds to whichever copy appears first, whose alleles\n")
  cat("    may differ from the model's.\n")
  cat("    Deduplicate with: plink2 --set-all-var-ids '@:#' --rm-dup exclude-all\n\n")
}

rm(bim_ids, bim_pos)

# Read FAM to get full sample list without loading genotypes.
# colClasses = "character" is required: without it, numeric-looking sample IDs
# are read as numbers, so "007" becomes 7 and the --keep file no longer matches
# the .fam. plink2 then keeps 0 samples and the chunk fails.
fam_data    <- read.table(fam.fn, header = FALSE,
                          colClasses = "character", stringsAsFactors = FALSE)
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

# Output directory is created up front so a failed-chunk report can be written
# even if the run ends early. dir.create is a no-op when it already exists;
# showWarnings = FALSE keeps concurrent loci from each logging a spurious
# "already exists" warning as they race to create the shared KIR/ directory.
out_dir <- file.path(output, "KIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output directory: ", out_dir)

chunk_results <- vector("list", n_chunks)
chunk_sample_ids <- vector("list", n_chunks)   # for the failed-chunk report

for (i in seq_len(n_chunks)) {
  idx_start <- (i - 1) * CHUNK_SIZE + 1
  idx_end   <- min(i * CHUNK_SIZE, n_samples)
  chunk_ids <- all_samples[idx_start:idx_end]
  chunk_sample_ids[[i]] <- chunk_ids
  
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
    
    # Rebuild IDs from position, exactly as the model's were, so match() pairs
    # the two on coordinate rather than on whatever the input called things.
    # In place: same length, same order.
    if (length(genotype$snp.position) != length(genotype$snp.id))
      stop("genotype snp.position and snp.id are not aligned for chunk ", i)
    genotype$snp.id <- rebuild_ids(genotype$snp.position)
    
    chunk_geno <- hlaGenoSubsetFlank(genotype, locus,
                                     region * 5000,
                                     assembly = assembly)
    
    # Authoritative overlap: this is the exact set kirPredict will match
    # against, whatever hlaBED2Geno did to the IDs on the way in. Reported once
    # so the number in the log is the one that actually governs prediction.
    if (i == 1) {
      n_geno_hit <- sum(model_ids %in% chunk_geno$snp.id)
      cat(sprintf("Model SNPs present in genotype object: %d / %d (%.2f%%)\n",
                  n_geno_hit, length(model_ids),
                  100 * n_geno_hit / length(model_ids)))
      if (n_geno_hit < length(model_ids)) {
        cat("    ", length(model_ids) - n_geno_hit,
            " model positions are not in the data and count as missing.\n", sep = "")
      }
    }
    
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
  # Record which samples were lost. A warning alone is easy to miss when many
  # loci are logging at once, and a silently short results file is worse than
  # a loud one.
  failed_ids   <- unlist(chunk_sample_ids[failed], use.names = FALSE)
  failed_file  <- file.path(out_dir, paste0(locus, "_failed_samples.txt"))
  writeLines(as.character(failed_ids), failed_file)
  
  cat("\n*** WARNING:", length(failed), "of", n_chunks, "chunks failed ***\n")
  cat("    Failed chunks:  ", paste(failed, collapse = ", "), "\n")
  cat("    Samples lost:   ", length(failed_ids), "of", n_samples, "\n")
  cat("    Sample IDs in:  ", failed_file, "\n")
  cat("    Results below are INCOMPLETE.\n\n")
  
  warning("Failed chunks: ", paste(failed, collapse = ", "),
          " — excluded from final results (", length(failed_ids),
          " samples; see ", failed_file, ")")
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

# Completeness check: flag a short result even when no chunk reported an error,
# so samples dropped inside kirPredict do not pass unnoticed.
n_out <- nrow(pred.guess$value)
if (length(failed) == 0 && !is.null(n_out) && n_out != n_samples) {
  cat("\n*** NOTE: results contain", n_out, "of", n_samples,
      "samples although no chunk reported an error.\n")
  cat("    Samples may have been dropped during prediction.\n\n")
}

# =============================================================================
# SAVE RESULTS
# =============================================================================
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