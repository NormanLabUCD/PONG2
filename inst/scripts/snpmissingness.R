#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(PONG2))

# =============================================================================
# ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: snpmissingness.R <input> <output> <assembly> <locus> <filter> <PONG2_root> [model_path]")
}

input      <- args[1]
output     <- args[2]
assembly   <- args[3]
locus      <- args[4]
filter     <- as.numeric(args[5])
PONG2_root <- args[6]
# An unset $MODEL_PATH arrives as "", which must mean "use the built-in model",
# not "open the file called ''".
model_path <- if (length(args) >= 7 && !is.na(args[7]) && nzchar(args[7]) &&
                  !args[7] %in% c("Null", "NULL", "NA", "None")) args[7] else NULL

# =============================================================================
# VALIDATE ARGUMENTS
# =============================================================================
stopifnot(
  "input is missing"      = !is.na(input)      && nzchar(input),
  "output is missing"     = !is.na(output)     && nzchar(output),
  "locus is missing"      = !is.na(locus)      && nzchar(locus),
  "assembly is missing"   = !is.na(assembly)   && nzchar(assembly),
  "filter is invalid"     = !is.na(filter),
  "PONG2_root is missing" = !is.na(PONG2_root) && nzchar(PONG2_root)
)

# =============================================================================
# VALIDATE ASSEMBLY & LOCUS
# =============================================================================
assembly <- match.arg(assembly, choices = c("hg19", "hg38"))
locus    <- toupper(locus)

# =============================================================================
# LOAD MODEL — custom or built-in
# =============================================================================
if (!is.null(model_path)) {
  # ── User-supplied model ───────────────────────────────────────────────────
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
  
} else {
  # ── Built-in pre-trained model ────────────────────────────────────────────
  
  # Resolve the model store exactly as predict.R does. The previous
  # system.file("data", "Rdata.rds", ...) returns "" when the file is not in
  # the installed package, and readRDS("") then fails with
  #   cannot open compressed file ''
  # rather than saying the model store is missing. .get_model_path() also
  # covers the user-cache location that system.file never sees.
  # ::: so this works whether or not the helper is exported.
  rds_path <- PONG2:::.get_model_path()
  if (!length(rds_path) || !nzchar(rds_path) || !file.exists(rds_path)) {
    stop("PONG2 model store not found (.get_model_path() returned ",
         if (length(rds_path)) shQuote(rds_path) else "nothing", ").\n",
         "  Pass a model explicitly with --model, or reinstall PONG2 so the\n",
         "  model store is present.")
  }
  
  # readRDS returns the model list itself; predict.R uses it directly. The old
  # get(object$models) indirection belongs to an earlier RDS layout and errors
  # on the current one.
  getObject <- readRDS(rds_path)
  
  # Validate filter
  valid_filters <- c(0, 0.01, 0.005)
  if (!filter %in% valid_filters) {
    stop("filter must be one of: ", paste(valid_filters, collapse = ", "),
         " — received: ", filter)
  }
  
  filter_key <- switch(as.character(filter),
                       "0"     = "allele_fileter_00",
                       "0.01"  = "allele_fileter_001",
                       "0.005" = "allele_fileter_0005"
  )
  
  # Validate locus — derived directly from model object
  available_filter_keys <- names(getObject[[assembly]])
  
  supported_loci <- unique(unlist(
    lapply(available_filter_keys, function(fk) names(getObject[[assembly]][[fk]]))
  ))
  
  if (!locus %in% supported_loci) {
    stop(
      "Locus '", locus, "' is not available for assembly '", assembly, "'.\n",
      "Supported loci for ", assembly, ": ",
      paste(sort(supported_loci), collapse = ", ")
    )
  }
  
  # Check locus exists for the requested filter level specifically
  if (is.null(getObject[[assembly]][[filter_key]][[locus]])) {
    available_for_locus <- available_filter_keys[
      sapply(available_filter_keys, function(fk)
        !is.null(getObject[[assembly]][[fk]][[locus]])
      )
    ]
    stop(
      "Locus '", locus, "' is not available for filter=", filter,
      " in assembly '", assembly, "'.\n",
      "Available filter levels for this locus: ",
      paste(available_for_locus, collapse = ", ")
    )
  }
  
  mobj <- getObject[[assembly]][[filter_key]][[locus]]
}

# =============================================================================
# COMPUTE SNP OVERLAP
# =============================================================================
model <- hlaModelFromObj(mobj)

# KIR is chromosome 19 only. CHR must be part of the filter: a bare BP range
# matches the same coordinate on every other chromosome, which inflates the
# overlap — up to 100% — whenever the input is not already chr19-only.
KIR_CHR <- "19"

bim <- read.table(
  paste0(input, ".bim"),
  header           = FALSE,
  sep              = "",
  col.names        = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
  colClasses       = c("character", "character", "numeric",
                       "integer", "character", "character"),
  stringsAsFactors = FALSE
)

# Accept "19" or "chr19" in the .bim without silently matching nothing
bim$CHR <- sub("^chr", "", bim$CHR)
bim <- bim[bim$CHR == KIR_CHR, , drop = FALSE]

if (nrow(bim) == 0) {
  stop("no chromosome ", KIR_CHR, " variants in ", input, ".bim")
}

# Count model SNPs, not distinct positions, and divide by model SNPs — the
# same ratio predict.R reports, so the two can no longer disagree. Counting
# distinct positions in the numerator against all model SNPs in the
# denominator mixes two different quantities.
mk_id      <- function(p) paste0(KIR_CHR, ":", as.integer(p))
model_ids  <- mk_id(model$snp.position)
data_ids   <- mk_id(bim$BP)
match_rate <- sum(model_ids %in% data_ids) / length(model_ids)

# Full precision: the caller multiplies by 100 and formats. format(nsmall = 2)
# sets only a minimum number of decimals, so it was never the rounding step.
cat(sprintf("%.6f", match_rate), "\n")