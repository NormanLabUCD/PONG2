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
model_path <- if (length(args) >= 7 && args[7] != "Null") args[7] else NULL

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

  # Load model object once
  rds_path  <- system.file("data", "Rdata.rds", package = "PONG2")
  object    <- readRDS(rds_path)
  getObject <- get(object$models)

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

min_pos <- min(model$snp.position)
max_pos <- max(model$snp.position)

# Read bim file
bim <- read.table(
  paste0(input, ".bim"),
  header           = FALSE,
  sep              = "",
  col.names        = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
  stringsAsFactors = FALSE
)

# Subset bim to model region
subset_bim <- subset(bim, BP >= min_pos & BP <= max_pos)

# Compute overlap
overlap    <- unique(intersect(model$snp.position, subset_bim$BP))
match_rate <- length(overlap) / length(model$snp.position)

cat(format(match_rate, nsmall = 2), "\n")
