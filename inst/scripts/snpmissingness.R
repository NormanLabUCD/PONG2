#!/usr/bin/env Rscript
library(PONG2)

# =============================================================================
# ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: snpmissingness.R <input> <output> <assembly> <locus> <filter> <PONG2_root>")
}

input      <- args[1]
output     <- args[2]
assembly   <- args[3]
locus      <- args[4]
filter     <- as.numeric(args[5])
PONG2_root <- args[6]

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
# LOAD MODEL OBJECT ONCE
# =============================================================================
rds_path  <- system.file("data", "Rdata.rds", package = "PONG2")
object    <- readRDS(rds_path)
getObject <- get(object$models)

# =============================================================================
# VALIDATE ASSEMBLY
# =============================================================================
assembly <- match.arg(assembly, choices = c("hg19", "hg38"))
locus    <- toupper(locus)

# =============================================================================
# VALIDATE FILTER
# =============================================================================
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

# =============================================================================
# VALIDATE LOCUS — derived directly from model object
# =============================================================================

# Get available filter keys for this assembly
available_filter_keys <- names(getObject[[assembly]])

# Get all loci available across ALL filter levels for this assembly
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

# =============================================================================
# LOAD MODEL AND COMPUTE SNP OVERLAP
# =============================================================================
mobj  <- getObject[[assembly]][[filter_key]][[locus]]
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
