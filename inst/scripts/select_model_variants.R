#!/usr/bin/env Rscript
# select_model_variants.R
#
# Writes a plink2 --extract list that resolves duplicated positions using the
# model's own alleles instead of discarding them.
#
# --rm-dup exclude-all drops every record at a position shared by more than one
# variant. At a multi-allelic site that destroys the record the model needs:
# KIR3DL1 position 54805984 has model allele T/C while the data carries both
# A/C and T/C, and exclude-all deletes both rows. The model states which row is
# correct, so the right action is to keep that one.
#
# Positions with a single record are kept unchanged. Positions with several are
# kept only when exactly one record's allele pair matches the model's; an
# unresolvable position is still dropped, because binding a model SNP to an
# arbitrary record would feed HIBAG the wrong genotypes silently.
#
# Input .bim must carry per-record unique IDs (plink2 --set-all-var-ids
# '@:#:$r:$a'), since the IDs are what plink2 --extract matches on.
#
# Usage:
#   select_model_variants.R <tagged.bim> <out_keep.txt> <assembly> <locus> \
#                           <filter> [model_path]

suppressPackageStartupMessages(library(PONG2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("usage: select_model_variants.R <tagged.bim> <out_keep.txt> <assembly> <locus> <filter> [model_path]")
}

BIM        <- args[1]
OUT        <- args[2]
assembly   <- args[3]
locus      <- toupper(args[4])
filter     <- as.numeric(args[5])
model_path <- if (length(args) >= 6 && !is.na(args[6]) && nzchar(args[6]) &&
                  !args[6] %in% c("Null", "NULL", "NA", "None")) args[6] else NULL

stopifnot(file.exists(BIM))
assembly <- match.arg(assembly, choices = c("hg19", "hg38"))

# ---- model ------------------------------------------------------------------
if (!is.null(model_path)) {
  if (!file.exists(model_path)) stop("Custom model file not found: ", model_path)
  e <- new.env()
  load(model_path, envir = e)
  if (!"mobj" %in% ls(e)) stop("Custom model file must contain an object named 'mobj'")
  mobj <- e$mobj
} else {
  rds_path <- PONG2:::.get_model_path()
  if (!length(rds_path) || !nzchar(rds_path) || !file.exists(rds_path)) {
    stop("PONG2 model store not found (.get_model_path() returned nothing usable)")
  }
  getObject  <- readRDS(rds_path)
  filter_key <- switch(as.character(filter),
                       "0"     = "allele_fileter_00",
                       "0.01"  = "allele_fileter_001",
                       "0.005" = "allele_fileter_0005",
                       stop("filter must be one of 0, 0.01, 0.005 — received: ", filter))
  mobj <- getObject[[assembly]][[filter_key]][[locus]]
  if (is.null(mobj)) {
    stop("No model for locus ", locus, " | assembly ", assembly, " | filter ", filter)
  }
}

m_pos <- as.integer(mobj$snp.position)
m_al  <- as.character(mobj$snp.allele)

# ---- data -------------------------------------------------------------------
bim <- read.table(BIM, header = FALSE, sep = "", stringsAsFactors = FALSE,
                  colClasses = c("character", "character", "numeric",
                                 "integer", "character", "character"))
names(bim) <- c("chr", "id", "cm", "pos", "a1", "a2")

# Two records sharing position AND allele pair produce the same tagged ID.
# They are genuinely indistinguishable, so keep the first and drop the rest
# rather than failing the run; --extract would match both anyway.
n_dup_id <- sum(duplicated(bim$id))
if (n_dup_id) {
  message("[select] records with a non-unique tagged ID dropped: ", n_dup_id)
  bim <- bim[!duplicated(bim$id), , drop = FALSE]
}

# allele pair as an order-independent key, so A/C and C/A compare equal
pair_key <- function(x, y) {
  x <- toupper(x); y <- toupper(y)
  ifelse(x < y, paste0(x, "/", y), paste0(y, "/", x))
}

bim$key <- pair_key(bim$a1, bim$a2)

split_al  <- strsplit(m_al, "/", fixed = TRUE)
ok        <- lengths(split_al) == 2L
model_key <- rep(NA_character_, length(m_al))
model_key[ok] <- pair_key(vapply(split_al[ok], `[`, "", 1L),
                          vapply(split_al[ok], `[`, "", 2L))

# model allele keyed by position; a position the model lists twice is ambiguous
# from the model's side too, so it is excluded from the lookup
dup_model_pos <- m_pos[duplicated(m_pos)]
usable        <- !is.na(model_key) & !(m_pos %in% dup_model_pos)
model_lookup  <- setNames(model_key[usable], as.character(m_pos[usable]))

# ---- select -----------------------------------------------------------------
n_at_pos <- ave(seq_len(nrow(bim)), bim$pos, FUN = length)

single <- n_at_pos == 1L
multi  <- !single

want_key <- model_lookup[as.character(bim$pos)]
matches  <- multi & !is.na(want_key) & bim$key == want_key

# keep one record per position even if two share position AND allele pair
first_of_pos <- !duplicated(bim$pos[matches])
resolved_ids <- bim$id[matches][first_of_pos]

keep_ids <- c(bim$id[single], resolved_ids)

writeLines(keep_ids, OUT)

# ---- report (stderr; stdout stays clean for callers) ------------------------
multi_pos     <- unique(bim$pos[multi])
resolved_pos  <- unique(bim$pos[matches])
model_multi   <- intersect(multi_pos, m_pos)


