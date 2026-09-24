#!/usr/bin/env Rscript
# Build the inputs for GRS scoring from a SNP input file:
#   - snp_list.txt           : one rsID per line (for PLINK2 --extract)
#   - copa_unweighted_score.tsv : SNP  A1  1        (always)
#   - copa_weighted_score.tsv   : SNP  A1  WEIGHT   (only when mode needs it)
#   - copa_score_variants_qc.tsv: the parsed/filtered table (QC record)
#
# Input may be an .xlsx (e.g. the GWAS weight workbook) or a .csv/.tsv/.txt
# candidate list, optionally gzipped (e.g. a PGS Catalog scoring file). Required columns (auto-detected): an rsID/SNP column and a risk/effect
# allele column. A weight column is only required for weighted scoring.
#
# Usage:
#   Rscript prepare_copa_score_files.R --input <file> [--mode both|weighted|unweighted]
#       [--outdir DIR] [--sheet NAME] [--skip N]

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

# --weights kept as an alias for backward compatibility.
input  <- get_arg("--input", get_arg("--weights", NULL))
if (is.null(input)) stop("Provide --input <xlsx|csv|tsv>")
mode   <- tolower(get_arg("--mode", "both"))
if (!mode %in% c("both", "weighted", "unweighted")) {
  stop("--mode must be one of: both, weighted, unweighted")
}
outdir <- get_arg("--outdir", "score_files")
sheet_arg <- get_arg("--sheet", NULL)
skip_n <- suppressWarnings(as.integer(get_arg("--skip", NA)))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read input ----
is_xlsx <- grepl("\\.xlsx?$", input, ignore.case = TRUE)
if (is_xlsx) {
  library(readxl)
  if (is.na(skip_n)) skip_n <- 2  # the supp. table has 2 header rows above the column names
  sheets <- excel_sheets(input)
  sheet <- if (!is.null(sheet_arg)) {
    sheet_arg
  } else if ("Supp. Table 28 RiskScoreWeights" %in% sheets) {
    "Supp. Table 28 RiskScoreWeights"
  } else {
    sheets[[1]]
  }
  raw <- read_excel(input, sheet = sheet, skip = skip_n)
} else {
  # Default: skip leading "#" metadata lines (e.g. PGS Catalog scoring-file
  # headers). file() reads .gz transparently.
  if (is.na(skip_n)) {
    head_lines <- readLines(input, n = 1000, warn = FALSE)
    skip_n <- match(FALSE, startsWith(head_lines, "#"), nomatch = length(head_lines) + 1) - 1
  }
  raw <- if (grepl("\\.csv(\\.gz)?$", input, ignore.case = TRUE)) {
    read_csv(input, skip = skip_n, show_col_types = FALSE)
  } else {
    read_tsv(input, skip = skip_n, show_col_types = FALSE)
  }
}

# ---- detect columns ----
norm <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
pick <- function(nms, patterns) {
  n <- norm(nms)
  for (p in patterns) {                      # exact normalized match first
    hit <- which(n == p)
    if (length(hit)) return(nms[hit[1]])
  }
  for (p in patterns) {                      # then substring
    hit <- grep(p, n)
    if (length(hit)) return(nms[hit[1]])
  }
  NULL
}

nms <- names(raw)
id_col  <- pick(nms, c("markername", "rsid", "snpid", "snp", "id"))
a1_col  <- pick(nms, c("riskallele", "effectallele", "a1", "ea"))
wt_col  <- pick(nms, c("weight", "beta"))
chr_col <- pick(nms, c("chrom", "chr", "chromosome"))
pos_col <- pick(nms, c("positionb37", "positionb38", "position", "pos", "bp"))

if (is.null(id_col)) {
  stop("Could not find an rsID/SNP column. Columns present: ", paste(nms, collapse = ", "))
}
if (is.null(a1_col)) {
  stop("Could not find a risk/effect allele column. Columns present: ", paste(nms, collapse = ", "))
}

need_weighted <- mode %in% c("both", "weighted")
if (need_weighted && is.null(wt_col)) {
  stop("Mode '", mode, "' requires weights, but no weight/beta column was found. ",
       "Use --mode unweighted, or provide a weight column. Columns present: ",
       paste(nms, collapse = ", "))
}

# ---- build + filter ----
df <- tibble(
  SNP = as.character(raw[[id_col]]),
  A1  = toupper(as.character(raw[[a1_col]]))
)
df$CHR <- if (!is.null(chr_col)) sub("^chr", "", as.character(raw[[chr_col]]), ignore.case = TRUE) else NA_character_
df$POS <- if (!is.null(pos_col)) as.character(raw[[pos_col]]) else NA_character_
if (!is.null(wt_col)) df$WEIGHT <- suppressWarnings(as.numeric(raw[[wt_col]]))

df <- df %>%
  filter(!is.na(SNP), SNP != "", A1 %in% c("A", "C", "G", "T"))
if (need_weighted) df <- df %>% filter(!is.na(WEIGHT))
df <- df %>% distinct(SNP, .keep_all = TRUE)

if (nrow(df) == 0) stop("No valid SNPs after filtering - check the input columns.")

# ---- write outputs ----
write_lines(df$SNP, file.path(outdir, "snp_list.txt"))
write_tsv(tibble(SNP = df$SNP, A1 = df$A1, WEIGHT = 1),
          file.path(outdir, "copa_unweighted_score.tsv"))
if (need_weighted) {
  write_tsv(tibble(SNP = df$SNP, A1 = df$A1, WEIGHT = df$WEIGHT),
            file.path(outdir, "copa_weighted_score.tsv"))
}
# QC table (also the "expected" list for the coverage report): SNP, CHR, POS, A1[, WEIGHT]
write_tsv(df %>% select(SNP, CHR, POS, A1, any_of("WEIGHT")),
          file.path(outdir, "copa_score_variants_qc.tsv"))

message(sprintf("Detected columns: rsID='%s', allele='%s'%s",
                id_col, a1_col,
                if (is.null(wt_col)) "" else sprintf(", weight='%s'", wt_col)))
message(sprintf("Wrote %d variants to %s (mode=%s)", nrow(df), outdir, mode))
