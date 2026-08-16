#!/usr/bin/env bash
# Shared configuration for the GRS pipeline.
# Every value can be overridden from the environment (e.g. via submit.local.sh).
# Do NOT hardcode personal paths / allocation here - keep those in submit.local.sh.
set -euo pipefail

CONFIG_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$CONFIG_DIR}"

# --- inputs ---
# Root holding the genotype data and the SNP input. Override in submit.local.sh.
DATA_ROOT="${DATA_ROOT:-/path/to/your/grs_data}"
# SNP input: an .xlsx (e.g. the GWAS weight workbook) OR a .csv/.tsv candidate
# list. Must contain at least an rsID column and a risk/effect allele column;
# a weight column is only needed for weighted scoring.
SNP_INPUT="${SNP_INPUT:-$DATA_ROOT/41588_2018_321_MOESM3_ESM_suptable_28.xlsx}"
# Per-chromosome imputed dose VCFs: $IMPUTE_DIR/chr_<N>/chr<N>.dose.vcf.gz
IMPUTE_DIR="${IMPUTE_DIR:-$DATA_ROOT/imputation}"

# --- tools ---
PLINK2="${PLINK2:-plink2}"   # binary name on PATH, or an absolute path

# --- options ---
SCORE_MODE="${SCORE_MODE:-both}"   # both | weighted | unweighted
R2_THRESH="${R2_THRESH:-0}"        # 0 = keep all; e.g. 0.8 to filter at scoring time

# --- output ---
RUN_BASE="${RUN_BASE:-$PROJECT_ROOT/run_output}"
# Auto-prune SNP_extract_* run directories older than this many days on each
# new submit.sh invocation. 0 disables pruning.
RETENTION_DAYS="${RETENTION_DAYS:-30}"

if [[ "${DEBUG_CONFIG:-0}" == "1" ]]; then
  for v in PROJECT_ROOT DATA_ROOT SNP_INPUT IMPUTE_DIR PLINK2 SCORE_MODE R2_THRESH RUN_BASE RETENTION_DAYS; do
    echo "[CONFIG] $v=${!v}"
  done
fi
