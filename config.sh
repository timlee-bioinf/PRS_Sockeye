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

# --- multi-ancestry (optional; default off, classic single-list behavior unchanged) ---
# RUN_ANCESTRY=1 turns on: (a) per-ancestry SNP/weight lists instead of one SNP_INPUT,
# (b) ancestry inference (PCA projection + population classification) in step2b,
# (c) ancestry-routed scoring + discrete/continuous ancestry-normalized scores in step3.
# See scripts/setup_ancestry_reference.sh for the required one-time cache build.
RUN_ANCESTRY="${RUN_ANCESTRY:-0}"
# TSV, no header: <ancestry_label><TAB><snp_input path>. snp_input is the same
# xlsx/csv/tsv format prepare_copa_score_files.R already accepts, one file per
# ancestry. ancestry_label values MUST exactly match the population labels used
# in the reference panel's ANCESTRY_REF_LABEL_COL (e.g. EUR, EAS, SAS, AFR, AMR).
ANCESTRY_MAP="${ANCESTRY_MAP:-$DATA_ROOT/ancestry_input.tsv}"
# Reference ancestry panel (e.g. your "g2k" data): PLINK2 pgen/pvar/psam, path
# WITHOUT extension. Must carry rsIDs in its .pvar (matched by rsID, same as
# the rest of this pipeline - no build/liftover requirement as long as both
# sides use rsIDs).
ANCESTRY_REF_PFILE="${ANCESTRY_REF_PFILE:-$DATA_ROOT/ancestry_ref/ref}"
# Column in the reference .psam holding each sample's population label.
ANCESTRY_REF_LABEL_COL="${ANCESTRY_REF_LABEL_COL:-SuperPop}"
# Optional: sample-ID file (one per line) of related/duplicate reference
# samples to exclude before computing PCA + reference score distributions
# (e.g. a *.king.cutoff.out.id file). Leave empty to use all reference samples.
ANCESTRY_REF_EXCLUDE="${ANCESTRY_REF_EXCLUDE:-}"
# One-time setup output cache (see scripts/setup_ancestry_reference.sh). Not
# tied to RUN_ID - built once, reused by every pipeline run.
ANCESTRY_REF_CACHE="${ANCESTRY_REF_CACHE:-$DATA_ROOT/ancestry_ref_cache}"
ANCESTRY_N_PCS="${ANCESTRY_N_PCS:-10}"            # PCs computed for the reference PCA basis
ANCESTRY_N_POPCOMP="${ANCESTRY_N_POPCOMP:-5}"     # top PCs used for population classification (Mahalanobis)
ANCESTRY_N_NORM="${ANCESTRY_N_NORM:-4}"           # top PCs used for continuous PC-regression normalization (Z_norm1)
# Mahalanobis chi-sq p-value floor below which a call is flagged LowConfidence
# (samples are still assigned + scored - this never rejects, only flags).
ANCESTRY_PVAL_THRESH="${ANCESTRY_PVAL_THRESH:-1e-10}"

if [[ "${DEBUG_CONFIG:-0}" == "1" ]]; then
  for v in PROJECT_ROOT DATA_ROOT SNP_INPUT IMPUTE_DIR PLINK2 SCORE_MODE R2_THRESH RUN_BASE RETENTION_DAYS \
           RUN_ANCESTRY ANCESTRY_MAP ANCESTRY_REF_PFILE ANCESTRY_REF_LABEL_COL ANCESTRY_REF_EXCLUDE \
           ANCESTRY_REF_CACHE ANCESTRY_N_PCS ANCESTRY_N_POPCOMP ANCESTRY_N_NORM ANCESTRY_PVAL_THRESH; do
    echo "[CONFIG] $v=${!v}"
  done
fi
